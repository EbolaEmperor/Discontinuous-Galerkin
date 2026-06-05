#include "HybridHDG.h"
#include "Quadrature.h"

#include <vector>
#include <map>
#include <cmath>
#include <chrono>
#include <iostream>

#include <Eigen/SparseCholesky>
#ifdef EIGEN_USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::RowVectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::PartialPivLU;

namespace {

double cross2d(const Vector2d& a, const Vector2d& b) {
    return a(0) * b(1) - a(1) * b(0);
}

// Local edges of a triangle, following the repo convention used elsewhere
// (FEM::getConformingDOF, PostProcess): edge ee joins local vertices (a,b).
constexpr int kEdgeNodes[3][2] = {{1, 2}, {2, 0}, {0, 1}};

// ---- 1D Lagrange trace basis on the global edge parameter t in [0,1] -------
// Nodes run from the edge's min-index node (t=0) to its max-index node (t=1),
// so two elements sharing a face evaluate identical basis functions.
struct TraceBasis {
    int k;
    int n;                       // k + 1
    std::vector<double> nodes;   // n equispaced nodes on [0,1]

    explicit TraceBasis(int order) : k(order), n(order + 1) {
        nodes.resize(n);
        if (k == 0) {
            nodes[0] = 0.5;
        } else {
            for (int i = 0; i < n; ++i) nodes[i] = static_cast<double>(i) / k;
        }
    }

    // Lagrange basis values at parameter t -> length-n row.
    RowVectorXd eval(double t) const {
        RowVectorXd v = RowVectorXd::Ones(n);
        for (int i = 0; i < n; ++i) {
            double val = 1.0;
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                val *= (t - nodes[j]) / (nodes[i] - nodes[j]);
            }
            v(i) = val;
        }
        return v;
    }
};

// Everything needed for both the global assembly and the local reconstruction
// of one element.
struct LocalHDG {
    PartialPivLU<MatrixXd> lu;   // factorization of Aloc
    MatrixXd R;                  // (nq+nu) x nl
    VectorXd f0;                 // (nq+nu),  [0 ; (f,phi)]
    MatrixXd AK;                 // nl x nl,  element Schur contribution
    VectorXd bK;                 // nl,       element rhs contribution
    std::vector<int> gtrace;     // nl global trace dof indices
    int nq = 0, nu = 0, nl = 0;
};

// Assemble the local HDG operators for element t and statically condense.
LocalHDG computeLocal(int t, VecFEM& femQ, FEM& femU, const Mesh& mesh,
                      double tau, const ExactSolution& sol,
                      const std::map<std::pair<int, int>, int>& edgeMap,
                      const TraceBasis& trace) {
    LocalHDG L;
    const int nq = femQ.locDof;
    const int nu = femU.locDof;
    const int kt = trace.n;          // k+1 trace dofs per edge
    const int nl = 3 * kt;
    L.nq = nq; L.nu = nu; L.nl = nl;

    MatrixXd A_qq = MatrixXd::Zero(nq, nq);
    MatrixXd Dm   = MatrixXd::Zero(nq, nu);   // D[i,j] = (div vi, uj)
    MatrixXd Huu  = MatrixXd::Zero(nu, nu);
    MatrixXd E    = MatrixXd::Zero(nq, nl);   // <vi.n, mu_l>
    MatrixXd Ful  = MatrixXd::Zero(nu, nl);   // tau <ui, mu_l>
    MatrixXd Mll  = MatrixXd::Zero(nl, nl);   // tau <mu_l, mu_m>
    VectorXd fvec = VectorXd::Zero(nu);

    const double area = femU.area(t);
    Vector2d p1 = mesh.node.row(mesh.elem(t, 0)).transpose();
    Vector2d p2 = mesh.node.row(mesh.elem(t, 1)).transpose();
    Vector2d p3 = mesh.node.row(mesh.elem(t, 2)).transpose();

    // ---- volume integrals ----
    MatrixXd q2; VectorXd w2;
    Quadrature::quadpts2_my(2 * femQ.ord + 2, q2, w2);
    for (int q = 0; q < w2.size(); ++q) {
        Vector3d lam = q2.row(q).transpose();
        MatrixXd lamMat(1, 3); lamMat << lam.transpose();

        MatrixXd phiVec = femQ.computeBasisValue_all(lam);        // 2 x nq
        RowVectorXd div = femQ.computeBasisDiv_all(t, lam);       // 1 x nq
        RowVectorXd phiU = femU.computeBasisValue_all(lamMat).row(0); // 1 x nu

        double wgt = w2(q) * area;
        A_qq.noalias() += wgt * (phiVec.transpose() * phiVec);
        Dm.noalias()   += wgt * (div.transpose() * phiU);

        Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
        MatrixXd ptMat(1, 2); ptMat << pt.transpose();
        double fval = -sol.laplace_u_exact(ptMat)(0);            // f = -Lap u
        fvec.noalias() += (wgt * fval) * phiU.transpose();
    }

    // ---- edge integrals ----
    MatrixXd q1; VectorXd w1;
    Quadrature::quadpts1_my(2 * femQ.ord + 2, q1, w1);
    for (int ee = 0; ee < 3; ++ee) {
        int a = kEdgeNodes[ee][0];
        int b = kEdgeNodes[ee][1];
        int na = mesh.elem(t, a);
        int nb = mesh.elem(t, b);
        Vector2d Pa = mesh.node.row(na).transpose();
        Vector2d Pb = mesh.node.row(nb).transpose();
        Vector2d evec = Pb - Pa;
        double he = evec.norm();
        Vector2d nrm = Vector2d(evec(1), -evec(0)) / he;
        int wid = 3 - (a + b);
        Vector2d Pw = mesh.node.row(mesh.elem(t, wid)).transpose();
        if (cross2d(evec, Pw - Pa) < 0) nrm = -nrm;     // outward normal

        int ge = edgeMap.at({std::min(na, nb), std::max(na, nb)});
        bool maxIsB = (nb > na);
        int col0 = ee * kt;

        for (int q = 0; q < w1.size(); ++q) {
            double lamA = q1(q, 0);
            double lamB = q1(q, 1);
            double tg = maxIsB ? lamB : lamA;           // global edge parameter
            RowVectorXd mu = trace.eval(tg);            // 1 x (k+1)

            Vector3d lam = Vector3d::Zero();
            lam(a) = lamA; lam(b) = lamB;
            MatrixXd lamMat(1, 3); lamMat << lam.transpose();

            MatrixXd phiVec = femQ.computeBasisValue_all(lam);    // 2 x nq
            RowVectorXd nphi = nrm.transpose() * phiVec;          // 1 x nq  (v.n)
            RowVectorXd phiU = femU.computeBasisValue_all(lamMat).row(0); // 1 x nu

            double wgt = w1(q) * he;

            // E(:, edge block) += wgt * (v.n)^T * mu
            E.block(0, col0, nq, kt).noalias() += wgt * (nphi.transpose() * mu);
            // Ful(:, edge block) += wgt * tau * u^T * mu
            Ful.block(0, col0, nu, kt).noalias() += (wgt * tau) * (phiU.transpose() * mu);
            // Huu += wgt * tau * u^T u
            Huu.noalias() += (wgt * tau) * (phiU.transpose() * phiU);
            // Mll(edge block) += wgt * tau * mu^T mu
            Mll.block(col0, col0, kt, kt).noalias() += (wgt * tau) * (mu.transpose() * mu);
        }

        for (int l = 0; l < kt; ++l) L.gtrace.push_back(ge * kt + l);
    }

    // ---- assemble local saddle operator and statically condense ----
    const int m = nq + nu;
    MatrixXd Aloc(m, m);
    Aloc.topLeftCorner(nq, nq)     = A_qq;
    Aloc.topRightCorner(nq, nu)    = -Dm;
    Aloc.bottomLeftCorner(nu, nq)  = Dm.transpose();
    Aloc.bottomRightCorner(nu, nu) = Huu;

    L.R = MatrixXd::Zero(m, nl);
    L.R.topRows(nq)    = -E;
    L.R.bottomRows(nu) = Ful;

    MatrixXd S(nl, m);
    S.leftCols(nq)  = E.transpose();
    S.rightCols(nu) = Ful.transpose();

    L.f0 = VectorXd::Zero(m);
    L.f0.tail(nu) = fvec;

    L.lu.compute(Aloc);
    MatrixXd AinvR = L.lu.solve(L.R);          // m x nl
    VectorXd Ainvf = L.lu.solve(L.f0);         // m

    L.AK = Mll - S * AinvR;                    // SPD element contribution
    L.bK = S * Ainvf;
    return L;
}

// L2 projection of g = u_exact onto P_k of a single (boundary) edge.
VectorXd projectBoundaryTrace(int nmin, int nmax, const Mesh& mesh,
                              const ExactSolution& sol, const TraceBasis& trace) {
    int kt = trace.n;
    Vector2d Pmin = mesh.node.row(nmin).transpose();
    Vector2d Pmax = mesh.node.row(nmax).transpose();

    MatrixXd q1; VectorXd w1;
    Quadrature::quadpts1_my(2 * trace.k + 2, q1, w1);

    MatrixXd Me = MatrixXd::Zero(kt, kt);
    VectorXd rhs = VectorXd::Zero(kt);
    for (int q = 0; q < w1.size(); ++q) {
        double tg = q1(q, 1);                       // parameter along min->max
        RowVectorXd mu = trace.eval(tg);
        Vector2d pt = (1.0 - tg) * Pmin + tg * Pmax;
        MatrixXd ptMat(1, 2); ptMat << pt.transpose();
        double gval = sol.u_exact(ptMat)(0);
        double w = w1(q);                           // he cancels in M^{-1} rhs
        Me.noalias()  += w * (mu.transpose() * mu);
        rhs.noalias() += (w * gval) * mu.transpose();
    }
    return Me.ldlt().solve(rhs);
}

template <typename Solver>
VectorXd cholSolve(const SparseMatrix<double>& A, const VectorXd& b,
                   const char* name, std::string& outName) {
    Solver solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
        std::cerr << "[HDG] " << name << " factorization failed\n";
    VectorXd x = solver.solve(b);
    outName = name;
    return x;
}

} // namespace

HybridHDGResult solveHybridHDG(VecFEM& femQ, FEM& femU, Mesh& mesh,
                               const MatrixXi& elem2dofSigma,
                               const MatrixXi& elem2dofU,
                               const MatrixXi& edge,
                               const MatrixXi& edge2side,
                               double tau, const ExactSolution& sol,
                               int solver_type) {
    using clock = std::chrono::high_resolution_clock;

    const int NT = mesh.elem.rows();
    const int NE = edge.rows();
    const int kt = femQ.ord + 1;            // trace dofs per edge
    const int nLambda = NE * kt;

    TraceBasis trace(femQ.ord);

    // sorted node-pair -> global edge index
    std::map<std::pair<int, int>, int> edgeMap;
    for (int e = 0; e < NE; ++e)
        edgeMap[{edge(e, 0), edge(e, 1)}] = e;

    // ---- boundary trace dofs (Dirichlet) and free numbering ----
    std::vector<char> isDir(nLambda, 0);
    VectorXd lambdaVal = VectorXd::Zero(nLambda);
    for (int e = 0; e < NE; ++e) {
        bool boundary = (edge2side(e, 0) == -1) || (edge2side(e, 1) == -1);
        if (!boundary) continue;
        VectorXd le = projectBoundaryTrace(edge(e, 0), edge(e, 1), mesh, sol, trace);
        for (int l = 0; l < kt; ++l) {
            isDir[e * kt + l] = 1;
            lambdaVal(e * kt + l) = le(l);
        }
    }
    std::vector<int> redIdx(nLambda, -1);
    int nFree = 0;
    for (int g = 0; g < nLambda; ++g)
        if (!isDir[g]) redIdx[g] = nFree++;

    // ---- assemble SPD reduced system A_ff Lambda_f = b_f ----
    auto t0 = clock::now();
    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(NT) * 9 * kt * kt);
    VectorXd bred = VectorXd::Zero(nFree);

    for (int t = 0; t < NT; ++t) {
        LocalHDG L = computeLocal(t, femQ, femU, mesh, tau, sol, edgeMap, trace);
        for (int i = 0; i < L.nl; ++i) {
            int gi = L.gtrace[i];
            int ri = redIdx[gi];
            if (ri < 0) continue;                 // Dirichlet row: skip
            bred(ri) += L.bK(i);
            for (int j = 0; j < L.nl; ++j) {
                int gj = L.gtrace[j];
                double v = L.AK(i, j);
                int rj = redIdx[gj];
                if (rj >= 0) trip.emplace_back(ri, rj, v);
                else         bred(ri) -= v * lambdaVal(gj);   // eliminate known
            }
        }
    }
    SparseMatrix<double> A(nFree, nFree);
    A.setFromTriplets(trip.begin(), trip.end());
    auto t1 = clock::now();

    // ---- solve (default: CHOLMOD supernodal LLT) ----
    VectorXd lamFree;
    std::string solverName;
    if (solver_type == 0) {
        lamFree = cholSolve<Eigen::SimplicialLDLT<SparseMatrix<double>>>(A, bred, "SimplicialLDLT", solverName);
    } else if (solver_type == 1) {
        lamFree = cholSolve<Eigen::SimplicialLLT<SparseMatrix<double>>>(A, bred, "SimplicialLLT", solverName);
    } else if (solver_type == 2) {
        lamFree = cholSolve<Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper>>(A, bred, "ConjugateGradient", solverName);
    } else {
#ifdef EIGEN_USE_CHOLMOD
        lamFree = cholSolve<Eigen::CholmodSupernodalLLT<SparseMatrix<double>>>(A, bred, "CholmodSupernodalLLT", solverName);
#else
        lamFree = cholSolve<Eigen::SimplicialLLT<SparseMatrix<double>>>(A, bred, "SimplicialLLT (CHOLMOD unavailable)", solverName);
#endif
    }
    auto t2 = clock::now();

    for (int g = 0; g < nLambda; ++g)
        if (redIdx[g] >= 0) lambdaVal(g) = lamFree(redIdx[g]);

    // ---- local reconstruction of (q,u); sigma = grad u = -q ----
    auto t3 = clock::now();
    int nSig = elem2dofSigma.maxCoeff() + 1;
    int nU   = elem2dofU.maxCoeff() + 1;
    HybridHDGResult res;
    res.sigmah = VectorXd::Zero(nSig);
    res.uh     = VectorXd::Zero(nU);

    for (int t = 0; t < NT; ++t) {
        LocalHDG L = computeLocal(t, femQ, femU, mesh, tau, sol, edgeMap, trace);
        VectorXd lamK(L.nl);
        for (int i = 0; i < L.nl; ++i) lamK(i) = lambdaVal(L.gtrace[i]);

        VectorXd x = L.lu.solve(L.f0 + L.R * lamK);   // [Q ; U]
        VectorXd Q = x.head(L.nq);
        VectorXd U = x.tail(L.nu);

        for (int i = 0; i < L.nq; ++i) res.sigmah(elem2dofSigma(t, i)) = -Q(i); // sigma = -q
        for (int i = 0; i < L.nu; ++i) res.uh(elem2dofU(t, i)) = U(i);
    }

    auto t4 = clock::now();
    res.lambda      = lambdaVal;
    res.nLambda     = nLambda;
    res.nLambdaFree = nFree;
    res.tAssemble   = std::chrono::duration<double>(t1 - t0).count();
    res.tSolve      = std::chrono::duration<double>(t2 - t1).count();
    res.tRecover    = std::chrono::duration<double>(t4 - t3).count();
    res.solverName  = solverName;
    return res;
}
