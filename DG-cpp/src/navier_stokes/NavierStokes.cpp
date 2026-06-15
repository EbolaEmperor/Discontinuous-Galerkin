#include "NavierStokes.h"
#include "DGAssembly.h"     // assembleK_Poi2D, assembleIP_Poi2D (SIPG -Lap volume + interior)
#include "Quadrature.h"
#include "utils.h"          // numSplit3 (monomial multi-indices) for fast field rendering

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

namespace ns {

// ===========================================================================
// Small geometric / basis helpers (shared by the assemblers)
// ===========================================================================
namespace {

// Edge type 0=(0,1), 1=(1,2), 2=(2,0); dir 0=forward, 1=reverse, from the two
// local vertex indices of the edge inside an element (matches DGAssembly.cpp).
std::pair<int, int> edgeTypeDir(int a, int b) {
    if ((a == 0 && b == 1) || (a == 1 && b == 0)) return {0, (a > b) ? 1 : 0};
    if ((a == 1 && b == 2) || (a == 2 && b == 1)) return {1, (a > b) ? 1 : 0};
    return {2, (a == 0 && b == 2) ? 1 : 0};
}

// Per-(edge_type,dir,quad) reference-trace basis data, evaluated at the 1-D rule.
struct EdgeBasis {
    int nq = 0;
    std::vector<std::vector<std::vector<RowVectorXd>>> phi;     // [3][2][nq] : 1 x locDof
    std::vector<std::vector<std::vector<MatrixXd>>>    dphi;    // [3][2][nq] : 3 x locDof
    std::vector<std::vector<std::vector<MatrixXd>>>    Hlam;    // [3][2][nq] : 6 x locDof
};

EdgeBasis makeEdgeBasis(FEM& fem, const MatrixXd& quad1d) {
    EdgeBasis E;
    E.nq = static_cast<int>(quad1d.rows());
    E.phi.assign(3, std::vector<std::vector<RowVectorXd>>(2, std::vector<RowVectorXd>(E.nq)));
    E.dphi.assign(3, std::vector<std::vector<MatrixXd>>(2, std::vector<MatrixXd>(E.nq)));
    E.Hlam.assign(3, std::vector<std::vector<MatrixXd>>(2, std::vector<MatrixXd>(E.nq)));
    for (int et = 0; et < 3; ++et) {
        int i1 = 0, i2 = 1;
        if (et == 1) { i1 = 1; i2 = 2; }
        else if (et == 2) { i1 = 2; i2 = 0; }
        for (int dir = 0; dir < 2; ++dir)
            for (int q = 0; q < E.nq; ++q) {
                double l1 = quad1d(q, 0), l2 = quad1d(q, 1);
                Vector3d lam = Vector3d::Zero();
                if (dir == 0) { lam(i1) = l1; lam(i2) = l2; }
                else          { lam(i1) = l2; lam(i2) = l1; }
                E.phi[et][dir][q]  = fem.computeBasisValue_all(lam.transpose()).row(0);
                E.dphi[et][dir][q] = fem.computeBasisDlam_all(lam);
                E.Hlam[et][dir][q] = fem.computeBasisHlam_all(lam);
            }
    }
    return E;
}

// For element t and the global edge (n1,n2): local indices, (edge_type,dir),
// outward unit normal of t on that edge, edge length, and the physical points
// at the 1-D quad nodes (param running n1 -> n2).
struct EdgeOnElem {
    int et, dir;
    Vector2d nout;
    double he;
};
EdgeOnElem edgeOnElem(const Mesh& mesh, int t, int n1, int n2) {
    int a = -1, b = -1, w = -1;
    for (int k = 0; k < 3; ++k) {
        int v = mesh.elem(t, k);
        if (v == n1) a = k; else if (v == n2) b = k; else w = k;
    }
    Vector2d p1 = mesh.node.row(n1), p2 = mesh.node.row(n2), pw = mesh.node.row(mesh.elem(t, w));
    Vector2d evec = p2 - p1;
    double he = evec.norm();
    Vector2d nrm(evec.y() / he, -evec.x() / he);
    if (nrm.dot(pw - 0.5 * (p1 + p2)) > 0) nrm = -nrm;   // point away from interior
    auto td = edgeTypeDir(a, b);
    return {td.first, td.second, nrm, he};
}

void quadDeg(int deg, MatrixXd& quadL, VectorXd& w) {
    Quadrature::quadpts2_my(std::min(14, std::max(1, deg)), quadL, w);
}

} // namespace

// ===========================================================================
// Scalar DG mass matrix
// ===========================================================================
SparseMatrix<double> assembleScalarMassDG(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows(), locDof = fem.locDof;
    int nDof = elem2dof.maxCoeff() + 1;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);

    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(NT) * locDof * locDof);
    MatrixXd Me(locDof, locDof);
    for (int t = 0; t < NT; ++t) {
        Me.setZero();
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) Me.noalias() += (w(q) * area) * (phi[q].transpose() * phi[q]);
        for (int i = 0; i < locDof; ++i)
            for (int j = 0; j < locDof; ++j) trip.emplace_back(elem2dof(t, i), elem2dof(t, j), Me(i, j));
    }
    SparseMatrix<double> M(nDof, nDof);
    M.setFromTriplets(trip.begin(), trip.end());
    return M;
}

// ===========================================================================
// Boundary Nitsche (Dirichlet) contribution to the SIPG -Laplacian
// ===========================================================================
SparseMatrix<double> assembleNitscheDirichlet(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                              const MatrixXi& edge, const MatrixXi& edge2side,
                                              const VectorXi& isDir, double sigma, double beta) {
    int NE = edge.rows(), locDof = fem.locDof;
    int nDof = elem2dof.maxCoeff() + 1;
    MatrixXd quad1d; VectorXd w1d; fem.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem, quad1d);

    std::vector<Triplet<double>> trip;
    MatrixXd Mloc(locDof, locDof);
    for (int e = 0; e < NE; ++e) {
        if (!isDir(e)) continue;
        int t = (edge2side(e, 0) != -1) ? edge2side(e, 0) : edge2side(e, 1);
        int n1 = edge(e, 0), n2 = edge(e, 1);
        EdgeOnElem ei = edgeOnElem(mesh, t, n1, n2);
        double pen = sigma / ei.he;
        Mloc.setZero();
        for (int q = 0; q < E.nq; ++q) {
            const RowVectorXd& phi = E.phi[ei.et][ei.dir][q];      // 1 x locDof
            MatrixXd g = fem.Dlam[t] * E.dphi[ei.et][ei.dir][q];   // 2 x locDof
            RowVectorXd dn = ei.nout.transpose() * g;              // grad.n , 1 x locDof
            double whe = w1d(q) * ei.he;
            // -(grad u.n) v - beta (grad v.n) u + (sigma/h) u v   (i=test, j=trial)
            Mloc.noalias() += whe * (-phi.transpose() * dn
                                     - beta * dn.transpose() * phi
                                     + pen * phi.transpose() * phi);
        }
        for (int i = 0; i < locDof; ++i)
            for (int j = 0; j < locDof; ++j) trip.emplace_back(elem2dof(t, i), elem2dof(t, j), Mloc(i, j));
    }
    SparseMatrix<double> B(nDof, nDof);
    B.setFromTriplets(trip.begin(), trip.end());
    return B;
}

// ===========================================================================
// DG weak first-derivative operators Gx, Gy  (central interior flux, single-
// sided boundary trace).  (Gm p)_i = \int (d_m p) phi_i ;  Gm * const == 0.
// ===========================================================================
void assembleWeakGrad(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                      const MatrixXi& edge, const MatrixXi& edge2side,
                      SparseMatrix<double>& Gx, SparseMatrix<double>& Gy) {
    int NT = mesh.elem.rows(), NE = edge.rows(), locDof = fem.locDof;
    int nDof = elem2dof.maxCoeff() + 1;

    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phiV(nq);
    std::vector<MatrixXd> dphiV(nq);
    for (int q = 0; q < nq; ++q) {
        phiV[q]  = fem.computeBasisValue_all(quadL.row(q)).row(0);
        dphiV[q] = fem.computeBasisDlam_all(quadL.row(q));
    }
    std::vector<Triplet<double>> tx, ty;

    // --- volume:  -\int_K phi_j (d_m phi_i) ---
    for (int t = 0; t < NT; ++t) {
        double area = fem.area(t);
        MatrixXd Vx = MatrixXd::Zero(locDof, locDof), Vy = MatrixXd::Zero(locDof, locDof);
        for (int q = 0; q < nq; ++q) {
            MatrixXd gp = fem.Dlam[t] * dphiV[q];           // 2 x locDof : grad phi_i
            double wa = w(q) * area;
            Vx.noalias() -= wa * gp.row(0).transpose() * phiV[q];   // (i,j) -= wa*dphi_x_i*phi_j
            Vy.noalias() -= wa * gp.row(1).transpose() * phiV[q];
        }
        for (int i = 0; i < locDof; ++i)
            for (int j = 0; j < locDof; ++j) {
                tx.emplace_back(elem2dof(t, i), elem2dof(t, j), Vx(i, j));
                ty.emplace_back(elem2dof(t, i), elem2dof(t, j), Vy(i, j));
            }
    }

    // --- faces ---
    MatrixXd quad1d; VectorXd w1d; fem.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem, quad1d);
    auto addBlock = [&](std::vector<Triplet<double>>& T, const MatrixXi& e2d,
                        int tr, int tc, const MatrixXd& B) {
        for (int i = 0; i < locDof; ++i)
            for (int j = 0; j < locDof; ++j) T.emplace_back(e2d(tr, i), e2d(tc, j), B(i, j));
    };
    for (int e = 0; e < NE; ++e) {
        int t1 = edge2side(e, 0), t2 = edge2side(e, 1);
        int n1 = edge(e, 0), n2 = edge(e, 1);
        if (t1 != -1 && t2 != -1) {
            EdgeOnElem e1 = edgeOnElem(mesh, t1, n1, n2);   // normal outward from t1
            EdgeOnElem e2 = edgeOnElem(mesh, t2, n1, n2);
            Vector2d n = e1.nout;                            // from t1 -> t2
            double he = e1.he;
            MatrixXd B11x = MatrixXd::Zero(locDof, locDof), B12x = B11x, B21x = B11x, B22x = B11x;
            MatrixXd B11y = B11x, B12y = B11x, B21y = B11x, B22y = B11x;
            for (int q = 0; q < E.nq; ++q) {
                const RowVectorXd& p1 = E.phi[e1.et][e1.dir][q];
                const RowVectorXd& p2 = E.phi[e2.et][e2.dir][q];
                double whe = w1d(q) * he;
                double hx = 0.5 * whe * n.x(), hy = 0.5 * whe * n.y();
                // block (test K, trial side): {phi_j} n_m^K phi_i^K
                B11x.noalias() += hx * p1.transpose() * p1;  B11y.noalias() += hy * p1.transpose() * p1;
                B12x.noalias() += hx * p1.transpose() * p2;  B12y.noalias() += hy * p1.transpose() * p2;
                B21x.noalias() -= hx * p2.transpose() * p1;  B21y.noalias() -= hy * p2.transpose() * p1;
                B22x.noalias() -= hx * p2.transpose() * p2;  B22y.noalias() -= hy * p2.transpose() * p2;
            }
            addBlock(tx, elem2dof, t1, t1, B11x); addBlock(tx, elem2dof, t1, t2, B12x);
            addBlock(tx, elem2dof, t2, t1, B21x); addBlock(tx, elem2dof, t2, t2, B22x);
            addBlock(ty, elem2dof, t1, t1, B11y); addBlock(ty, elem2dof, t1, t2, B12y);
            addBlock(ty, elem2dof, t2, t1, B21y); addBlock(ty, elem2dof, t2, t2, B22y);
        } else {
            int t = (t1 != -1) ? t1 : t2;
            EdgeOnElem ei = edgeOnElem(mesh, t, n1, n2);
            MatrixXd Bx = MatrixXd::Zero(locDof, locDof), By = Bx;
            for (int q = 0; q < E.nq; ++q) {
                const RowVectorXd& phi = E.phi[ei.et][ei.dir][q];   // single-sided trace
                double whe = w1d(q) * ei.he;
                Bx.noalias() += (whe * ei.nout.x()) * phi.transpose() * phi;
                By.noalias() += (whe * ei.nout.y()) * phi.transpose() * phi;
            }
            addBlock(tx, elem2dof, t, t, Bx);
            addBlock(ty, elem2dof, t, t, By);
        }
    }

    Gx.resize(nDof, nDof); Gx.setFromTriplets(tx.begin(), tx.end());
    Gy.resize(nDof, nDof); Gy.setFromTriplets(ty.begin(), ty.end());
}

// Optional coupled grad-div stabilisation:
//   gd * int_K (du_x/dx + du_y/dy) (dw_x/dx + dw_y/dy).
// The default fast path keeps gd=0 and instead damps divergence in the PPE.
static SparseMatrix<double> assembleVolumeGradDivDG(FEM& fem, const Mesh& mesh,
                                                     const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows(), locDof = fem.locDof;
    int nDof = elem2dof.maxCoeff() + 1;

    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<MatrixXd> dphiV(nq);
    for (int q = 0; q < nq; ++q) dphiV[q] = fem.computeBasisDlam_all(quadL.row(q));

    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(NT) * 4 * locDof * locDof);
    MatrixXd Kxx(locDof, locDof), Kxy(locDof, locDof), Kyx(locDof, locDof), Kyy(locDof, locDof);
    for (int t = 0; t < NT; ++t) {
        Kxx.setZero(); Kxy.setZero(); Kyx.setZero(); Kyy.setZero();
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) {
            MatrixXd gp = fem.Dlam[t] * dphiV[q];
            double wa = w(q) * area;
            const RowVectorXd gx = gp.row(0);
            const RowVectorXd gy = gp.row(1);
            Kxx.noalias() += wa * gx.transpose() * gx;
            Kxy.noalias() += wa * gx.transpose() * gy;
            Kyx.noalias() += wa * gy.transpose() * gx;
            Kyy.noalias() += wa * gy.transpose() * gy;
        }
        for (int i = 0; i < locDof; ++i) {
            int ri = elem2dof(t, i);
            for (int j = 0; j < locDof; ++j) {
                int cj = elem2dof(t, j);
                trip.emplace_back(ri,        cj,        Kxx(i, j));
                trip.emplace_back(ri,        nDof + cj, Kxy(i, j));
                trip.emplace_back(nDof + ri, cj,        Kyx(i, j));
                trip.emplace_back(nDof + ri, nDof + cj, Kyy(i, j));
            }
        }
    }
    SparseMatrix<double> GD(2 * nDof, 2 * nDof);
    GD.setFromTriplets(trip.begin(), trip.end());
    return GD;
}

static SparseMatrix<double> coupledVelocityMatrix(const SparseMatrix<double>& Hu,
                                                   const SparseMatrix<double>& Hv,
                                                   const SparseMatrix<double>& GD,
                                                   double gradDiv) {
    const int n = Hu.rows();
    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(Hu.nonZeros() + Hv.nonZeros() + GD.nonZeros()));
    for (int k = 0; k < Hu.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(Hu, k); it; ++it)
            trip.emplace_back(it.row(), it.col(), it.value());
    for (int k = 0; k < Hv.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(Hv, k); it; ++it)
            trip.emplace_back(n + it.row(), n + it.col(), it.value());
    if (gradDiv != 0.0)
        for (int k = 0; k < GD.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(GD, k); it; ++it)
                trip.emplace_back(it.row(), it.col(), gradDiv * it.value());
    SparseMatrix<double> H(2 * n, 2 * n);
    H.setFromTriplets(trip.begin(), trip.end());
    return H;
}

// ===========================================================================
// NSIntegrator
// ===========================================================================
NSIntegrator::NSIntegrator(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
                           const MatrixXi& edge, const MatrixXi& edge2side,
                           const BCData& bc, double nu, double dt, double sigma, double beta,
                           double gradDiv, int pressureMode, double ppeDivDamping)
    : fem_(fem), mesh_(mesh), elem2dof_(elem2dof), edge_(edge), edge2side_(edge2side),
      bc_(bc), nu_(nu), dt_(dt), sigma_(sigma), beta_(beta),
      gradDiv_(gradDiv), ppeDivDamping_(ppeDivDamping), pressureMode_(pressureMode),
      nDof_(elem2dof.maxCoeff() + 1), n_(0) {
    buildOperators();
}

void NSIntegrator::buildOperators() {
    M_ = assembleScalarMassDG(fem_, mesh_, elem2dof_);
    SparseMatrix<double> K  = assembleK_Poi2D(fem_, mesh_, elem2dof_);
    SparseMatrix<double> IP = assembleIP_Poi2D(fem_, mesh_, elem2dof_, edge_, edge2side_, sigma_, beta_);
    SparseMatrix<double> base = K + IP;   // SIPG -Lap, interior only (natural Neumann)

    auto dirMask = [&](const VectorXi& bccomp) {
        VectorXi m = VectorXi::Zero(edge_.rows());
        for (int e = 0; e < edge_.rows(); ++e) if (bccomp(e) == BCN_DIRICHLET) m(e) = 1;
        return m;
    };
    Au_ = base + assembleNitscheDirichlet(fem_, mesh_, elem2dof_, edge_, edge2side_, dirMask(bc_.bcU), sigma_, beta_);
    Av_ = base + assembleNitscheDirichlet(fem_, mesh_, elem2dof_, edge_, edge2side_, dirMask(bc_.bcV), sigma_, beta_);
    Ap_ = base + assembleNitscheDirichlet(fem_, mesh_, elem2dof_, edge_, edge2side_, dirMask(bc_.bcP), sigma_, beta_);
    bool hasPressureDirichlet = false;
    for (int e = 0; e < bc_.bcP.size(); ++e)
        if (bc_.bcP(e) == BCN_DIRICHLET) { hasPressureDirichlet = true; break; }
    if (!hasPressureDirichlet) {
        // Pure Neumann PPE determines pressure up to a constant.  A tiny mass
        // gauge fixes the nullspace without changing the pressure gradient at
        // truncation-error scale.
        Ap_ = Ap_ + 1e-12 * M_;
    }

    assembleWeakGrad(fem_, mesh_, elem2dof_, edge_, edge2side_, Gx_, Gy_);

    Hu_  = (1.5 / dt_) * M_ + nu_ * Au_;     // BDF2 viscous (gamma0 = 3/2)
    Hv_  = (1.5 / dt_) * M_ + nu_ * Av_;
    Hu1_ = (1.0 / dt_) * M_ + nu_ * Au_;     // BDF1 bootstrap (gamma0 = 1)
    Hv1_ = (1.0 / dt_) * M_ + nu_ * Av_;
    SparseMatrix<double> GD;
    if (gradDiv_ > 0.0) {
        GD = assembleVolumeGradDivDG(fem_, mesh_, elem2dof_);
        Huv_  = coupledVelocityMatrix(Hu_,  Hv_,  GD, gradDiv_);
        Huv1_ = coupledVelocityMatrix(Hu1_, Hv1_, GD, gradDiv_);
    }

    luM_.compute(M_);
    luAp_.compute(Ap_);
    if (gradDiv_ > 0.0) {
        luHuv_.compute(Huv_);
        luHuv1_.compute(Huv1_);
    } else {
        luHu_.compute(Hu_);   luHv_.compute(Hv_);
        luHu1_.compute(Hu1_); luHv1_.compute(Hv1_);
    }
    if (luM_.info() != Success || luAp_.info() != Success ||
        (gradDiv_ > 0.0 && (luHuv_.info() != Success || luHuv1_.info() != Success)) ||
        (gradDiv_ <= 0.0 && (luHu_.info() != Success || luHv_.info() != Success ||
                             luHu1_.info() != Success || luHv1_.info() != Success)))
        std::cerr << "NSIntegrator: a factorisation failed (try a larger penalty sigma)\n";
}

void NSIntegrator::setDt(double newDt) {
    dt_ = newDt;
    Hu_  = (1.5 / dt_) * M_ + nu_ * Au_;
    Hv_  = (1.5 / dt_) * M_ + nu_ * Av_;
    Hu1_ = (1.0 / dt_) * M_ + nu_ * Au_;
    Hv1_ = (1.0 / dt_) * M_ + nu_ * Av_;
    if (gradDiv_ > 0.0) {
        SparseMatrix<double> GD = assembleVolumeGradDivDG(fem_, mesh_, elem2dof_);
        Huv_  = coupledVelocityMatrix(Hu_,  Hv_,  GD, gradDiv_);
        Huv1_ = coupledVelocityMatrix(Hu1_, Hv1_, GD, gradDiv_);
        luHuv_.compute(Huv_);
        luHuv1_.compute(Huv1_);
    } else {
        luHu_.compute(Hu_);   luHv_.compute(Hv_);
        luHu1_.compute(Hu1_); luHv1_.compute(Hv1_);
    }
    // Reset to BDF1 so the next step doesn't use stale history at the old dt.
    n_ = 0;
    uPrev_ = u_; vPrev_ = v_;
    cxPrev_ = VectorXd::Zero(nDof_); cyPrev_ = VectorXd::Zero(nDof_);
}

void NSIntegrator::setInitial(const VectorXd& u0, const VectorXd& v0, const VectorXd& p0) {
    u_ = u0; v_ = v0;
    uPrev_ = u0; vPrev_ = v0;
    p_ = (p0.size() == nDof_) ? p0 : VectorXd::Zero(nDof_);
    cxPrev_ = VectorXd::Zero(nDof_); cyPrev_ = VectorXd::Zero(nDof_);
    n_ = 0;
}

VectorXd NSIntegrator::massSolve(const VectorXd& b) const { return luM_.solve(b); }

bool NSIntegrator::step(double tEnd) {
    return stepWithBodyForce(tEnd, VectorXd(), VectorXd());
}

bool NSIntegrator::stepWithBodyForce(double tEnd,
                                     const VectorXd& loadFx, const VectorXd& loadFy) {
    if (filterStrength > 0.0 && !filterBuilt_) buildModalFilter_();
    if (avBeta > 0.0 && !avBuilt_) buildAV_();
    if (hvFac > 0.0 && !hvBuilt_) buildHV_();
    const double tN = tEnd - dt_;       // t^n (explicit terms live here)
    const bool bdf2 = (n_ >= 1);
    const double g0 = bdf2 ? 1.5 : 1.0;

    // --- explicit convection load at t^n ---
    VectorXd cx, cy;
    if (includeConvection) assembleConvection(u_, v_, tN, cx, cy);
    else { cx = VectorXd::Zero(nDof_); cy = VectorXd::Zero(nDof_); }

    // history mass term  ha = (BDF2) 2 u^n - 1/2 u^{n-1}  | (BDF1) u^n
    VectorXd haU = bdf2 ? (2.0 * u_ - 0.5 * uPrev_) : u_;
    VectorXd haV = bdf2 ? (2.0 * v_ - 0.5 * vPrev_) : v_;
    // extrapolated convection load  c* = (BDF2) 2 c^n - c^{n-1} | (BDF1) c^n
    VectorXd cxs = bdf2 ? (2.0 * cx - cxPrev_) : cx;
    VectorXd cys = bdf2 ? (2.0 * cy - cyPrev_) : cy;

    // intermediate field  uhat = ha - dt * N*   (N* = M^{-1} c*)
    VectorXd nxStar = massSolve(cxs);
    VectorXd nyStar = massSolve(cys);
    VectorXd uhat = haU - dt_ * nxStar;
    VectorXd vhat = haV - dt_ * nyStar;

    // --- pressure Poisson ---
    VectorXd bp = assemblePressureNeumann(u_, v_, tEnd) + presDirichletLift(tEnd);
    if (pressureMode_ == NSPRESSURE_DIRECT_PPE) {
        // Consistent PPE:  Delta p = -div N*, so the SIPG -Lap RHS is div N*.
        // This removes the amplified discrete-divergence history term present in
        // the projection RHS and gives a cleaner pressure gradient at the boundary.
        bp += Gx_ * nxStar + Gy_ * nyStar;
        if (ppeDivDamping_ > 0.0) {
            VectorXd uEx = bdf2 ? (2.0 * u_ - uPrev_) : u_;
            VectorXd vEx = bdf2 ? (2.0 * v_ - vPrev_) : v_;
            bp += -ppeDivDamping_ * (Gx_ * uEx + Gy_ * vEx);
        }
    } else {
        // Projection PPE: Delta p = (1/dt) div(ha - dt N*).
        bp += -(1.0 / dt_) * (Gx_ * uhat + Gy_ * vhat);
    }
    p_ = luAp_.solve(bp);
    if (luAp_.info() != Success) return false;

    // --- viscous Helmholtz (per component) ---
    //   H u^{n+1} = (1/dt) M ha - c* - G p + nu * DirichletLift  + IB body force
    VectorXd rhsU = (1.0 / dt_) * (M_ * haU) - cxs - Gx_ * p_ + nu_ * velDirichletLift(0, tEnd);
    VectorXd rhsV = (1.0 / dt_) * (M_ * haV) - cys - Gy_ * p_ + nu_ * velDirichletLift(1, tEnd);
    if (loadFx.size() == nDof_) rhsU += loadFx;
    if (loadFy.size() == nDof_) rhsV += loadFy;
    VectorXd uNew, vNew;
    if (gradDiv_ > 0.0) {
        VectorXd rhs(2 * nDof_);
        rhs.head(nDof_) = rhsU;
        rhs.tail(nDof_) = rhsV;
        VectorXd uv = bdf2 ? luHuv_.solve(rhs) : luHuv1_.solve(rhs);
        if ((bdf2 ? luHuv_.info() : luHuv1_.info()) != Success) return false;
        uNew = uv.head(nDof_);
        vNew = uv.tail(nDof_);
    } else {
        uNew = bdf2 ? luHu_.solve(rhsU) : luHu1_.solve(rhsU);
        vNew = bdf2 ? luHv_.solve(rhsV) : luHv1_.solve(rhsV);
        if ((bdf2 ? luHu_.info() : luHu1_.info()) != Success ||
            (bdf2 ? luHv_.info() : luHv1_.info()) != Success) return false;
    }
    (void)g0;

    // --- semi-implicit immersed-boundary no-slip constraint (Schur correction) ---
    // Impose u_h(X_k) = V_k at the armed markers as a Lagrange-multiplier
    // correction to the viscous velocity, solved by CG reusing luHu_/luHv_.
    // Unconditionally stable in the constraint strength (no explicit forcing).
    if (ib_.active && gradDiv_ <= 0.0) {
        int itu = 0, itv = 0;
        applyIBConstraint_(uNew, ib_.V.col(0), bdf2 ? luHu_ : luHu1_, itu);
        applyIBConstraint_(vNew, ib_.V.col(1), bdf2 ? luHv_ : luHv1_, itv);
        ib_.lastIters = std::max(itu, itv);
        // EXPERIMENTAL constraint<->projection sub-iterations: the IB correction
        // above injects a local divergence error (its force never entered the
        // PPE).  Project the corrected velocity back to ~div-free, which
        // slightly violates the marker no-slip, and re-apply the constraint;
        // 1-3 rounds approximate the coupled (pressure, IB-force) saddle point.
        for (int it = 0; it < ibSubIters; ++it) {
            // projection:  Ap q = -(g0/dt) (Gx u + Gy v)  =>  Lap q = (g0/dt) div u
            VectorXd bq = -(g0 / dt_) * (Gx_ * uNew + Gy_ * vNew);
            VectorXd q = luAp_.solve(bq);
            uNew.noalias() -= (dt_ / g0) * massSolve(Gx_ * q);
            vNew.noalias() -= (dt_ / g0) * massSolve(Gy_ * q);
            p_.noalias() += q;
            if (it == ibSubIters - 1 && ibEndWithProjection) break;
            applyIBConstraint_(uNew, ib_.V.col(0), bdf2 ? luHu_ : luHu1_, itu);
            applyIBConstraint_(vNew, ib_.V.col(1), bdf2 ? luHv_ : luHv1_, itv);
            ib_.lastIters = std::max({ib_.lastIters, itu, itv});
        }
        ib_.lastResid = std::max((ibInterp_(uNew) - ib_.V.col(0)).cwiseAbs().maxCoeff(),
                                 (ibInterp_(vNew) - ib_.V.col(1)).cwiseAbs().maxCoeff());
    }

    // --- per-element modal filter (grid-scale / checkerboard dissipation) ---
    // The only fluid-side damping of grid-scale DG modes (no slope limiter).
    // Sensor-gated so smooth elements keep full P_k accuracy; troubled (oscillatory)
    // elements have their top-degree modal content attenuated, which removes shed
    // checkerboard packets before they persist/amplify in the wake.
    if (filterStrength > 0.0) {
        applyModalFilter_(uNew);
        applyModalFilter_(vNew);
    }
    // Sensor-gated artificial viscosity: damp INTER-element checkerboard (the
    // element-to-element sign alternation the modal filter cannot reach), which is
    // the only fluid-side dissipation of shed grid-scale wake oscillations.
    if (avBeta > 0.0) applyAV_(uNew, vNew);
    // Global biharmonic hyperviscosity: scale-selective grid-scale dissipation
    // (k^4) -- removes amplifying wake checkerboard without touching resolved scales.
    if (hvFac > 0.0) applyHV_(uNew, vNew);

    uPrev_ = u_; vPrev_ = v_;
    u_ = uNew;   v_ = vNew;
    cxPrev_ = cx; cyPrev_ = cy;
    ++n_;
    return true;
}

// ---------------------------------------------------------------------------
// Per-element modal filter (spectral vanishing-viscosity-style, sensor gated).
// ---------------------------------------------------------------------------
void NSIntegrator::buildModalFilter_() {
    if (filterBuilt_) return;
    const int nd = fem_.locDof, ord = fem_.ord;
    // Graded monomials in affine coords (b,c) = (lam2, lam3): all b^i c^j with
    // i+j <= ord, ordered by total degree d=i+j (then by i).  Count == nd.
    std::vector<std::pair<int,int>> exps;
    for (int d = 0; d <= ord; ++d)
        for (int i = d; i >= 0; --i) exps.emplace_back(i, d - i);
    filtDeg_.resize(nd);
    filtMaxDeg_ = ord;
    for (int k = 0; k < nd; ++k) filtDeg_[k] = exps[k].first + exps[k].second;

    MatrixXd quadL; VectorXd qw; fem_.quad2d(quadL, qw);
    const int nq = (int)qw.size();
    auto monoAt = [&](double b, double c, int k) {
        return std::pow(b, exps[k].first) * std::pow(c, exps[k].second);
    };
    // Mmono (nd x nq) monomials at quad pts; Phi (nq x nd) basis at quad pts.
    MatrixXd Mm(nd, nq);
    for (int q = 0; q < nq; ++q) {
        double b = quadL(q, 1), c = quadL(q, 2);
        for (int k = 0; k < nd; ++k) Mm(k, q) = monoAt(b, c, k);
    }
    MatrixXd Phi = fem_.computeBasisValue_all(quadL);   // nq x nd
    MatrixXd Wq = qw.asDiagonal();
    MatrixXd B = Mm * Wq * Mm.transpose();              // nd x nd : <m_j,m_l>
    MatrixXd C = Mm * Wq * Phi;                          // nd x nd : <m_j,phi_i>
    // Orthonormal graded modal basis psi_k = sum_j T_kj m_j with T B T^T = I,
    // T lower-triangular (graded) -> T = R^{-1}, B = R R^T (Cholesky).
    Eigen::LLT<MatrixXd> llt(B);
    MatrixXd R = llt.matrixL();
    MatrixXd T = R.inverse();
    MatrixXd nodes = fem_.lagrangeNodes();              // nd x 3
    MatrixXd G(nd, nd);                                  // G(i,k) = m_k(node_i)
    for (int i = 0; i < nd; ++i) {
        double b = nodes(i, 1), c = nodes(i, 2);
        for (int k = 0; k < nd; ++k) G(i, k) = monoAt(b, c, k);
    }
    filtV_ = G * T.transpose();    // V(i,k) = psi_k(node_i)         (modal -> nodal)
    filtN_ = T * C;                // N(k,i) = <psi_k, phi_i>        (nodal -> modal)
    filterBuilt_ = true;
}

void NSIntegrator::applyModalFilter_(VectorXd& w) const {
    // buildModalFilter_ is const-incompatible (mutates members) so it is invoked
    // from step(); guarantee it ran.
    const int nd = fem_.locDof, NT = mesh_.elem.rows();
    if (!filterBuilt_ || filtMaxDeg_ <= 0) return;
    const double lo = filterSensorLo, hi = std::max(filterSensorHi, filterSensorLo + 1e-9);
    VectorXd fe(nd), c(nd);
    for (int e = 0; e < NT; ++e) {
        for (int i = 0; i < nd; ++i) fe(i) = w(elem2dof_(e, i));
        c.noalias() = filtN_ * fe;
        double tot = 0.0, top = 0.0;
        for (int k = 0; k < nd; ++k) {
            double e2 = c(k) * c(k);
            tot += e2;
            if (filtDeg_[k] == filtMaxDeg_) top += e2;
        }
        if (tot < 1e-300) continue;
        double frac = top / tot;
        if (frac <= lo) continue;                       // smooth element -> untouched
        double g = (frac - lo) / (hi - lo);
        g = g < 0 ? 0 : (g > 1 ? 1 : g);
        g = g * g * (3.0 - 2.0 * g);                    // smoothstep
        double damp = filterStrength * g;
        if (damp < 1e-12) continue;
        const double scale = 1.0 - damp;
        for (int k = 0; k < nd; ++k) if (filtDeg_[k] == filtMaxDeg_) c(k) *= scale;
        fe.noalias() = filtV_ * c;
        for (int i = 0; i < nd; ++i) w(elem2dof_(e, i)) = fe(i);
    }
}

// ---------------------------------------------------------------------------
// Sensor-gated artificial viscosity for inter-element checkerboard modes.
// ---------------------------------------------------------------------------
void NSIntegrator::buildAV_() {
    if (avBuilt_) return;
    const int NT = mesh_.elem.rows();
    // Face-neighbour lists from edge2side (interior edges only).
    faceNbr_.assign(NT, {-1, -1, -1});
    std::vector<int> cnt(NT, 0);
    for (int e = 0; e < edge_.rows(); ++e) {
        int t1 = edge2side_(e, 0), t2 = edge2side_(e, 1);
        if (t1 >= 0 && t2 >= 0) {
            if (cnt[t1] < 3) faceNbr_[t1][cnt[t1]++] = t2;
            if (cnt[t2] < 3) faceNbr_[t2][cnt[t2]++] = t1;
        }
    }
    // Representative h^2 from the median element area (area ~ 0.43 h^2 -> h^2 ~ 2.3 area).
    std::vector<double> ar(NT);
    for (int e = 0; e < NT; ++e) ar[e] = fem_.area(e);
    std::nth_element(ar.begin(), ar.begin() + NT / 2, ar.end());
    avH2_ = 2.31 * ar[NT / 2];
    const double beta = avBeta * avH2_;
    // (M + beta * A) for each velocity component (BCs of Au_/Av_), factorised once.
    SparseMatrix<double> Su = M_ + beta * Au_;
    SparseMatrix<double> Sv = M_ + beta * Av_;
    luAVu_.compute(Su);
    luAVv_.compute(Sv);
    avBuilt_ = true;
}

void NSIntegrator::applyAV_(VectorXd& u, VectorXd& v) const {
    if (!avBuilt_) return;
    const int NT = mesh_.elem.rows(), nd = fem_.locDof;
    // One implicit diffusion step (unconditionally stable).  u_diff smooths grid-
    // scale content of ANY kind (inter-element checkerboard, intra-element high
    // modes, gradient-alternation) while leaving smooth modes ~unchanged.
    VectorXd ud = luAVu_.solve(M_ * u);
    VectorXd vd = luAVv_.solve(M_ * v);
    // Self-sensor: the per-element diffusion residual ||u - u_diff|| normalised by
    // the global velocity scale.  ~0 where the field is smooth (so smooth elements
    // are NEVER blended -> no accuracy loss, no compounding over a long run), large
    // exactly where there is grid-scale content the diffusion removes -- which the
    // old element-MEAN Laplacian sensor missed for intra-element / gradient-
    // alternation modes (e.g. the blade-tip withdrawal checkerboard).
    double vrms = 0.0;
    for (int i = 0; i < nDof_; ++i) vrms += u(i) * u(i) + v(i) * v(i);
    vrms = std::sqrt(vrms / std::max(1, nDof_)) + 1e-12;
    const double lo = avSensorLo, hi = std::max(avSensorHi, avSensorLo + 1e-9);
    for (int e = 0; e < NT; ++e) {
        double r2 = 0.0;
        for (int i = 0; i < nd; ++i) {
            int d = elem2dof_(e, i);
            double dxu = u(d) - ud(d), dxv = v(d) - vd(d);
            r2 += dxu * dxu + dxv * dxv;
        }
        double resid = std::sqrt(r2 / nd) / vrms;       // diffusion-residual sensor
        double g = 0.0;
        if (resid > lo) {
            g = (resid - lo) / (hi - lo);
            g = g < 0 ? 0 : (g > 1 ? 1 : g);
            g = g * g * (3.0 - 2.0 * g);                // smoothstep
        }
        if (avGlobalFloor > g) g = avGlobalFloor;       // ungated floor (free-decay)
        if (g < 1e-9) continue;
        for (int i = 0; i < nd; ++i) {
            int d = elem2dof_(e, i);
            u(d) = (1.0 - g) * u(d) + g * ud(d);
            v(d) = (1.0 - g) * v(d) + g * vd(d);
        }
    }
}

// ---------------------------------------------------------------------------
// Global implicit hyperviscosity (biharmonic grid-scale filter).
// ---------------------------------------------------------------------------
void NSIntegrator::buildHV_() {
    if (hvBuilt_) return;
    const int NT = mesh_.elem.rows(), nd = fem_.locDof;
    // Block-diagonal inverse mass M^{-1} (DG mass is block-diagonal per element).
    std::vector<Eigen::Triplet<double>> tr;
    tr.reserve((size_t)NT * nd * nd);
    MatrixXd Mb(nd, nd);
    for (int e = 0; e < NT; ++e) {
        for (int i = 0; i < nd; ++i)
            for (int j = 0; j < nd; ++j)
                Mb(i, j) = M_.coeff(elem2dof_(e, i), elem2dof_(e, j));
        MatrixXd Mi = Mb.inverse();
        for (int i = 0; i < nd; ++i)
            for (int j = 0; j < nd; ++j)
                tr.emplace_back(elem2dof_(e, i), elem2dof_(e, j), Mi(i, j));
    }
    SparseMatrix<double> Minv(nDof_, nDof_);
    Minv.setFromTriplets(tr.begin(), tr.end());
    // A M^{-1} A ~ M*nabla^4 (SPD).  Representative h^4 from the median element area.
    std::vector<double> ar(NT);
    for (int e = 0; e < NT; ++e) ar[e] = fem_.area(e);
    std::nth_element(ar.begin(), ar.begin() + NT / 2, ar.end());
    double h2 = 2.31 * ar[NT / 2];
    double beta = hvFac * h2 * h2;
    SparseMatrix<double> AMAu = (Au_ * Minv * Au_).pruned();
    SparseMatrix<double> AMAv = (Av_ * Minv * Av_).pruned();
    luHVu_.compute((M_ + beta * AMAu).pruned());
    luHVv_.compute((M_ + beta * AMAv).pruned());
    hvBuilt_ = true;
}

void NSIntegrator::applyHV_(VectorXd& u, VectorXd& v) const {
    if (!hvBuilt_) return;
    u = luHVu_.solve(M_ * u);
    v = luHVv_.solve(M_ * v);
}

// ---------------------------------------------------------------------------
// Semi-implicit immersed-boundary constraint (Schur-complement / IBPM).
// ---------------------------------------------------------------------------
void NSIntegrator::buildIBKernelAux_() {
    if (ibAux_.built) return;
    int NT = mesh_.elem.rows();
    ibAux_.cent.resize(NT);
    ibAux_.rad.resize(NT);
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    ibAux_.radMax = 0.0;
    for (int t = 0; t < NT; ++t) {
        Vector2d c(0, 0);
        for (int k = 0; k < 3; ++k) c += mesh_.node.row(mesh_.elem(t, k)).transpose();
        c /= 3.0;
        double r = 0.0;
        for (int k = 0; k < 3; ++k)
            r = std::max(r, (mesh_.node.row(mesh_.elem(t, k)).transpose() - c).norm());
        ibAux_.cent[t] = c; ibAux_.rad[t] = r;
        ibAux_.radMax = std::max(ibAux_.radMax, r);
        xmin = std::min(xmin, c.x()); xmax = std::max(xmax, c.x());
        ymin = std::min(ymin, c.y()); ymax = std::max(ymax, c.y());
    }
    // Bucket cell large enough that a marker's kernel support (radius delta)
    // plus the largest element circumradius is covered by the 3x3 neighborhood.
    ibAux_.cell = std::max(1e-12, 2.0 * ibAux_.radMax);
    ibAux_.xmin = xmin; ibAux_.ymin = ymin;
    ibAux_.nx = std::max(1, (int)std::ceil((xmax - xmin) / ibAux_.cell));
    ibAux_.ny = std::max(1, (int)std::ceil((ymax - ymin) / ibAux_.cell));
    ibAux_.buckets.assign((size_t)ibAux_.nx * ibAux_.ny, {});
    for (int t = 0; t < NT; ++t) {
        int i = std::min(ibAux_.nx - 1, std::max(0, (int)((ibAux_.cent[t].x() - xmin) / ibAux_.cell)));
        int j = std::min(ibAux_.ny - 1, std::max(0, (int)((ibAux_.cent[t].y() - ymin) / ibAux_.cell)));
        ibAux_.buckets[(size_t)j * ibAux_.nx + i].push_back(t);
    }
    fem_.quad2d(ibAux_.quadL, ibAux_.quadW);
    int nq = (int)ibAux_.quadW.size();
    ibAux_.quadPhi.resize(nq);
    for (int q = 0; q < nq; ++q)
        ibAux_.quadPhi[q] = fem_.computeBasisValue_all(ibAux_.quadL.row(q)).row(0);
    ibAux_.built = true;
}

void NSIntegrator::setIBConstraint(const MatrixXd& X, const MatrixXd& V,
                                   const MeshLocator& loc, double eps,
                                   int maxCG, double tol, double kernelDelta) {
    int m = (int)X.rows();
    ib_.eps = eps; ib_.maxCG = maxCG; ib_.tol = tol;
    ib_.delta = std::max(0.0, kernelDelta);
    ib_.elem.clear(); ib_.phi.clear(); ib_.rows.clear();
    std::vector<int> keep;
    keep.reserve(m);

    if (ib_.delta <= 0.0) {
        // Legacy pointwise transfer: locate each marker, sample the broken basis.
        for (int k = 0; k < m; ++k) {
            double lam[3];
            int e = loc.locate(mesh_, X(k, 0), X(k, 1), -1, lam);
            if (e < 0) continue;             // marker outside the mesh -> drop it
            ib_.elem.push_back(e);
            Vector3d L(lam[0], lam[1], lam[2]);
            ib_.phi.push_back(fem_.computeBasisValue_all(L.transpose()).row(0));
            keep.push_back(k);
        }
    } else {
        // Mollified transfer: normalised Wendland-C2 kernel rows per marker,
        //   row_e(i) = \int_e psi(|x - X_k|) phi_i dx / m_k,  m_k = sum_e \int_e psi.
        buildIBKernelAux_();
        const double delta = ib_.delta;
        const int nq = (int)ibAux_.quadW.size(), locDof = fem_.locDof;
        const double reach = delta + ibAux_.radMax;
        const int span = (int)std::ceil(reach / ibAux_.cell);
        for (int k = 0; k < m; ++k) {
            Vector2d Xk(X(k, 0), X(k, 1));
            int ci = (int)((Xk.x() - ibAux_.xmin) / ibAux_.cell);
            int cj = (int)((Xk.y() - ibAux_.ymin) / ibAux_.cell);
            std::vector<std::pair<int, RowVectorXd>> rows;
            double mass = 0.0;
            for (int j = cj - span; j <= cj + span; ++j) {
                if (j < 0 || j >= ibAux_.ny) continue;
                for (int i = ci - span; i <= ci + span; ++i) {
                    if (i < 0 || i >= ibAux_.nx) continue;
                    for (int t : ibAux_.buckets[(size_t)j * ibAux_.nx + i]) {
                        if ((ibAux_.cent[t] - Xk).norm() > delta + ibAux_.rad[t]) continue;
                        Vector2d P1 = mesh_.node.row(mesh_.elem(t, 0));
                        Vector2d P2 = mesh_.node.row(mesh_.elem(t, 1));
                        Vector2d P3 = mesh_.node.row(mesh_.elem(t, 2));
                        double area = fem_.area(t);
                        RowVectorXd row = RowVectorXd::Zero(locDof);
                        for (int q = 0; q < nq; ++q) {
                            Vector2d xq = ibAux_.quadL(q, 0) * P1 + ibAux_.quadL(q, 1) * P2
                                        + ibAux_.quadL(q, 2) * P3;
                            double r = (xq - Xk).norm();
                            if (r >= delta) continue;
                            double s = 1.0 - r / delta;          // Wendland C2: (1-q)^4 (4q+1)
                            double psi = s * s * s * s * (4.0 * r / delta + 1.0);
                            double wq = ibAux_.quadW(q) * area * psi;
                            row.noalias() += wq * ibAux_.quadPhi[q];
                            mass += wq;
                        }
                        if (row.cwiseAbs().sum() > 0.0) rows.emplace_back(t, row);
                    }
                }
            }
            if (rows.empty() || mass <= 1e-300) continue;        // kernel sees no mesh -> drop
            for (auto& pr : rows) pr.second /= mass;             // normalise: I_k(1) = 1
            ib_.rows.push_back(std::move(rows));
            keep.push_back(k);
        }
    }
    int mk = (int)keep.size();
    ib_.V.resize(mk, 2);
    for (int i = 0; i < mk; ++i) ib_.V.row(i) = V.row(keep[i]);
    ib_.active = (mk > 0);
}

VectorXd NSIntegrator::ibInterp_(const VectorXd& w) const {
    int locDof = fem_.locDof;
    if (ib_.delta > 0.0) {
        int mk = (int)ib_.rows.size();
        VectorXd r(mk);
        for (int k = 0; k < mk; ++k) {
            double s = 0.0;
            for (const auto& pr : ib_.rows[k]) {
                const RowVectorXd& ph = pr.second;
                for (int i = 0; i < locDof; ++i) s += ph(i) * w(elem2dof_(pr.first, i));
            }
            r(k) = s;
        }
        return r;
    }
    int mk = (int)ib_.elem.size();
    VectorXd r(mk);
    for (int k = 0; k < mk; ++k) {
        int e = ib_.elem[k];
        const RowVectorXd& ph = ib_.phi[k];
        double s = 0.0;
        for (int i = 0; i < locDof; ++i) s += ph(i) * w(elem2dof_(e, i));
        r(k) = s;
    }
    return r;
}

VectorXd NSIntegrator::ibScatter_(const VectorXd& h) const {
    int locDof = fem_.locDof;
    VectorXd load = VectorXd::Zero(nDof_);
    if (ib_.delta > 0.0) {
        int mk = (int)ib_.rows.size();
        for (int k = 0; k < mk; ++k)
            for (const auto& pr : ib_.rows[k]) {
                const RowVectorXd& ph = pr.second;
                for (int i = 0; i < locDof; ++i) load(elem2dof_(pr.first, i)) += ph(i) * h(k);
            }
        return load;
    }
    int mk = (int)ib_.elem.size();
    for (int k = 0; k < mk; ++k) {
        int e = ib_.elem[k];
        const RowVectorXd& ph = ib_.phi[k];
        for (int i = 0; i < locDof; ++i) load(elem2dof_(e, i)) += ph(i) * h(k);
    }
    return load;
}

void NSIntegrator::applyIBConstraint_(VectorXd& w, const VectorXd& target,
                                      const SimplicialLDLT<SparseMatrix<double>>& luH,
                                      int& iters) {
    int mk = (ib_.delta > 0.0) ? (int)ib_.rows.size() : (int)ib_.elem.size();
    iters = 0;
    if (mk == 0) return;
    // CG on  (G + eps I) h = r0,  G = I H^{-1} I^T,  r0 = I w - target.
    VectorXd r0 = ibInterp_(w) - target;
    double bnorm = r0.norm();
    if (bnorm < 1e-300) return;
    VectorXd h = VectorXd::Zero(mk);
    VectorXd res = r0, p = r0;
    double rs = res.squaredNorm();
    for (iters = 0; iters < ib_.maxCG; ) {
        VectorXd Gp = ibInterp_(luH.solve(ibScatter_(p)));    // I H^{-1} I^T p
        if (ib_.eps > 0.0) Gp.noalias() += ib_.eps * p;
        double denom = p.dot(Gp);
        ++iters;
        if (denom <= 0.0) break;                              // G is SPD: safety only
        double alpha = rs / denom;
        h.noalias()   += alpha * p;
        res.noalias() -= alpha * Gp;
        double rsNew = res.squaredNorm();
        if (std::sqrt(rsNew) <= ib_.tol * bnorm) break;
        p = res + (rsNew / rs) * p;
        rs = rsNew;
    }
    // velocity correction:  w <- w - H^{-1}(I^T h)
    w.noalias() -= luH.solve(ibScatter_(h));
}

// --- explicit Lax-Friedrichs convection load  c_m = \int N_m(u) phi --------
void NSIntegrator::assembleConvection(const VectorXd& uu, const VectorXd& vv,
                                      double t, VectorXd& cx, VectorXd& cy) const {
    int NT = mesh_.elem.rows(), NE = edge_.rows(), locDof = fem_.locDof;
    cx = VectorXd::Zero(nDof_); cy = VectorXd::Zero(nDof_);

    MatrixXd quadL; VectorXd w; fem_.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phiV(nq); std::vector<MatrixXd> dphiV(nq);
    for (int q = 0; q < nq; ++q) {
        phiV[q]  = fem_.computeBasisValue_all(quadL.row(q)).row(0);
        dphiV[q] = fem_.computeBasisDlam_all(quadL.row(q));
    }
    auto getElem = [&](const VectorXd& f, int tt, VectorXd& fe) {
        fe.resize(locDof); for (int i = 0; i < locDof; ++i) fe(i) = f(elem2dof_(tt, i));
    };

    // volume:  c_m -= \int (u_m * U) . grad phi_i
    VectorXd ue(locDof), ve(locDof);
    for (int t0 = 0; t0 < NT; ++t0) {
        getElem(uu, t0, ue); getElem(vv, t0, ve);
        double area = fem_.area(t0);
        VectorXd cex = VectorXd::Zero(locDof), cey = VectorXd::Zero(locDof);
        for (int q = 0; q < nq; ++q) {
            double U = phiV[q].dot(ue), V = phiV[q].dot(ve);
            MatrixXd gp = fem_.Dlam[t0] * dphiV[q];     // 2 x locDof
            double wa = w(q) * area;
            // F_x = U*(U,V), F_y = V*(U,V)
            RowVectorXd Fx_dot = (U * U) * gp.row(0) + (U * V) * gp.row(1);
            RowVectorXd Fy_dot = (V * U) * gp.row(0) + (V * V) * gp.row(1);
            cex.noalias() -= wa * Fx_dot.transpose();
            cey.noalias() -= wa * Fy_dot.transpose();
        }
        for (int i = 0; i < locDof; ++i) { cx(elem2dof_(t0, i)) += cex(i); cy(elem2dof_(t0, i)) += cey(i); }
    }

    // faces: local Lax-Friedrichs flux
    MatrixXd quad1d; VectorXd w1d; fem_.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem_, quad1d);
    for (int e = 0; e < NE; ++e) {
        int t1 = edge2side_(e, 0), t2 = edge2side_(e, 1);
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        int ta = (t1 != -1) ? t1 : t2;             // "interior" side of this edge
        EdgeOnElem ea = edgeOnElem(mesh_, ta, n1, n2);
        Vector2d n = ea.nout; double he = ea.he;   // normal outward from ta

        VectorXd u1, v1; getElem(uu, ta, u1); getElem(vv, ta, v1);
        bool interior = (t1 != -1 && t2 != -1);
        VectorXd cAx = VectorXd::Zero(locDof), cAy = VectorXd::Zero(locDof);
        VectorXd cBx = VectorXd::Zero(locDof), cBy = VectorXd::Zero(locDof);
        EdgeOnElem eb{};
        VectorXd u2, v2;
        if (interior) { getElem(uu, t2, u2); getElem(vv, t2, v2); eb = edgeOnElem(mesh_, t2, n1, n2); }

        // BC kind on this boundary edge (if any)
        int bu = bc_.bcU(e), bv = bc_.bcV(e);
        bool isOutflow = (!interior && bu == BCN_NEUMANN && bv == BCN_NEUMANN);
        bool isWall    = (!interior && bu == BCN_NEUMANN && bv == BCN_DIRICHLET);
        bool isDirVel  = (!interior && bu == BCN_DIRICHLET && bv == BCN_DIRICHLET);

        for (int q = 0; q < E.nq; ++q) {
            const RowVectorXd& pa = E.phi[ea.et][ea.dir][q];
            double Um = pa.dot(u1), Vm = pa.dot(v1);           // interior trace
            double Up, Vp;                                      // exterior trace
            const RowVectorXd* pb = nullptr;
            if (interior) {
                pb = &E.phi[eb.et][eb.dir][q];
                Up = pb->dot(u2); Vp = pb->dot(v2);
            } else if (isDirVel) {
                double l1 = quad1d(q, 0), l2 = quad1d(q, 1);
                Vector2d P = l1 * mesh_.node.row(n1).transpose() + l2 * mesh_.node.row(n2).transpose();
                Vector2d ub = bc_.velDir(P.x(), P.y(), t);
                Up = ub.x(); Vp = ub.y();
            } else { // outflow / wall: exterior = interior (handled specially for wall below)
                Up = Um; Vp = Vm;
            }

            // numerical normal flux for momentum components
            double fAx, fAy;
            if (isWall) {            // slip wall: u.n = 0 -> no convective flux through
                fAx = 0.0; fAy = 0.0;
            } else {
                double FnA_x = (Um * Um) * n.x() + (Um * Vm) * n.y();   // F_x^- . n
                double FnA_y = (Vm * Um) * n.x() + (Vm * Vm) * n.y();
                double FnB_x = (Up * Up) * n.x() + (Up * Vp) * n.y();   // F_x^+ . n
                double FnB_y = (Vp * Up) * n.x() + (Vp * Vp) * n.y();
                double lam = std::max(std::abs(Um * n.x() + Vm * n.y()),
                                      std::abs(Up * n.x() + Vp * n.y()));
                fAx = 0.5 * (FnA_x + FnB_x) + 0.5 * lam * (Um - Up);
                fAy = 0.5 * (FnA_y + FnB_y) + 0.5 * lam * (Vm - Vp);
            }
            double whe = w1d(q) * he;
            cAx.noalias() += (whe * fAx) * pa.transpose();
            cAy.noalias() += (whe * fAy) * pa.transpose();
            if (interior) {          // opposite flux on the t2 side
                cBx.noalias() += (whe * (-fAx)) * pb->transpose();
                cBy.noalias() += (whe * (-fAy)) * pb->transpose();
            }
        }
        for (int i = 0; i < locDof; ++i) { cx(elem2dof_(ta, i)) += cAx(i); cy(elem2dof_(ta, i)) += cAy(i); }
        if (interior)
            for (int i = 0; i < locDof; ++i) { cx(elem2dof_(t2, i)) += cBx(i); cy(elem2dof_(t2, i)) += cBy(i); }
    }
}

// --- high-order pressure Neumann load on bcP == 2 edges --------------------
VectorXd NSIntegrator::assemblePressureNeumann(const VectorXd& us, const VectorXd& vs,
                                               double tEnd) const {
    int NE = edge_.rows(), locDof = fem_.locDof;
    VectorXd b = VectorXd::Zero(nDof_);
    const bool bdf2 = (n_ >= 1);
    MatrixXd quad1d; VectorXd w1d; fem_.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem_, quad1d);

    const double invSqrt2 = 1.0 / std::sqrt(2.0);
    const bool rot = rotationalPressureBC;
    // Returns the convective term N=(u.grad)u and the viscous operator Vis used in
    // the pressure BC (the caller multiplies Vis by nu).  With `rot` (default) Vis
    // is the rotational form  -curl(omega) = (u_yy - v_xy, v_xx - u_xy); otherwise
    // the Laplacian form  Lap u  (the latter only reaches 1st order -- see header).
    auto traceQuantities = [&](int t, const EdgeOnElem& ei, int q,
                               const VectorXd& uf, const VectorXd& vf,
                               double& U, double& V, double& Nx, double& Ny,
                               double& VisU, double& VisV) {
        VectorXd ue(locDof), ve(locDof);
        for (int i = 0; i < locDof; ++i) { ue(i) = uf(elem2dof_(t, i)); ve(i) = vf(elem2dof_(t, i)); }
        const RowVectorXd& phi = E.phi[ei.et][ei.dir][q];
        MatrixXd g = fem_.Dlam[t] * E.dphi[ei.et][ei.dir][q];     // 2 x locDof
        MatrixXd H = fem_.R[t] * E.Hlam[ei.et][ei.dir][q];        // 3 x locDof Voigt [xx,yy,sqrt2 xy]
        U = phi.dot(ue); V = phi.dot(ve);
        double ux = g.row(0).dot(ue), uy = g.row(1).dot(ue);
        double vx = g.row(0).dot(ve), vy = g.row(1).dot(ve);
        Nx = U * ux + V * uy;  Ny = U * vx + V * vy;              // (u.grad)u
        double uxx = H.row(0).dot(ue), uyy = H.row(1).dot(ue), uxy = invSqrt2 * H.row(2).dot(ue);
        double vxx = H.row(0).dot(ve), vyy = H.row(1).dot(ve), vxy = invSqrt2 * H.row(2).dot(ve);
        if (rot) { VisU = uyy - vxy; VisV = vxx - uxy; }          // -curl(omega)
        else     { VisU = uxx + uyy; VisV = vxx + vyy; }          // Lap u
    };

    for (int e = 0; e < NE; ++e) {
        if (bc_.bcP(e) != BCN_NEUMANN_HO) continue;
        int t = (edge2side_(e, 0) != -1) ? edge2side_(e, 0) : edge2side_(e, 1);
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        EdgeOnElem ei = edgeOnElem(mesh_, t, n1, n2);
        VectorXd be = VectorXd::Zero(locDof);
        for (int q = 0; q < E.nq; ++q) {
            double Un, Vn, Nxn, Nyn, LuN, LvN;
            traceQuantities(t, ei, q, us, vs, Un, Vn, Nxn, Nyn, LuN, LvN);
            double Nx = Nxn, Ny = Nyn, Lu = LuN, Lv = LvN;
            if (bdf2) {
                double Uo, Vo, Nxo, Nyo, LuO, LvO;
                traceQuantities(t, ei, q, uPrev_, vPrev_, Uo, Vo, Nxo, Nyo, LuO, LvO);
                Nx = 2 * Nxn - Nxo; Ny = 2 * Nyn - Nyo;
                Lu = 2 * LuN - LuO; Lv = 2 * LvN - LvO;
            }
            if (!includeConvection) { Nx = Ny = 0; }   // diagnostic: unsteady Stokes
            double l1 = quad1d(q, 0), l2 = quad1d(q, 1);
            Vector2d P = l1 * mesh_.node.row(n1).transpose() + l2 * mesh_.node.row(n2).transpose();
            double gN;
            if (bc_.presGrad) {                 // diagnostic: exact dp/dn
                Vector2d gp = bc_.presGrad(P.x(), P.y(), tEnd);
                gN = ei.nout.x() * gp.x() + ei.nout.y() * gp.y();
            } else {
                Vector2d ab = bc_.velAcc ? bc_.velAcc(P.x(), P.y(), tEnd) : Vector2d(0, 0);
                double gNx = -Nx + nu_ * Lu - ab.x();
                double gNy = -Ny + nu_ * Lv - ab.y();
                gN = ei.nout.x() * gNx + ei.nout.y() * gNy;
            }
            be.noalias() += (w1d(q) * ei.he * gN) * E.phi[ei.et][ei.dir][q].transpose();
        }
        for (int i = 0; i < locDof; ++i) b(elem2dof_(t, i)) += be(i);
    }
    return b;
}

// --- Nitsche Dirichlet RHS lift for a velocity component -------------------
VectorXd NSIntegrator::velDirichletLift(int comp, double tEnd) const {
    int NE = edge_.rows(), locDof = fem_.locDof;
    VectorXd b = VectorXd::Zero(nDof_);
    const VectorXi& bccomp = (comp == 0) ? bc_.bcU : bc_.bcV;
    MatrixXd quad1d; VectorXd w1d; fem_.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem_, quad1d);
    for (int e = 0; e < NE; ++e) {
        if (bccomp(e) != BCN_DIRICHLET) continue;
        int t = (edge2side_(e, 0) != -1) ? edge2side_(e, 0) : edge2side_(e, 1);
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        EdgeOnElem ei = edgeOnElem(mesh_, t, n1, n2);
        double pen = sigma_ / ei.he;
        VectorXd be = VectorXd::Zero(locDof);
        for (int q = 0; q < E.nq; ++q) {
            const RowVectorXd& phi = E.phi[ei.et][ei.dir][q];
            MatrixXd g = fem_.Dlam[t] * E.dphi[ei.et][ei.dir][q];
            RowVectorXd dn = ei.nout.transpose() * g;
            double l1 = quad1d(q, 0), l2 = quad1d(q, 1);
            Vector2d P = l1 * mesh_.node.row(n1).transpose() + l2 * mesh_.node.row(n2).transpose();
            double gval = bc_.velDir(P.x(), P.y(), tEnd)(comp);
            double whe = w1d(q) * ei.he;
            be.noalias() += whe * gval * (pen * phi.transpose() - beta_ * dn.transpose());
        }
        for (int i = 0; i < locDof; ++i) b(elem2dof_(t, i)) += be(i);
    }
    return b;
}

// --- Nitsche Dirichlet RHS lift for pressure (bcP == 1) --------------------
VectorXd NSIntegrator::presDirichletLift(double tEnd) const {
    int NE = edge_.rows(), locDof = fem_.locDof;
    VectorXd b = VectorXd::Zero(nDof_);
    MatrixXd quad1d; VectorXd w1d; fem_.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem_, quad1d);
    for (int e = 0; e < NE; ++e) {
        if (bc_.bcP(e) != BCN_DIRICHLET) continue;
        int t = (edge2side_(e, 0) != -1) ? edge2side_(e, 0) : edge2side_(e, 1);
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        EdgeOnElem ei = edgeOnElem(mesh_, t, n1, n2);
        double pen = sigma_ / ei.he;
        VectorXd be = VectorXd::Zero(locDof);
        for (int q = 0; q < E.nq; ++q) {
            const RowVectorXd& phi = E.phi[ei.et][ei.dir][q];
            MatrixXd g = fem_.Dlam[t] * E.dphi[ei.et][ei.dir][q];
            RowVectorXd dn = ei.nout.transpose() * g;
            double l1 = quad1d(q, 0), l2 = quad1d(q, 1);
            Vector2d P = l1 * mesh_.node.row(n1).transpose() + l2 * mesh_.node.row(n2).transpose();
            double gval = bc_.presDir ? bc_.presDir(P.x(), P.y(), tEnd) : 0.0;
            double whe = w1d(q) * ei.he;
            be.noalias() += whe * gval * (pen * phi.transpose() - beta_ * dn.transpose());
        }
        for (int i = 0; i < locDof; ++i) b(elem2dof_(t, i)) += be(i);
    }
    return b;
}

// --- post-processing -------------------------------------------------------
VectorXd NSIntegrator::vorticity() const {
    // omega = dv/dx - du/dy  projected:  M omega = Gx v - Gy u
    return luM_.solve(Gx_ * v_ - Gy_ * u_);
}

VectorXd NSIntegrator::divergence() const {
    // div u = du/dx + dv/dy  projected:  M d = Gx u + Gy v
    return luM_.solve(Gx_ * u_ + Gy_ * v_);
}

VectorXd NSIntegrator::speed() const {
    VectorXd s(nDof_);
    for (int i = 0; i < nDof_; ++i) s(i) = std::hypot(u_(i), v_(i));
    return s;
}

void NSIntegrator::cylinderForce(const VectorXi& isCylEdge, double& Fx, double& Fy) const {
    int NE = edge_.rows(), locDof = fem_.locDof;
    Fx = 0.0; Fy = 0.0;
    MatrixXd quad1d; VectorXd w1d; fem_.quad1d(quad1d, w1d);
    EdgeBasis E = makeEdgeBasis(fem_, quad1d);
    for (int e = 0; e < NE; ++e) {
        if (!isCylEdge(e)) continue;
        int t = (edge2side_(e, 0) != -1) ? edge2side_(e, 0) : edge2side_(e, 1);
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        EdgeOnElem ei = edgeOnElem(mesh_, t, n1, n2);
        VectorXd ue(locDof), ve(locDof), pe(locDof);
        for (int i = 0; i < locDof; ++i) { ue(i) = u_(elem2dof_(t, i)); ve(i) = v_(elem2dof_(t, i)); pe(i) = p_(elem2dof_(t, i)); }
        for (int q = 0; q < E.nq; ++q) {
            const RowVectorXd& phi = E.phi[ei.et][ei.dir][q];
            MatrixXd g = fem_.Dlam[t] * E.dphi[ei.et][ei.dir][q];
            double ux = g.row(0).dot(ue), uy = g.row(1).dot(ue);
            double vx = g.row(0).dot(ve), vy = g.row(1).dot(ve);
            double pv = phi.dot(pe);
            // ei.nout is the fluid element's outward normal, which on a cylinder
            // edge points INTO the cylinder; the force on the body integrates the
            // traction against the cylinder's outward (into-fluid) normal = -ei.nout.
            double nx = -ei.nout.x(), ny = -ei.nout.y();
            // traction t = sigma.n, sigma = -p I + nu(grad u + grad u^T)
            double tx = -pv * nx + nu_ * ((2 * ux) * nx + (uy + vx) * ny);
            double ty = -pv * ny + nu_ * ((uy + vx) * nx + (2 * vy) * ny);
            double whe = w1d(q) * ei.he;
            Fx += whe * tx; Fy += whe * ty;
        }
    }
}

// ===========================================================================
// Visualisation
// ===========================================================================
namespace {
void colormapRGB(double t, Colormap cm, unsigned char& R, unsigned char& G, unsigned char& B) {
    t = std::min(1.0, std::max(0.0, t));
    double col[3];
    if (cm == CM_COOLWARM) {
        const double lo[3] = {59, 76, 192}, mid[3] = {242, 242, 242}, hi[3] = {180, 4, 38};
        if (t < 0.5) { double s = t / 0.5; for (int k = 0; k < 3; ++k) col[k] = lo[k] + s * (mid[k] - lo[k]); }
        else { double s = (t - 0.5) / 0.5; for (int k = 0; k < 3; ++k) col[k] = mid[k] + s * (hi[k] - mid[k]); }
    } else { // viridis-ish 4-stop
        const double c0[3]={68,1,84}, c1[3]={59,82,139}, c2[3]={33,145,140}, c3[3]={94,201,98}, c4[3]={253,231,37};
        const double* C[5] = {c0,c1,c2,c3,c4};
        double s = t * 4.0; int k = std::min(3, (int)s); double f = s - k;
        for (int j = 0; j < 3; ++j) col[j] = C[k][j] + f * (C[k+1][j] - C[k][j]);
    }
    R = (unsigned char)std::lround(col[0]); G = (unsigned char)std::lround(col[1]); B = (unsigned char)std::lround(col[2]);
}
}

void writeFieldPPM(const std::string& path, FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                   const VectorXd& field, int Wpix, int Hpix,
                   double xmin, double xmax, double ymin, double ymax,
                   double vmin, double vmax, Colormap cmap,
                   const std::function<bool(double, double)>& inDomain) {
    const int W = Wpix, H = Hpix;
    int locDof = fem.locDof;
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    std::vector<unsigned char> img((size_t)W * H * 3, 255);
    int NT = mesh.elem.rows();
    // Represent the element field over the monomial (barycentric) basis once:
    //   field(lam) = sum_j mco_j * lam1^a_j lam2^b_j lam3^c_j ,  mco = coef * fe.
    // Then each pixel is a cheap, allocation-free polynomial evaluation (smooth dP_k
    // interior shading, not one colour per element).
    MatrixXi midx = numSplit3(fem.ord);   // 3 x locDof exponents
    VectorXd fe(locDof), mco(locDof);
    for (int t = 0; t < NT; ++t) {
        Vector2d P1 = mesh.node.row(mesh.elem(t, 0));
        Vector2d P2 = mesh.node.row(mesh.elem(t, 1));
        Vector2d P3 = mesh.node.row(mesh.elem(t, 2));
        for (int i = 0; i < locDof; ++i) fe(i) = field(elem2dof(t, i));   // element field coeffs
        mco.noalias() = fem.coef * fe;
        Vector2d v0 = P2 - P1, v1 = P3 - P1;
        double det = v0.x() * v1.y() - v1.x() * v0.y();
        if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;
        double txmin = std::min({P1.x(), P2.x(), P3.x()}), txmax = std::max({P1.x(), P2.x(), P3.x()});
        double tymin = std::min({P1.y(), P2.y(), P3.y()}), tymax = std::max({P1.y(), P2.y(), P3.y()});
        int col0 = std::max(0, (int)std::floor((txmin - xmin) / dx) - 1);
        int col1 = std::min(W - 1, (int)std::ceil((txmax - xmin) / dx) + 1);
        int row0 = std::max(0, (int)std::floor((ymax - tymax) / dy) - 1);
        int row1 = std::min(H - 1, (int)std::ceil((ymax - tymin) / dy) + 1);
        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double x = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double y = ymax - (row + 0.5) * dy;
                double px = x - P1.x(), py = y - P1.y();
                double l2 = (px * v1.y() - v1.x() * py) * invd;
                double l3 = (v0.x() * py - px * v0.y()) * invd;
                double l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                if (inDomain && !inDomain(x, y)) continue;
                double val = 0.0;                              // evaluate the dP_k polynomial
                for (int j = 0; j < locDof; ++j) {
                    double s = mco(j);
                    for (int e = 0; e < midx(0, j); ++e) s *= l1;
                    for (int e = 0; e < midx(1, j); ++e) s *= l2;
                    for (int e = 0; e < midx(2, j); ++e) s *= l3;
                    val += s;
                }
                double s = (val - vmin) / (vmax - vmin);
                unsigned char R, G, B; colormapRGB(s, cmap, R, G, B);
                size_t idx = ((size_t)row * W + col) * 3;
                img[idx] = R; img[idx + 1] = G; img[idx + 2] = B;
            }
        }
    }
    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "writeFieldPPM: cannot open " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

// ===========================================================================
// Tracer-particle visualisation: O(1)-amortised mesh locator + RK2 advection
// + composited frame writer.  See NavierStokes.h for the API contract.
// ===========================================================================

namespace {

// Triangle barycentric coords for (x,y) in element t; tol allows mild outside slack.
inline bool baryInElem(const Mesh& mesh, int t, double x, double y, double lam[3], double tol) {
    Vector2d P1 = mesh.node.row(mesh.elem(t, 0));
    Vector2d P2 = mesh.node.row(mesh.elem(t, 1));
    Vector2d P3 = mesh.node.row(mesh.elem(t, 2));
    Vector2d v0 = P2 - P1, v1 = P3 - P1;
    double det = v0.x() * v1.y() - v1.x() * v0.y();
    if (std::abs(det) < 1e-300) return false;
    double invd = 1.0 / det;
    double px = x - P1.x(), py = y - P1.y();
    double l2 = (px * v1.y() - v1.x() * py) * invd;
    double l3 = (v0.x() * py - px * v0.y()) * invd;
    double l1 = 1.0 - l2 - l3;
    if (l1 < tol || l2 < tol || l3 < tol) return false;
    lam[0] = l1; lam[1] = l2; lam[2] = l3;
    return true;
}

// Evaluate a scalar dP_k field at barycentric `lam` of element `t`.
inline double evalFieldAt(FEM& fem, const MatrixXi& elem2dof, const VectorXd& f,
                          int t, const double lam[3]) {
    int locDof = fem.locDof;
    Vector3d L(lam[0], lam[1], lam[2]);
    RowVectorXd phi = fem.computeBasisValue_all(L.transpose()).row(0);
    double s = 0.0;
    for (int i = 0; i < locDof; ++i) s += phi(i) * f(elem2dof(t, i));
    return s;
}

// Evaluate the velocity (u,v) at point (x,y); returns false if outside the mesh.
inline bool evalVelocity(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                         const VectorXd& uF, const VectorXd& vF,
                         const MeshLocator& loc, double x, double y,
                         int& hint, double& u, double& v) {
    double lam[3];
    int t = loc.locate(mesh, x, y, hint, lam);
    if (t < 0) return false;
    hint = t;
    u = evalFieldAt(fem, elem2dof, uF, t, lam);
    v = evalFieldAt(fem, elem2dof, vF, t, lam);
    return true;
}

inline std::uint32_t xorshift32(std::uint32_t& s) {
    s ^= s << 13; s ^= s >> 17; s ^= s << 5;
    return s ? s : (s = 0x9E3779B9u);
}
inline double frand(std::uint32_t& s) {
    return (xorshift32(s) >> 8) * (1.0 / 16777216.0); // [0,1)
}

} // anon

void MeshLocator::build(const Mesh& mesh, int nx) {
    int NT = mesh.elem.rows();
    if (NT == 0) { nx_ = ny_ = 0; cells_.clear(); return; }
    double xMin = mesh.node.col(0).minCoeff(), xMax = mesh.node.col(0).maxCoeff();
    double yMin = mesh.node.col(1).minCoeff(), yMax = mesh.node.col(1).maxCoeff();
    // tiny pad to keep boundary points strictly inside
    double pad = 1e-9 * std::max(xMax - xMin, yMax - yMin) + 1e-12;
    xMin -= pad; xMax += pad; yMin -= pad; yMax += pad;
    if (nx <= 0) nx = std::max(8, (int)std::round(std::sqrt((double)NT)));
    double aspect = (yMax - yMin) / std::max(1e-300, (xMax - xMin));
    int ny = std::max(4, (int)std::round(nx * aspect));
    nx_ = nx; ny_ = ny;
    xmin_ = xMin; ymin_ = yMin;
    invDx_ = nx_ / (xMax - xMin);
    invDy_ = ny_ / (yMax - yMin);
    cells_.assign((size_t)nx_ * ny_, {});
    for (int t = 0; t < NT; ++t) {
        Vector2d P1 = mesh.node.row(mesh.elem(t, 0));
        Vector2d P2 = mesh.node.row(mesh.elem(t, 1));
        Vector2d P3 = mesh.node.row(mesh.elem(t, 2));
        double txmin = std::min({P1.x(), P2.x(), P3.x()});
        double txmax = std::max({P1.x(), P2.x(), P3.x()});
        double tymin = std::min({P1.y(), P2.y(), P3.y()});
        double tymax = std::max({P1.y(), P2.y(), P3.y()});
        int c0 = std::max(0, (int)std::floor((txmin - xmin_) * invDx_));
        int c1 = std::min(nx_ - 1, (int)std::floor((txmax - xmin_) * invDx_));
        int r0 = std::max(0, (int)std::floor((tymin - ymin_) * invDy_));
        int r1 = std::min(ny_ - 1, (int)std::floor((tymax - ymin_) * invDy_));
        for (int r = r0; r <= r1; ++r)
            for (int c = c0; c <= c1; ++c)
                cells_[(size_t)r * nx_ + c].push_back(t);
    }
}

int MeshLocator::locate(const Mesh& mesh, double x, double y, int hint, double lam[3]) const {
    const double tol = -1e-10;
    if (hint >= 0 && hint < (int)mesh.elem.rows() && baryInElem(mesh, hint, x, y, lam, tol))
        return hint;
    if (nx_ <= 0 || ny_ <= 0) return -1;
    int c = (int)std::floor((x - xmin_) * invDx_);
    int r = (int)std::floor((y - ymin_) * invDy_);
    if (c < 0 || c >= nx_ || r < 0 || r >= ny_) return -1;
    const std::vector<int>& bucket = cells_[(size_t)r * nx_ + c];
    for (int t : bucket)
        if (baryInElem(mesh, t, x, y, lam, tol)) return t;
    return -1;
}

void ParticleTracer::reset() {
    int L = std::max(1, trailLen);
    x.assign(N, 0.0); y.assign(N, 0.0);
    age.assign(N, 0); alive.assign(N, 0);
    hint.assign(N, -1);
    trailX.assign((size_t)N * L, 0.0);
    trailY.assign((size_t)N * L, 0.0);
    trailValid.assign((size_t)N * L, 0);
    trailHead.assign(N, 0);
    strideTick = 0;
    // Scatter uniformly inside the window; reject the cylinder hole.
    for (int i = 0; i < N; ++i) {
        for (int tries = 0; tries < 32; ++tries) {
            double rx = xa + (xb - xa) * frand(rng);
            double ry = ya + (yb - ya) * frand(rng);
            if (!inDomain || inDomain(rx, ry)) {
                x[i] = rx; y[i] = ry;
                trailX[(size_t)i * L] = rx; trailY[(size_t)i * L] = ry;
                trailValid[(size_t)i * L] = 1;
                trailHead[i] = 0;
                alive[i] = 1;
                age[i] = (int)(frand(rng) * 30.0); // staggered ages -> trail variety
                break;
            }
        }
    }
}

void ParticleTracer::advance(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                             const VectorXd& uField, const VectorXd& vField,
                             const MeshLocator& loc, double dt) {
    const double pad = 1e-6;
    const int L = std::max(1, trailLen);
    const int stride = std::max(1, trailStride);
    bool record = ((++strideTick % stride) == 0);
    auto resetTrail = [&](int i, double rx, double ry) {
        size_t base = (size_t)i * L;
        for (int k = 0; k < L; ++k) trailValid[base + k] = 0;
        trailHead[i] = 0;
        trailX[base] = rx; trailY[base] = ry;
        trailValid[base] = 1;
    };
    for (int i = 0; i < N; ++i) {
        if (!alive[i]) {
            // Re-seed at the inflow rake (vertical line) at a random y.
            for (int tries = 0; tries < 16; ++tries) {
                double ry = ya + (yb - ya) * frand(rng);
                if (!inDomain || inDomain(inflowX, ry)) {
                    x[i] = inflowX; y[i] = ry;
                    resetTrail(i, inflowX, ry);
                    alive[i] = 1; age[i] = 0; hint[i] = -1;
                    break;
                }
            }
            if (!alive[i]) continue;
        }
        // RK2: k1 = u(p), k2 = u(p + dt*k1); p_new = p + 0.5*dt*(k1+k2).
        double u1, v1; int h = hint[i];
        if (!evalVelocity(fem, mesh, elem2dof, uField, vField, loc, x[i], y[i], h, u1, v1)) {
            alive[i] = 0; continue;
        }
        double xm = x[i] + dt * u1, ym = y[i] + dt * v1;
        double u2, v2; int h2 = h;
        if (!evalVelocity(fem, mesh, elem2dof, uField, vField, loc, xm, ym, h2, u2, v2)) {
            // mid-point fell outside -> Euler fallback, but keep going
            u2 = u1; v2 = v1; h2 = h;
        }
        x[i] += 0.5 * dt * (u1 + u2);
        y[i] += 0.5 * dt * (v1 + v2);
        hint[i] = h2;
        ++age[i];
        // Kill on leaving the render window or the fluid domain.
        if (x[i] < xa - pad || x[i] > xb + pad || y[i] < ya - pad || y[i] > yb + pad) {
            alive[i] = 0;
        } else if (inDomain && !inDomain(x[i], y[i])) {
            alive[i] = 0;
        }
        if (alive[i] && record) {
            int head = (trailHead[i] + 1) % L;
            trailHead[i] = head;
            size_t idx = (size_t)i * L + head;
            trailX[idx] = x[i]; trailY[idx] = y[i];
            trailValid[idx] = 1;
        }
    }
}

void writeParticlesPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                       const MatrixXi& elem2dof, const VectorXd& speedField,
                       int Wpix, int Hpix,
                       double xmin, double xmax, double ymin, double ymax,
                       double speedMin, double speedMax,
                       const ParticleTracer& tracer,
                       const std::function<bool(double, double)>& inDomain,
                       double bgDim) {
    const int W = Wpix, H = Hpix;
    int locDof = fem.locDof;
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    std::vector<unsigned char> img((size_t)W * H * 3, 255);
    if (bgDim < 0) bgDim = 0; if (bgDim > 1) bgDim = 1;

    // ---- 1) deeply dimmed |u| background ----
    int NT = mesh.elem.rows();
    MatrixXi midx = numSplit3(fem.ord);
    VectorXd fe(locDof), mco(locDof);
    for (int t = 0; t < NT; ++t) {
        Vector2d P1 = mesh.node.row(mesh.elem(t, 0));
        Vector2d P2 = mesh.node.row(mesh.elem(t, 1));
        Vector2d P3 = mesh.node.row(mesh.elem(t, 2));
        for (int i = 0; i < locDof; ++i) fe(i) = speedField(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        Vector2d v0 = P2 - P1, v1 = P3 - P1;
        double det = v0.x() * v1.y() - v1.x() * v0.y();
        if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;
        double txmin = std::min({P1.x(), P2.x(), P3.x()}), txmax = std::max({P1.x(), P2.x(), P3.x()});
        double tymin = std::min({P1.y(), P2.y(), P3.y()}), tymax = std::max({P1.y(), P2.y(), P3.y()});
        int col0 = std::max(0, (int)std::floor((txmin - xmin) / dx) - 1);
        int col1 = std::min(W - 1, (int)std::ceil((txmax - xmin) / dx) + 1);
        int row0 = std::max(0, (int)std::floor((ymax - tymax) / dy) - 1);
        int row1 = std::min(H - 1, (int)std::ceil((ymax - tymin) / dy) + 1);
        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double xq = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double yq = ymax - (row + 0.5) * dy;
                double px = xq - P1.x(), py = yq - P1.y();
                double l2 = (px * v1.y() - v1.x() * py) * invd;
                double l3 = (v0.x() * py - px * v0.y()) * invd;
                double l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                if (inDomain && !inDomain(xq, yq)) continue;
                double val = 0.0;
                for (int j = 0; j < locDof; ++j) {
                    double s = mco(j);
                    for (int e = 0; e < midx(0, j); ++e) s *= l1;
                    for (int e = 0; e < midx(1, j); ++e) s *= l2;
                    for (int e = 0; e < midx(2, j); ++e) s *= l3;
                    val += s;
                }
                double s = (val - speedMin) / (speedMax - speedMin);
                unsigned char R, G, B; colormapRGB(s, CM_VIRIDIS, R, G, B);
                size_t idx = ((size_t)row * W + col) * 3;
                img[idx]     = (unsigned char)(R * bgDim);
                img[idx + 1] = (unsigned char)(G * bgDim);
                img[idx + 2] = (unsigned char)(B * bgDim);
            }
        }
    }

    // ---- 2) draw fading ribbon trails ----
    auto blendPx = [&](int col, int row, double a, unsigned char ink) {
        if (col < 0 || col >= W || row < 0 || row >= H) return;
        if (a < 0) a = 0; if (a > 1) a = 1;
        size_t idx = ((size_t)row * W + col) * 3;
        for (int c = 0; c < 3; ++c) {
            unsigned v = (unsigned)((1.0 - a) * img[idx + c] + a * ink + 0.5);
            img[idx + c] = (unsigned char)(v > 255 ? 255 : v);
        }
    };
    auto blendLine = [&](double x0, double y0, double x1, double y1, double a, unsigned char ink) {
        // Bresenham, alpha-blended.
        int c0 = (int)std::lround((x0 - xmin) / dx - 0.5);
        int rr0 = (int)std::lround((ymax - y0) / dy - 0.5);
        int c1 = (int)std::lround((x1 - xmin) / dx - 0.5);
        int rr1 = (int)std::lround((ymax - y1) / dy - 0.5);
        int dxp = std::abs(c1 - c0), sx = c0 < c1 ? 1 : -1;
        int dyp = -std::abs(rr1 - rr0), sy = rr0 < rr1 ? 1 : -1;
        int err = dxp + dyp;
        int cur_c = c0, cur_r = rr0;
        for (int safety = 0; safety < 4 * (dxp - dyp + 2); ++safety) {
            blendPx(cur_c, cur_r, a, ink);
            if (cur_c == c1 && cur_r == rr1) break;
            int e2 = 2 * err;
            if (e2 >= dyp) { err += dyp; cur_c += sx; }
            if (e2 <= dxp) { err += dxp; cur_r += sy; }
        }
    };

    const int L = std::max(1, tracer.trailLen);
    for (int i = 0; i < tracer.N; ++i) {
        if (!tracer.alive[i]) continue;
        // newest-first walk through the ring; segments fade out toward the tail.
        int head = tracer.trailHead[i];
        size_t base = (size_t)i * L;
        // Spawn fade-in keeps freshly seeded particles from popping.
        double spawnFade = std::min(1.0, tracer.age[i] / 6.0);
        // Head dot (anti-aliased a touch via two-pixel cluster)
        {
            double xp = tracer.x[i], yp = tracer.y[i];
            int col = (int)std::lround((xp - xmin) / dx - 0.5);
            int row = (int)std::lround((ymax - yp) / dy - 0.5);
            blendPx(col, row, 0.95 * spawnFade, 255);
            blendPx(col + 1, row, 0.55 * spawnFade, 255);
            blendPx(col, row + 1, 0.55 * spawnFade, 255);
        }
        // Walk the trail backwards (newest -> oldest) drawing fading segments.
        double xa1 = tracer.x[i], ya1 = tracer.y[i];
        for (int k = 0; k < L; ++k) {
            int slot = (head - k + L * 2) % L;
            if (!tracer.trailValid[base + slot]) break;
            double xb1 = tracer.trailX[base + slot];
            double yb1 = tracer.trailY[base + slot];
            // segment from sample k-1 (closer to head) to sample k (older) -> fade by depth k
            double depth = (double)k / (double)(L - 1 < 1 ? 1 : L - 1);
            double a = (1.0 - depth) * 0.85 * spawnFade;
            blendLine(xa1, ya1, xb1, yb1, a, 255);
            xa1 = xb1; ya1 = yb1;
        }
    }

    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "writeParticlesPPM: cannot open " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

} // namespace ns
