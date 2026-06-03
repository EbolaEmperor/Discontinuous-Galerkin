#include "Argyris.h"
#include "FEM.h"        // gradbasis_my
#include "utils.h"      // polyBasisHomo3D / Grad / Hess
#include "Quadrature.h"

#include <map>
#include <array>
#include <vector>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::RowVectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseMatrix;
using Eigen::Triplet;

namespace {

// local edge ee joins local vertices (a,b); same convention as the rest of the repo
constexpr int kEdgeNodes[3][2] = {{1, 2}, {2, 0}, {0, 1}};

// quadrature order used everywhere (Hessian^2 is degree 6; load uses a smooth f)
constexpr int kQuadOrder = 12;

// barycentric coordinates of the 6 nodal sites: 3 vertices then 3 edge midpoints
void nodalSites(std::array<Vector3d, 6>& s) {
    s[0] = Vector3d(1, 0, 0);
    s[1] = Vector3d(0, 1, 0);
    s[2] = Vector3d(0, 0, 1);
    for (int e = 0; e < 3; ++e) {
        int a = kEdgeNodes[e][0], b = kEdgeNodes[e][1];
        Vector3d m = Vector3d::Zero();
        m(a) = 0.5; m(b) = 0.5;
        s[3 + e] = m;
    }
}

// barycentric-Hessian (6) -> physical [Hxx;Hyy;Hxy] map for a given Dlam (2x3)
MatrixXd buildRplain(const MatrixXd& Dl) {
    double a1 = Dl(0, 0), a2 = Dl(0, 1), a3 = Dl(0, 2);
    double b1 = Dl(1, 0), b2 = Dl(1, 1), b3 = Dl(1, 2);
    MatrixXd R(3, 6);
    R(0, 0) = a1 * a1; R(0, 1) = a2 * a2; R(0, 2) = a3 * a3;
    R(0, 3) = 2 * a1 * a2; R(0, 4) = 2 * a1 * a3; R(0, 5) = 2 * a2 * a3;
    R(1, 0) = b1 * b1; R(1, 1) = b2 * b2; R(1, 2) = b3 * b3;
    R(1, 3) = 2 * b1 * b2; R(1, 4) = 2 * b1 * b3; R(1, 5) = 2 * b2 * b3;
    R(2, 0) = a1 * b1; R(2, 1) = a2 * b2; R(2, 2) = a3 * b3;
    R(2, 3) = a1 * b2 + a2 * b1;
    R(2, 4) = a1 * b3 + a3 * b1;
    R(2, 5) = a2 * b3 + a3 * b2;
    return R;
}

// globally fixed unit normal of the edge with endpoints (lo,hi) -- depends only
// on the sorted endpoints, so adjacent elements get the identical functional.
Vector2d edgeNormal(const Vector2d& Plo, const Vector2d& Phi) {
    Vector2d t = Phi - Plo;
    double len = t.norm();
    return Vector2d(t(1), -t(0)) / len;
}

} // namespace

ArgyrisFEM::ArgyrisFEM(Mesh& m) : mesh(m) {
    gradbasis_my(mesh, Dlam, area);
    int NT = mesh.elem.rows();
    Rplain.resize(NT);
    for (int t = 0; t < NT; ++t) Rplain[t] = buildRplain(Dlam[t]);
    buildDOFMap();
    buildLocalBases();
}

void ArgyrisFEM::buildDOFMap() {
    nNode = mesh.node.rows();
    mesh.getEdge2Side(edge, edge2side);
    nEdge = edge.rows();
    nDof = 6 * nNode + nEdge;

    std::map<std::pair<int, int>, int> edgeMap;
    for (int e = 0; e < nEdge; ++e)
        edgeMap[{edge(e, 0), edge(e, 1)}] = e;

    int NT = mesh.elem.rows();
    elem2dof.resize(NT, locDof);
    for (int t = 0; t < NT; ++t) {
        for (int v = 0; v < 3; ++v) {
            int gn = mesh.elem(t, v);
            for (int c = 0; c < 6; ++c) elem2dof(t, 6 * v + c) = 6 * gn + c;
        }
        for (int e = 0; e < 3; ++e) {
            int na = mesh.elem(t, kEdgeNodes[e][0]);
            int nb = mesh.elem(t, kEdgeNodes[e][1]);
            int ge = edgeMap.at({std::min(na, nb), std::max(na, nb)});
            elem2dof(t, 18 + e) = 6 * nNode + ge;
        }
    }
}

void ArgyrisFEM::buildLocalBases() {
    std::array<Vector3d, 6> site;
    nodalSites(site);

    // geometry-independent raw P5 spans at the 6 nodal sites
    std::array<RowVectorXd, 6> valSpan;       // 1 x 21
    std::array<MatrixXd, 6>    gradLam;       // 3 x 21
    std::array<MatrixXd, 6>    hessLam;       // 6 x 21
    for (int k = 0; k < 6; ++k) {
        MatrixXd p(1, 3); p << site[k].transpose();
        valSpan[k] = polyBasisHomo3D(5, p).row(0);
        gradLam[k] = polyBasisHomoGrad3D(5, site[k]);
        hessLam[k] = polyBasisHomoHess3D(5, site[k]);
    }

    int NT = mesh.elem.rows();
    coef.resize(NT);
    for (int t = 0; t < NT; ++t) {
        const MatrixXd& Dt = Dlam[t];
        MatrixXd V(locDof, locDof);

        for (int v = 0; v < 3; ++v) {
            MatrixXd gradPhys = Dt * gradLam[v];          // 2 x 21
            MatrixXd Hphys    = Rplain[t] * hessLam[v];   // 3 x 21 = [Hxx;Hyy;Hxy]
            V.row(6 * v + 0) = valSpan[v];
            V.row(6 * v + 1) = gradPhys.row(0);
            V.row(6 * v + 2) = gradPhys.row(1);
            V.row(6 * v + 3) = Hphys.row(0);              // Hxx
            V.row(6 * v + 4) = Hphys.row(2);              // Hxy
            V.row(6 * v + 5) = Hphys.row(1);              // Hyy
        }
        for (int e = 0; e < 3; ++e) {
            int na = mesh.elem(t, kEdgeNodes[e][0]);
            int nb = mesh.elem(t, kEdgeNodes[e][1]);
            int lo = std::min(na, nb), hi = std::max(na, nb);
            Vector2d Plo = mesh.node.row(lo).transpose();
            Vector2d Phi = mesh.node.row(hi).transpose();
            Vector2d ne  = edgeNormal(Plo, Phi);
            MatrixXd gradPhys = Dt * gradLam[3 + e];      // 2 x 21
            V.row(18 + e) = ne(0) * gradPhys.row(0) + ne(1) * gradPhys.row(1);
        }

        coef[t] = V.fullPivLu().inverse();                // columns = nodal basis
    }
}

SparseMatrix<double> ArgyrisFEM::assembleStiffness() const {
    int NT = mesh.elem.rows();
    MatrixXd quadL; VectorXd w;
    Quadrature::quadpts2_my(kQuadOrder, quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<MatrixXd> Hlam_q(nq);
    for (int q = 0; q < nq; ++q)
        Hlam_q[q] = polyBasisHomoHess3D(5, Vector3d(quadL.row(q).transpose())); // 6 x 21

    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(NT) * locDof * locDof);

    const double s2 = std::sqrt(2.0);
    for (int t = 0; t < NT; ++t) {
        MatrixXd At = MatrixXd::Zero(locDof, locDof);
        for (int q = 0; q < nq; ++q) {
            MatrixXd Hn = (Rplain[t] * Hlam_q[q]) * coef[t];   // 3 x 21 = [Hxx;Hyy;Hxy]
            Hn.row(2) *= s2;                                   // -> [Hxx;Hyy;sqrt2 Hxy]
            At.noalias() += (w(q) * area(t)) * (Hn.transpose() * Hn);
        }
        for (int i = 0; i < locDof; ++i)
            for (int j = 0; j < locDof; ++j)
                trip.emplace_back(elem2dof(t, i), elem2dof(t, j), At(i, j));
    }
    SparseMatrix<double> A(nDof, nDof);
    A.setFromTriplets(trip.begin(), trip.end());
    return A;
}

VectorXd ArgyrisFEM::assembleLoad(const ExactSolution& sol) const {
    int NT = mesh.elem.rows();
    MatrixXd quadL; VectorXd w;
    Quadrature::quadpts2_my(kQuadOrder, quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<RowVectorXd> valSpan_q(nq);
    for (int q = 0; q < nq; ++q)
        valSpan_q[q] = polyBasisHomo3D(5, quadL.row(q)).row(0); // 1 x 21

    VectorXd F = VectorXd::Zero(nDof);
    for (int t = 0; t < NT; ++t) {
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)).transpose();
        Vector2d p2 = mesh.node.row(mesh.elem(t, 1)).transpose();
        Vector2d p3 = mesh.node.row(mesh.elem(t, 2)).transpose();
        VectorXd Ft = VectorXd::Zero(locDof);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            RowVectorXd phi = valSpan_q[q] * coef[t];          // 1 x 21
            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            MatrixXd ptMat(1, 2); ptMat << pt.transpose();
            double fval = sol.lap_lap_u_exact(ptMat)(0);
            Ft.noalias() += (w(q) * area(t) * fval) * phi.transpose();
        }
        for (int i = 0; i < locDof; ++i) F(elem2dof(t, i)) += Ft(i);
    }
    return F;
}

void ArgyrisFEM::strongBC(const ExactSolution& sol, VectorXd& c,
                          std::vector<int>& freeDof) const {
    c = VectorXd::Zero(nDof);
    std::vector<char> fixed(nDof, 0);

    // Clamped plate on an axis-aligned square.  For the test space v=0 and
    // d_n v=0 on dOmega, the C^1 trace argument requires pinning, at every
    // boundary vertex: value, both gradient components, the mixed second
    // derivative u_xy, and the second *tangential* derivative.  The only DOF
    // left FREE is the pure second *normal* derivative d_nn u (on a horizontal
    // edge that is u_yy; on a vertical edge u_xx).  Leaving exactly that DOF
    // free is what preserves the Aubin-Nitsche L2 superconvergence (rate k+1);
    // at corners two normals meet, so all 6 vertex DOFs are pinned.
    const double tol = 1e-10;
    auto onBdry = [&](double z) { return std::abs(z) < tol || std::abs(z - 1.0) < tol; };
    for (int v = 0; v < nNode; ++v) {
        double x = mesh.node(v, 0), y = mesh.node(v, 1);
        bool vert = onBdry(x);   // on a vertical edge   -> normal x -> d_nn = u_xx (idx 3)
        bool horz = onBdry(y);   // on a horizontal edge -> normal y -> d_nn = u_yy (idx 5)
        if (!vert && !horz) continue;

        MatrixXd P(1, 2); P << x, y;
        double u = sol.u_exact(P)(0);
        MatrixXd g = sol.grad_u_exact(P);       // 1 x 2
        MatrixXd H = sol.hessian_u_exact(P);    // 1 x 3 : [uxx,uyy,uxy]
        double comp[6] = {u, g(0, 0), g(0, 1), H(0, 0), H(0, 2), H(0, 1)};
        //               val  u_x      u_y      u_xx     u_xy     u_yy

        bool pin[6] = {true, true, true, true, true, true};
        if (vert && !horz) pin[3] = false;      // free u_xx (= d_nn) on left/right edge
        else if (horz && !vert) pin[5] = false; // free u_yy (= d_nn) on bottom/top edge
        // corner (vert && horz): pin everything

        for (int cc = 0; cc < 6; ++cc) {
            if (!pin[cc]) continue;
            c(6 * v + cc) = comp[cc];
            fixed[6 * v + cc] = 1;
        }
    }
    for (int e = 0; e < nEdge; ++e) {
        bool boundary = (edge2side(e, 0) == -1) || (edge2side(e, 1) == -1);
        if (!boundary) continue;
        int lo = edge(e, 0), hi = edge(e, 1);
        Vector2d Plo = mesh.node.row(lo).transpose();
        Vector2d Phi = mesh.node.row(hi).transpose();
        Vector2d mid = 0.5 * (Plo + Phi);
        Vector2d ne  = edgeNormal(Plo, Phi);
        MatrixXd midMat(1, 2); midMat << mid.transpose();
        MatrixXd g = sol.grad_u_exact(midMat);
        int gd = 6 * nNode + e;
        c(gd) = g(0, 0) * ne(0) + g(0, 1) * ne(1);
        fixed[gd] = 1;
    }

    freeDof.clear();
    for (int g = 0; g < nDof; ++g) if (!fixed[g]) freeDof.push_back(g);
}

void ArgyrisFEM::errors(const VectorXd& c, const ExactSolution& sol,
                        double& errL2, double& errH1, double& errH2) const {
    int NT = mesh.elem.rows();
    MatrixXd quadL; VectorXd w;
    Quadrature::quadpts2_my(kQuadOrder, quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<RowVectorXd> valSpan_q(nq);
    std::vector<MatrixXd>     gradLam_q(nq);
    std::vector<MatrixXd>     hessLam_q(nq);
    for (int q = 0; q < nq; ++q) {
        valSpan_q[q] = polyBasisHomo3D(5, quadL.row(q)).row(0);
        gradLam_q[q] = polyBasisHomoGrad3D(5, Vector3d(quadL.row(q).transpose()));
        hessLam_q[q] = polyBasisHomoHess3D(5, Vector3d(quadL.row(q).transpose()));
    }

    double sL2 = 0, sH1 = 0, sH2 = 0;
    for (int t = 0; t < NT; ++t) {
        VectorXd uh(locDof);
        for (int i = 0; i < locDof; ++i) uh(i) = c(elem2dof(t, i));
        const MatrixXd& Dt = Dlam[t];
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)).transpose();
        Vector2d p2 = mesh.node.row(mesh.elem(t, 1)).transpose();
        Vector2d p3 = mesh.node.row(mesh.elem(t, 2)).transpose();

        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            RowVectorXd phi = valSpan_q[q] * coef[t];                 // 1 x 21
            MatrixXd gradPhys = Dt * (gradLam_q[q] * coef[t]);        // 2 x 21
            MatrixXd Hn = (Rplain[t] * hessLam_q[q]) * coef[t];       // 3 x 21

            double uh_v = phi.dot(uh);
            Vector2d gh = gradPhys * uh;                              // (ux,uy)
            Vector3d Hh = Hn * uh;                                    // (Hxx,Hyy,Hxy)

            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            MatrixXd ptMat(1, 2); ptMat << pt.transpose();
            double u_ex = sol.u_exact(ptMat)(0);
            MatrixXd g_ex = sol.grad_u_exact(ptMat);                  // 1 x 2
            MatrixXd H_ex = sol.hessian_u_exact(ptMat);               // 1 x 3 [uxx,uyy,uxy]

            double wgt = w(q) * area(t);
            double du = uh_v - u_ex;
            sL2 += wgt * du * du;
            double dgx = gh(0) - g_ex(0, 0), dgy = gh(1) - g_ex(0, 1);
            sH1 += wgt * (dgx * dgx + dgy * dgy);
            double dxx = Hh(0) - H_ex(0, 0);
            double dyy = Hh(1) - H_ex(0, 1);
            double dxy = Hh(2) - H_ex(0, 2);
            sH2 += wgt * (dxx * dxx + dyy * dyy + 2.0 * dxy * dxy);
        }
    }
    errL2 = std::sqrt(sL2);
    errH1 = std::sqrt(sL2 + sH1);
    errH2 = std::sqrt(sH2);
}
