#include "EulerDG.h"
#include "Quadrature.h"
#include "utils.h"          // numSplit3 (monomial multi-indices) for fast field rendering

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <thread>

namespace euler {

namespace {
// Portable std::thread parallel-for over disjoint index ranges.  worker(lo,hi)
// processes [lo,hi); each call gets its own range so threads never touch the
// same output rows (the element loops below are written so each element writes
// only its own DG block).  Allocate per-thread temporaries INSIDE the worker.
template <class Worker>
void parallel_ranges(int n, const Worker& worker) {
    unsigned hw = std::thread::hardware_concurrency();
    unsigned nt = std::max(1u, std::min<unsigned>(hw ? hw : 1u, (unsigned)std::max(1, n)));
    if (nt <= 1 || n < 64) { worker(0, n); return; }
    int chunk = (n + (int)nt - 1) / (int)nt;
    std::vector<std::thread> ths;
    for (unsigned t = 0; t < nt; ++t) {
        int lo = (int)t * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        ths.emplace_back([lo, hi, &worker] { worker(lo, hi); });
    }
    for (auto& th : ths) th.join();
}
} // namespace

// ===========================================================================
// Edge / basis trace helpers (same reference-trace bookkeeping as the NS DG
// solver: edge type 0=(0,1), 1=(1,2), 2=(2,0); dir 0=forward, 1=reverse).
// ===========================================================================
namespace {

std::pair<int, int> edgeTypeDir(int a, int b) {
    if ((a == 0 && b == 1) || (a == 1 && b == 0)) return {0, (a > b) ? 1 : 0};
    if ((a == 1 && b == 2) || (a == 2 && b == 1)) return {1, (a > b) ? 1 : 0};
    return {2, (a == 0 && b == 2) ? 1 : 0};
}

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

} // namespace

// ===========================================================================
// Opaque precomputed quadrature / trace tables.
// ===========================================================================
struct EulerDG::Impl {
    // volume rule
    MatrixXd quadL;                       // nqv x 3 barycentric
    VectorXd wv;                          // nqv
    std::vector<RowVectorXd> phiV;        // [nqv] : 1 x locDof
    std::vector<MatrixXd>    dphiV;       // [nqv] : 3 x locDof (d/dlam)
    // edge rule
    MatrixXd quad1d;                      // nqe x 2 (lam along edge)
    VectorXd w1d;                         // nqe
    int nqe = 0;
    std::vector<std::vector<std::vector<RowVectorXd>>> ephi;   // [3][2][nqe] : 1 x locDof
    std::vector<std::vector<std::vector<MatrixXd>>>    edphi;  // [3][2][nqe] : 3 x locDof
    // degree-(k-1) L2 projection residual operator on the reference element:
    // given element nodal coeffs fe (locDof), Pdrop * fe gives the part of f that
    // lies OUTSIDE P_{k-1} (used by the Persson-Peraire sensor). Built once.
    MatrixXd Pdrop;                       // locDof x locDof  (in nodal coeff space)
    MatrixXd Mref;                        // locDof x locDof  unit-area mass (mass = area*Mref)
    // precomputed per-edge geometry (topology is fixed) -- used by assembleAV.
    struct EG { bool interior; int ta, tb, et_a, dir_a, et_b, dir_b, n1, n2; double nx, ny, he; };
    std::vector<EG> eg;
    // per-element edge connectivity for the ELEMENT-CENTRIC (parallel, race-free)
    // residual: for each of element t's 3 edges, this element's trace (et,dir),
    // the neighbour (-1 = boundary) and its trace, the outward normal from t, the
    // edge length, the global edge index (for the BC tag), and the edge nodes.
    struct EE { int nb, et, dir, et_nb, dir_nb, ei, n1, n2; double nx, ny, he; };
    std::vector<std::array<EE, 3>> ee;
};

// ===========================================================================
// Construction + precompute
// ===========================================================================
EulerDG::EulerDG(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
                 const MatrixXi& edge, const MatrixXi& edge2side,
                 const VectorXi& edgeTag, const EulerConfig& cfg)
    : fem_(fem), mesh_(mesh), elem2dof_(elem2dof), edge_(edge), edge2side_(edge2side),
      edgeTag_(edgeTag), cfg_(cfg),
      NT_(static_cast<int>(mesh.elem.rows())), NE_(static_cast<int>(edge.rows())),
      locDof_(fem.locDof), nDof_(static_cast<int>(elem2dof.maxCoeff()) + 1), nstep_(0) {
    precompute();
}

void EulerDG::precompute() {
    P_ = std::make_shared<Impl>();
    Impl& P = *P_;

    // ---- volume rule + reference basis values/derivatives ----
    fem_.quad2d(P.quadL, P.wv);
    int nqv = static_cast<int>(P.wv.size());
    P.phiV.resize(nqv); P.dphiV.resize(nqv);
    for (int q = 0; q < nqv; ++q) {
        P.phiV[q]  = fem_.computeBasisValue_all(P.quadL.row(q)).row(0);
        P.dphiV[q] = fem_.computeBasisDlam_all(P.quadL.row(q).transpose());
    }

    // ---- reference unit-area mass matrix (M_e = area_e * Mref) and its inverse ----
    P.Mref = MatrixXd::Zero(locDof_, locDof_);
    for (int q = 0; q < nqv; ++q)
        P.Mref.noalias() += P.wv(q) * (P.phiV[q].transpose() * P.phiV[q]);
    MrefInv_ = P.Mref.inverse();

    // ---- edge rule + reference trace tables ----
    fem_.quad1d(P.quad1d, P.w1d);
    P.nqe = static_cast<int>(P.w1d.size());
    P.ephi.assign(3, std::vector<std::vector<RowVectorXd>>(2, std::vector<RowVectorXd>(P.nqe)));
    P.edphi.assign(3, std::vector<std::vector<MatrixXd>>(2, std::vector<MatrixXd>(P.nqe)));
    for (int et = 0; et < 3; ++et) {
        int i1 = 0, i2 = 1;
        if (et == 1) { i1 = 1; i2 = 2; } else if (et == 2) { i1 = 2; i2 = 0; }
        for (int dir = 0; dir < 2; ++dir)
            for (int q = 0; q < P.nqe; ++q) {
                double l1 = P.quad1d(q, 0), l2 = P.quad1d(q, 1);
                Vector3d lam = Vector3d::Zero();
                if (dir == 0) { lam(i1) = l1; lam(i2) = l2; }
                else          { lam(i1) = l2; lam(i2) = l1; }
                P.ephi[et][dir][q]  = fem_.computeBasisValue_all(lam.transpose()).row(0);
                P.edphi[et][dir][q] = fem_.computeBasisDlam_all(lam);
            }
    }

    // ---- per-element size (h_K = longest edge, the Persson-Peraire AV scale) ----
    hK_.resize(NT_);
    for (int t = 0; t < NT_; ++t) {
        Vector2d p0 = mesh_.node.row(mesh_.elem(t, 0)), p1 = mesh_.node.row(mesh_.elem(t, 1)),
                 p2 = mesh_.node.row(mesh_.elem(t, 2));
        hK_(t) = std::max({(p1 - p0).norm(), (p2 - p1).norm(), (p0 - p2).norm()});
    }
    hmin_ = hK_.minCoeff();

    // ---- block-diagonal DG mass matrix (only the EULERCHK debug path uses it; the
    // production solve applies M element-wise via Mref).  Building this NT x locDof^2
    // sparse matrix every (re)construction is pure overhead under AMR, so skip it
    // unless the partitioned-solve verification is enabled. ----
    static const bool CHK = std::getenv("EULERCHK") != nullptr;
    if (CHK) {
        std::vector<Triplet<double>> trip; trip.reserve(static_cast<size_t>(NT_) * locDof_ * locDof_);
        for (int t = 0; t < NT_; ++t) {
            MatrixXd Me = fem_.area(t) * P.Mref;
            for (int i = 0; i < locDof_; ++i)
                for (int j = 0; j < locDof_; ++j)
                    trip.emplace_back(elem2dof_(t, i), elem2dof_(t, j), Me(i, j));
        }
        Mblk_.resize(nDof_, nDof_);
        Mblk_.setFromTriplets(trip.begin(), trip.end());
    }

    // ---- degree-(k-1) projection-residual operator Pdrop (for the AV sensor) ----
    // In nodal coeff space: project f onto P_{k-1} (lower monomial modes), then the
    // residual coeffs are (I - Proj) f. We build Proj via the monomial split: drop the
    // top-degree homogeneous monomials. Using the L2-mass inner product on the element.
    {
        int k = fem_.ord;
        if (k >= 1) {
            // Build the L2 projection onto span{ phi values restricted to degree<=k-1 }.
            // Work in the monomial (numSplit3) basis: span monomials of total degree<=k-1.
            // coef maps nodal coeffs -> monomial coeffs (mono = coef * fe in the renderer),
            // here we instead use the basis-value matrix at quad points to do a discrete L2
            // projection onto the lower-degree polynomial space spanned by the first
            // locDof_{k-1} nodal basis of a degree (k-1) element. Simpler/robust approach:
            // least-squares fit the degree-k nodal field by a degree-(k-1) polynomial in
            // the SAME quadrature, using the monomial basis up to degree k-1.
            MatrixXi midxK   = numSplit3(k);     // 3 x locDof   (degree-k homogeneous in lam)
            // homogeneous-degree-k monomials in barycentric lam are a basis of P_k restricted
            // to the triangle; the subset with exponent on lam1 (or any) ">=1 in all" is messy.
            // Instead evaluate a full P_{k-1} monomial set {lam2^a lam3^b : a+b<=k-1} and fit.
            std::vector<std::pair<int,int>> lowmono;
            for (int a = 0; a <= k - 1; ++a)
                for (int b = 0; b <= k - 1 - a; ++b) lowmono.emplace_back(a, b);
            int m = static_cast<int>(lowmono.size());
            int nqv2 = static_cast<int>(P.wv.size());
            // Vandermonde of low monomials at quad points (nqv2 x m) and the nodal basis (nqv2 x locDof)
            MatrixXd Bmono(nqv2, m), Phi(nqv2, locDof_);
            for (int q = 0; q < nqv2; ++q) {
                double l2 = P.quadL(q, 1), l3 = P.quadL(q, 2);
                for (int j = 0; j < m; ++j)
                    Bmono(q, j) = std::pow(l2, lowmono[j].first) * std::pow(l3, lowmono[j].second);
                Phi.row(q) = P.phiV[q];
            }
            // weighted normal equations:  (B^T W B) cmono = B^T W (Phi fe);  projected values = B cmono.
            VectorXd Wq = P.wv;
            MatrixXd BtWB = Bmono.transpose() * Wq.asDiagonal() * Bmono;   // m x m
            MatrixXd BtW  = Bmono.transpose() * Wq.asDiagonal();           // m x nqv2
            // Proj_values (nqv2 x locDof) = B (BtWB)^{-1} BtW Phi  ; residual values = Phi - Proj_values.
            MatrixXd projVals = Bmono * BtWB.ldlt().solve(BtW * Phi);      // nqv2 x locDof
            MatrixXd resVals  = Phi - projVals;                           // nqv2 x locDof
            // residual nodal->residual L2 energy operator:  Pdrop(i,j) = sum_q w res_i res_j
            // so that  fe^T Pdrop fe = ||(I-Proj) f||^2_{L2}/area.
            P.Pdrop = resVals.transpose() * Wq.asDiagonal() * resVals;     // locDof x locDof (SPD)
        } else {
            P.Pdrop = MatrixXd::Zero(locDof_, locDof_);
        }
    }

    // ---- precompute per-edge geometry (topology fixed) ----
    P.eg.resize(NE_);
    for (int e = 0; e < NE_; ++e) {
        int t1 = edge2side_(e, 0), t2 = edge2side_(e, 1);
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        bool interior = (t1 != -1 && t2 != -1);
        int ta = (t1 != -1) ? t1 : t2;
        EdgeOnElem ea = edgeOnElem(mesh_, ta, n1, n2);
        Impl::EG g{};
        g.interior = interior; g.ta = ta; g.tb = interior ? t2 : -1;
        g.et_a = ea.et; g.dir_a = ea.dir; g.nx = ea.nout.x(); g.ny = ea.nout.y(); g.he = ea.he;
        g.n1 = n1; g.n2 = n2;
        if (interior) { EdgeOnElem eb = edgeOnElem(mesh_, t2, n1, n2); g.et_b = eb.et; g.dir_b = eb.dir; }
        P.eg[e] = g;
    }

    // ---- per-element edge connectivity (for the element-centric parallel residual) ----
    {
        // (min,max node) -> global edge index, via a flat hash (cheaper than std::map
        // when EulerDG is reconstructed every remesh under AMR).
        auto pkey = [](int a, int b) { return ((int64_t)std::min(a, b) << 32) | (uint32_t)std::max(a, b); };
        std::unordered_map<int64_t, int> emap;
        emap.reserve(static_cast<size_t>(NE_) * 2);
        for (int e = 0; e < NE_; ++e) emap[pkey(edge_(e, 0), edge_(e, 1))] = e;
        P.ee.assign(NT_, std::array<Impl::EE, 3>{});
        const int loc[3][2] = {{0, 1}, {1, 2}, {2, 0}};
        for (int t = 0; t < NT_; ++t)
            for (int k = 0; k < 3; ++k) {
                int va = mesh_.elem(t, loc[k][0]), vb = mesh_.elem(t, loc[k][1]);
                int e = emap[pkey(va, vb)];
                int n1 = edge_(e, 0), n2 = edge_(e, 1);            // global edge orientation
                EdgeOnElem es = edgeOnElem(mesh_, t, n1, n2);      // THIS element's trace + outward normal
                int nb = (edge2side_(e, 0) == t) ? edge2side_(e, 1) : edge2side_(e, 0);
                Impl::EE r{};
                r.nb = nb; r.et = es.et; r.dir = es.dir; r.ei = e; r.n1 = n1; r.n2 = n2;
                r.nx = es.nout.x(); r.ny = es.nout.y(); r.he = es.he;
                if (nb != -1) { EdgeOnElem en = edgeOnElem(mesh_, nb, n1, n2); r.et_nb = en.et; r.dir_nb = en.dir; }
                P.ee[t][k] = r;
            }
    }

    epsFrozen_ = VectorXd::Zero(NT_);
    U_ = MatrixXd::Zero(nDof_, 4);
}

void EulerDG::setState(const MatrixXd& U0) {
    U_ = U0;
    // The L2 projection of a discontinuous initial condition (e.g. the DMR shock)
    // has Gibbs over/undershoots that can make rho or p negative; clamp them with
    // the positivity limiter so the very first flux evaluation is well-defined.
    if (cfg_.use_positivity) positivityLimit(U_);
    nstep_ = 0;
    implicitReady_ = false;
    refreshCountdown_ = 0;
}

// ===========================================================================
// Inviscid spatial residual  R = int F:grad phi  -  oint Hhat phi   (nDof x 4)
//   dU/dt|_explicit = M^{-1} R
// ===========================================================================
void EulerDG::inviscidResidual(const MatrixXd& U, double t, const ExteriorStateFn& bc,
                               MatrixXd& R) const {
    const Impl& P = *P_;
    R.setZero(nDof_, 4);
    const int nqv = static_cast<int>(P.wv.size());
    const bool hllcFlux = cfg_.use_hllc;

    // ELEMENT-CENTRIC residual: each element computes its own volume term + its 3
    // face fluxes (reading neighbours / ghost states) and writes ONLY its own DG
    // block, so the loop parallelises over elements with no write races.  Interior
    // face fluxes are evaluated twice (once per adjacent element) -- cheap vs the
    // ~10x core speedup.
    parallel_ranges(NT_, [&](int lo, int hi) {
        Vector4d Fx, Fy;
        MatrixXd Ue(locDof_, 4), Unb(locDof_, 4), G(2, locDof_);
        Matrix<double, Dynamic, 4> Re(locDof_, 4);
        for (int tt = lo; tt < hi; ++tt) {
            for (int i = 0; i < locDof_; ++i) Ue.row(i) = U.row(elem2dof_(tt, i));
            Re.setZero();
            double area = fem_.area(tt);
            // volume:  + int_K (Fx d_x phi_i + Fy d_y phi_i)
            for (int q = 0; q < nqv; ++q) {
                Vector4d Uq = (P.phiV[q] * Ue).transpose();
                fluxes(Uq, Fx, Fy);
                G.noalias() = fem_.Dlam[tt] * P.dphiV[q];
                double wa = P.wv(q) * area;
                Re.noalias() += (wa * G.row(0).transpose()) * Fx.transpose();
                Re.noalias() += (wa * G.row(1).transpose()) * Fy.transpose();
            }
            // faces:  - oint Hhat phi   over this element's 3 edges
            for (int k = 0; k < 3; ++k) {
                const Impl::EE& r = P.ee[tt][k];
                bool interior = (r.nb != -1);
                if (interior) for (int i = 0; i < locDof_; ++i) Unb.row(i) = U.row(elem2dof_(r.nb, i));
                int tag = edgeTag_.size() ? edgeTag_(r.ei) : 0;
                for (int q = 0; q < P.nqe; ++q) {
                    const RowVectorXd& pa = P.ephi[r.et][r.dir][q];
                    Vector4d Um = (pa * Ue).transpose();
                    Vector4d Up;
                    if (interior) {
                        Up = (P.ephi[r.et_nb][r.dir_nb][q] * Unb).transpose();
                    } else {
                        double l1 = P.quad1d(q, 0), l2 = P.quad1d(q, 1);
                        Vector2d Pp = l1 * mesh_.node.row(r.n1).transpose() + l2 * mesh_.node.row(r.n2).transpose();
                        Up = bc ? bc(Pp.x(), Pp.y(), t, Um, r.nx, r.ny, tag) : Um;
                    }
                    Vector4d Hn = hllcFlux ? hllc(Um, Up, r.nx, r.ny) : rusanov(Um, Up, r.nx, r.ny);
                    double whe = P.w1d(q) * r.he;
                    Re.noalias() -= (whe * pa.transpose()) * Hn.transpose();
                }
            }
            for (int i = 0; i < locDof_; ++i) R.row(elem2dof_(tt, i)) = Re.row(i);
        }
    });
}

void EulerDG::applyMassInverse(MatrixXd& R) const {
    parallel_ranges(NT_, [&](int lo, int hi) {
      Matrix<double, Dynamic, 4> Re(locDof_, 4);
      for (int t = lo; t < hi; ++t) {
        for (int i = 0; i < locDof_; ++i) Re.row(i) = R.row(elem2dof_(t, i));
        Re = (MrefInv_ / fem_.area(t)) * Re;
        for (int i = 0; i < locDof_; ++i) R.row(elem2dof_(t, i)) = Re.row(i);
      }
    });
}

// ===========================================================================
// Diagnostics
// ===========================================================================
double EulerDG::maxWaveSpeedGlobal() const {
    double m = 0.0;
    for (int i = 0; i < nDof_; ++i) {
        Vector4d U = U_.row(i).transpose();
        if (U(0) <= 0) continue;
        double sp = std::hypot(U(1), U(2)) / U(0) + soundSpeed(U);
        m = std::max(m, sp);
    }
    return m;
}

void EulerDG::quadMin(double& rmin, double& pmin, int& worstElem, double& cx, double& cy) const {
    const Impl& P = *P_;
    int nqv = static_cast<int>(P.wv.size());
    rmin = 1e300; pmin = 1e300; worstElem = -1;
    for (int t = 0; t < NT_; ++t) {
        MatrixXd Ue(locDof_, 4);
        for (int i = 0; i < locDof_; ++i) Ue.row(i) = U_.row(elem2dof_(t, i));
        auto check = [&](const RowVectorXd& ph) {
            Vector4d Uq = (ph * Ue).transpose();
            double r = Uq(0), p = pressure(Uq);
            if (r < rmin) rmin = r;
            if (p < pmin) { pmin = p; }
            if (r < 1e-12 || p < 1e-12) {
                worstElem = t;
                Vector2d c = (mesh_.node.row(mesh_.elem(t,0)) + mesh_.node.row(mesh_.elem(t,1)) + mesh_.node.row(mesh_.elem(t,2))).transpose() / 3.0;
                cx = c.x(); cy = c.y();
            }
        };
        for (int q = 0; q < nqv; ++q) check(P.phiV[q]);
        for (int et = 0; et < 3; ++et) for (int dir = 0; dir < 2; ++dir)
            for (int q = 0; q < P.nqe; ++q) check(P.ephi[et][dir][q]);
    }
}

Vector4d EulerDG::conservedTotals() const {
    const Impl& P = *P_;
    int nqv = static_cast<int>(P.wv.size());
    Vector4d tot = Vector4d::Zero();
    for (int t = 0; t < NT_; ++t) {
        double area = fem_.area(t);
        for (int q = 0; q < nqv; ++q) {
            Vector4d Uq = Vector4d::Zero();
            for (int i = 0; i < locDof_; ++i) Uq += P.phiV[q](i) * U_.row(elem2dof_(t, i)).transpose();
            tot += (P.wv(q) * area) * Uq;
        }
    }
    return tot;
}

VectorXd EulerDG::pressureField() const {
    VectorXd f(nDof_);
    for (int i = 0; i < nDof_; ++i) f(i) = pressure(U_.row(i).transpose(), cfg_.rho_floor, cfg_.p_floor);
    return f;
}
VectorXd EulerDG::machField() const {
    VectorXd f(nDof_);
    for (int i = 0; i < nDof_; ++i) {
        Vector4d U = U_.row(i).transpose();
        double rho = std::max(U(0), cfg_.rho_floor);
        double sp  = std::hypot(U(1), U(2)) / rho;
        f(i) = sp / std::max(soundSpeed(U), 1e-300);
    }
    return f;
}
VectorXd EulerDG::velocityMagField() const {
    VectorXd f(nDof_);
    for (int i = 0; i < nDof_; ++i) {
        Vector4d U = U_.row(i).transpose();
        f(i) = std::hypot(U(1), U(2)) / std::max(U(0), cfg_.rho_floor);
    }
    return f;
}

VectorXd EulerDG::schlierenField() const {
    // numerical Schlieren:  exp(-beta |grad rho| / max|grad rho|), element-wise constant
    // gradient (cheap, dramatic for shock visualisation).  Broadcast to nodes.
    VectorXd gmag(NT_);
    for (int t = 0; t < NT_; ++t) {
        // least-squares element gradient of rho from nodal coeffs at the centroid
        Vector3d lamc(1.0/3, 1.0/3, 1.0/3);
        MatrixXd G = fem_.computeBasisGrad_all(t, lamc);   // 2 x locDof
        double gx = 0, gy = 0;
        for (int i = 0; i < locDof_; ++i) { gx += G(0,i) * U_(elem2dof_(t,i),0); gy += G(1,i) * U_(elem2dof_(t,i),0); }
        gmag(t) = std::hypot(gx, gy);
    }
    double gmax = gmag.maxCoeff(); if (gmax < 1e-30) gmax = 1.0;
    VectorXd f(nDof_);
    for (int t = 0; t < NT_; ++t) {
        double s = std::exp(-30.0 * gmag(t) / gmax);   // bright = smooth, dark = shock
        for (int i = 0; i < locDof_; ++i) f(elem2dof_(t, i)) = s;
    }
    return f;
}

VectorXd EulerDG::avField() const {
    VectorXd f(nDof_);
    for (int t = 0; t < NT_; ++t)
        for (int i = 0; i < locDof_; ++i) f(elem2dof_(t, i)) = (t < epsFrozen_.size()) ? epsFrozen_(t) : 0.0;
    return f;
}

Vector4d EulerDG::cellAverage(const MatrixXd& U, int t) const {
    const Impl& P = *P_;
    int nqv = static_cast<int>(P.wv.size());
    Vector4d avg = Vector4d::Zero();
    for (int q = 0; q < nqv; ++q) {
        Vector4d Uq = Vector4d::Zero();
        for (int i = 0; i < locDof_; ++i) Uq += P.phiV[q](i) * U.row(elem2dof_(t, i)).transpose();
        avg += P.wv(q) * Uq;     // weights sum to 1 -> average over the element
    }
    return avg;
}

void EulerDG::elementBlock(const MatrixXd& U, int t, Matrix<double, Dynamic, 4>& Ue) const {
    Ue.resize(locDof_, 4);
    for (int i = 0; i < locDof_; ++i) Ue.row(i) = U.row(elem2dof_(t, i));
}

// ===========================================================================
// L2 projection of an initial primitive field + L2 error  (free functions)
// ===========================================================================
void makeRectMesh(Mesh& mesh, double x0, double x1, double y0, double y1, int nx, int ny) {
    int nnx = nx + 1, nny = ny + 1;
    mesh.node.resize(nnx * nny, 2);
    double dx = (x1 - x0) / nx, dy = (y1 - y0) / ny;
    auto nid = [&](int i, int j){ return j * nnx + i; };
    for (int j = 0; j < nny; ++j)
        for (int i = 0; i < nnx; ++i) {
            mesh.node(nid(i, j), 0) = x0 + i * dx;
            mesh.node(nid(i, j), 1) = y0 + j * dy;
        }
    mesh.elem.resize(2 * nx * ny, 3);
    int t = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int a = nid(i, j), b = nid(i + 1, j), c = nid(i + 1, j + 1), d = nid(i, j + 1);
            mesh.elem.row(t++) << a, b, c;      // CCW
            mesh.elem.row(t++) << a, c, d;      // CCW
        }
}

MatrixXd projectInitial(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                        const std::function<Vector4d(double, double)>& primField) {
    int NT = static_cast<int>(mesh.elem.rows()), locDof = fem.locDof;
    int nDof = static_cast<int>(elem2dof.maxCoeff()) + 1;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    MatrixXd Mref = MatrixXd::Zero(locDof, locDof);
    for (int q = 0; q < nq; ++q) { phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0); Mref.noalias() += w(q) * (phi[q].transpose() * phi[q]); }
    MatrixXd MrefInv = Mref.inverse();
    MatrixXd U = MatrixXd::Zero(nDof, 4);
    for (int t = 0; t < NT; ++t) {
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)), p2 = mesh.node.row(mesh.elem(t, 1)), p3 = mesh.node.row(mesh.elem(t, 2));
        Matrix<double, Dynamic, 4> be = Matrix<double, Dynamic, 4>::Zero(locDof, 4);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d Pp = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            Vector4d pr = primField(Pp.x(), Pp.y());
            Vector4d Uc = primToCons(pr(0), pr(1), pr(2), pr(3));
            be.noalias() += w(q) * (phi[q].transpose() * Uc.transpose());   // (/area later cancels)
        }
        Matrix<double, Dynamic, 4> Ue = MrefInv * be;   // area cancels (be has no area, Mref unit-area)
        for (int i = 0; i < locDof; ++i) U.row(elem2dof(t, i)) = Ue.row(i);
    }
    return U;
}

double l2Error(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
               const MatrixXd& U, int comp,
               const std::function<Vector4d(double, double)>& exactCons) {
    int NT = static_cast<int>(mesh.elem.rows()), locDof = fem.locDof;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
    double err = 0.0;
    for (int t = 0; t < NT; ++t) {
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)), p2 = mesh.node.row(mesh.elem(t, 1)), p3 = mesh.node.row(mesh.elem(t, 2));
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d Pp = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            double uh = 0.0;
            for (int i = 0; i < locDof; ++i) uh += phi[q](i) * U(elem2dof(t, i), comp);
            double ue = exactCons(Pp.x(), Pp.y())(comp);
            err += w(q) * area * (uh - ue) * (uh - ue);
        }
    }
    return std::sqrt(err);
}

// ===========================================================================
// Visualisation (true dP_k per-pixel sampling, multiple colormaps)
// ===========================================================================
namespace {
void cmap(double t, Colormap cm, unsigned char& R, unsigned char& G, unsigned char& B) {
    t = std::min(1.0, std::max(0.0, t));
    double col[3];
    auto lerp = [](const double* a, const double* b, double s, double* o){ for (int k=0;k<3;++k) o[k]=a[k]+s*(b[k]-a[k]); };
    if (cm == CM_GRAY)      { col[0]=col[1]=col[2]=255.0*t; }
    else if (cm == CM_GRAY_INV) { col[0]=col[1]=col[2]=255.0*(1.0-t); }
    else if (cm == CM_COOLWARM) {
        const double lo[3]={59,76,192}, mid[3]={242,242,242}, hi[3]={180,4,38};
        if (t<0.5) lerp(lo,mid,t/0.5,col); else lerp(mid,hi,(t-0.5)/0.5,col);
    } else if (cm == CM_JET) {
        const double c0[3]={0,0,131},c1[3]={0,60,170},c2[3]={5,255,255},c3[3]={255,255,0},c4[3]={250,0,0};
        const double* C[5]={c0,c1,c2,c3,c4}; double s=t*4.0; int k=std::min(3,(int)s); lerp(C[k],C[k+1],s-k,col);
    } else if (cm == CM_INFERNO) {
        const double c0[3]={0,0,4},c1[3]={87,16,110},c2[3]={188,55,84},c3[3]={249,142,9},c4[3]={252,255,164};
        const double* C[5]={c0,c1,c2,c3,c4}; double s=t*4.0; int k=std::min(3,(int)s); lerp(C[k],C[k+1],s-k,col);
    } else { // viridis
        const double c0[3]={68,1,84},c1[3]={59,82,139},c2[3]={33,145,140},c3[3]={94,201,98},c4[3]={253,231,37};
        const double* C[5]={c0,c1,c2,c3,c4}; double s=t*4.0; int k=std::min(3,(int)s); lerp(C[k],C[k+1],s-k,col);
    }
    R=(unsigned char)std::lround(col[0]); G=(unsigned char)std::lround(col[1]); B=(unsigned char)std::lround(col[2]);
}
}

void writeScalarPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                    const MatrixXi& elem2dof, const VectorXd& field, int W, int H,
                    double xmin, double xmax, double ymin, double ymax,
                    double vmin, double vmax, Colormap cm,
                    const std::function<bool(double, double)>& inDomain) {
    int locDof = fem.locDof;
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    std::vector<unsigned char> img((size_t)W * H * 3, 255);
    int NT = static_cast<int>(mesh.elem.rows());
    MatrixXi midx = numSplit3(fem.ord);
    VectorXd fe(locDof), mco(locDof);
    for (int t = 0; t < NT; ++t) {
        Vector2d P1 = mesh.node.row(mesh.elem(t,0)), P2 = mesh.node.row(mesh.elem(t,1)), P3 = mesh.node.row(mesh.elem(t,2));
        for (int i = 0; i < locDof; ++i) fe(i) = field(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        Vector2d v0 = P2 - P1, v1 = P3 - P1;
        double det = v0.x()*v1.y() - v1.x()*v0.y();
        if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;
        double txmin=std::min({P1.x(),P2.x(),P3.x()}), txmax=std::max({P1.x(),P2.x(),P3.x()});
        double tymin=std::min({P1.y(),P2.y(),P3.y()}), tymax=std::max({P1.y(),P2.y(),P3.y()});
        int col0=std::max(0,(int)std::floor((txmin-xmin)/dx)-1), col1=std::min(W-1,(int)std::ceil((txmax-xmin)/dx)+1);
        int row0=std::max(0,(int)std::floor((ymax-tymax)/dy)-1), row1=std::min(H-1,(int)std::ceil((ymax-tymin)/dy)+1);
        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double x = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double y = ymax - (row + 0.5) * dy;
                double px = x - P1.x(), py = y - P1.y();
                double l2 = (px*v1.y() - v1.x()*py)*invd, l3 = (v0.x()*py - px*v0.y())*invd, l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                if (inDomain && !inDomain(x, y)) continue;
                double val = 0.0;
                for (int j = 0; j < locDof; ++j) {
                    double s = mco(j);
                    for (int e = 0; e < midx(0,j); ++e) s *= l1;
                    for (int e = 0; e < midx(1,j); ++e) s *= l2;
                    for (int e = 0; e < midx(2,j); ++e) s *= l3;
                    val += s;
                }
                double s = (val - vmin) / (vmax - vmin);
                unsigned char R,G,B; cmap(s, cm, R, G, B);
                size_t idx = ((size_t)row*W + col)*3;
                img[idx]=R; img[idx+1]=G; img[idx+2]=B;
            }
        }
    }
    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "writeScalarPPM: cannot open " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

void writeSchlierenPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                       const MatrixXi& elem2dof, const VectorXd& rho, int W, int H,
                       double xmin, double xmax, double ymin, double ymax, double beta) {
    int locDof = fem.locDof, NT = static_cast<int>(mesh.elem.rows());
    MatrixXi midx = numSplit3(fem.ord);                       // 3 x locDof exponents
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    auto inWin = [&](int t) {
        Vector2d c = (mesh.node.row(mesh.elem(t,0)) + mesh.node.row(mesh.elem(t,1)) + mesh.node.row(mesh.elem(t,2))).transpose() / 3.0;
        double m = 0.05 * std::max(xmax - xmin, ymax - ymin);
        return c.x() > xmin - m && c.x() < xmax + m && c.y() > ymin - m && c.y() < ymax + m;
    };
    // gradient of the density polynomial at barycentric (l1,l2,l3) for element t
    VectorXd fe(locDof), mco(locDof);
    int ord = fem.ord;
    // gradient of the density polynomial via integer-power TABLES (no std::pow in
    // the per-pixel hot loop): pw1[e]=l1^e etc., built once per pixel.
    auto gradAt = [&](const MatrixXd& Dl, double l1, double l2, double l3, double& gx, double& gy) {
        double pw1[16], pw2[16], pw3[16];
        pw1[0] = pw2[0] = pw3[0] = 1.0;
        for (int e = 1; e <= ord; ++e) { pw1[e] = pw1[e-1]*l1; pw2[e] = pw2[e-1]*l2; pw3[e] = pw3[e-1]*l3; }
        double d1 = 0, d2 = 0, d3 = 0;
        for (int j = 0; j < locDof; ++j) {
            int a = midx(0, j), b = midx(1, j), c = midx(2, j);
            if (a > 0) d1 += mco(j) * a * pw1[a-1] * pw2[b]   * pw3[c];
            if (b > 0) d2 += mco(j) * b * pw1[a]   * pw2[b-1] * pw3[c];
            if (c > 0) d3 += mco(j) * c * pw1[a]   * pw2[b]   * pw3[c-1];
        }
        gx = Dl(0,0)*d1 + Dl(0,1)*d2 + Dl(0,2)*d3;
        gy = Dl(1,0)*d1 + Dl(1,1)*d2 + Dl(1,2)*d3;
    };
    // pass 1: max |grad rho| over the window (centroid sample, cheap)
    double gmax = 1e-30;
    Vector3d lc(1.0/3, 1.0/3, 1.0/3);
    for (int t = 0; t < NT; ++t) {
        if (!inWin(t)) continue;
        for (int i = 0; i < locDof; ++i) fe(i) = rho(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        double gx, gy; gradAt(fem.Dlam[t], lc(0), lc(1), lc(2), gx, gy);
        gmax = std::max(gmax, std::hypot(gx, gy));
    }
    // pass 2: per-pixel |grad rho| via the true dP_k polynomial gradient
    std::vector<unsigned char> img((size_t)W * H * 3, 255);
    for (int t = 0; t < NT; ++t) {
        if (!inWin(t)) continue;
        Vector2d P1 = mesh.node.row(mesh.elem(t,0)), P2 = mesh.node.row(mesh.elem(t,1)), P3 = mesh.node.row(mesh.elem(t,2));
        for (int i = 0; i < locDof; ++i) fe(i) = rho(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        const MatrixXd& Dl = fem.Dlam[t];
        Vector2d v0 = P2 - P1, v1 = P3 - P1;
        double det = v0.x()*v1.y() - v1.x()*v0.y(); if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;
        double txmin=std::min({P1.x(),P2.x(),P3.x()}), txmax=std::max({P1.x(),P2.x(),P3.x()});
        double tymin=std::min({P1.y(),P2.y(),P3.y()}), tymax=std::max({P1.y(),P2.y(),P3.y()});
        int col0=std::max(0,(int)std::floor((txmin-xmin)/dx)-1), col1=std::min(W-1,(int)std::ceil((txmax-xmin)/dx)+1);
        int row0=std::max(0,(int)std::floor((ymax-tymax)/dy)-1), row1=std::min(H-1,(int)std::ceil((ymax-tymin)/dy)+1);
        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double x = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double y = ymax - (row + 0.5) * dy;
                double px = x - P1.x(), py = y - P1.y();
                double l2 = (px*v1.y() - v1.x()*py)*invd, l3 = (v0.x()*py - px*v0.y())*invd, l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                double gx, gy; gradAt(Dl, l1, l2, l3, gx, gy);
                double s = std::exp(-beta * std::hypot(gx, gy) / gmax);   // white smooth, dark steep
                unsigned char v = (unsigned char)std::lround(255.0 * std::min(1.0, std::max(0.0, s)));
                size_t idx = ((size_t)row*W + col)*3;
                img[idx]=v; img[idx+1]=v; img[idx+2]=v;
            }
        }
    }
    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "writeSchlierenPPM: cannot open " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

// ===========================================================================
// Persson-Peraire artificial-viscosity sensor  ->  per-element eps_K
//   S_K = ||u - Pi_{k-1} u||^2 / ||u||^2 (density indicator),  s_K = log10 S_K,
//   eps0 = c_AV (h_K/k) lambda_max(K),  smooth sinusoidal ramp about s0.
// ===========================================================================
void EulerDG::computeViscosity(const MatrixXd& U, VectorXd& epsK) const {
    const Impl& P = *P_;
    epsK.resize(NT_);
    int k = fem_.ord;
    double s0 = cfg_.av_s0_set ? cfg_.av_s0 : ((k >= 2) ? -4.0 * std::log10((double)k) : -3.0);
    double kappa = cfg_.av_kappa;
    if (k < 1) { epsK.setZero(); return; }
    for (int t = 0; t < NT_; ++t) {
        VectorXd u(locDof_);                                  // smoothness indicator (nodal)
        for (int i = 0; i < locDof_; ++i) {
            if (cfg_.av_indicator == 0) {
                u(i) = U(elem2dof_(t, i), 0);                 // density (senses contacts too)
            } else {
                Vector4d Ui = U.row(elem2dof_(t, i)).transpose();
                double rho = std::max(Ui(0), 1e-10);
                u(i) = (GAMMA - 1.0) * (Ui(3) - 0.5 * (Ui(1) * Ui(1) + Ui(2) * Ui(2)) / rho); // pressure
            }
        }
        double num = u.dot(P.Pdrop * u);     // ||(I-Pi_{k-1}) u||^2 / area
        double den = u.dot(P.Mref  * u);     // ||u||^2 / area
        double eps = 0.0;
        if (den > 1e-30 && num > 0.0) {
            double s = std::log10(num / den);
            double lammax = 0.0;             // max (|u|+c) over the element nodes
            for (int i = 0; i < locDof_; ++i) {
                Vector4d Ui = U.row(elem2dof_(t, i)).transpose();
                if (Ui(0) <= 0) continue;
                double sp = std::hypot(Ui(1), Ui(2)) / Ui(0) + soundSpeed(Ui);
                lammax = std::max(lammax, sp);
            }
            double eps0 = cfg_.av_c * (hK_(t) / k) * lammax;
            if (s >= s0 + kappa)            eps = eps0;
            else if (s > s0 - kappa)        eps = 0.5 * eps0 * (1.0 + std::sin(M_PI * (s - s0) / (2.0 * kappa)));
            // else eps stays 0
        }
        epsK(t) = eps;
    }
}

// ===========================================================================
// Variable-coefficient SIPG discretisation of  -div(eps grad U), eps_K const.
// Symmetric ({eps}-weighted) interior-penalty form (SPSD); boundary faces are
// do-nothing so the artificial viscosity never pollutes the inviscid BCs.
// Acts identically on each conserved component, so one scalar operator suffices.
// ===========================================================================
SparseMatrix<double> EulerDG::assembleAV(const VectorXd& epsK) const {
    const Impl& P = *P_;
    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(NT_) * locDof_ * locDof_ * 2);

    // ---- volume:  + eps_K \int_K grad U . grad phi ----
    int nqv = static_cast<int>(P.wv.size());
    for (int t = 0; t < NT_; ++t) {
        double e = epsK(t);
        if (e <= 0.0) continue;
        double area = fem_.area(t);
        MatrixXd At = MatrixXd::Zero(locDof_, locDof_);
        for (int q = 0; q < nqv; ++q) {
            MatrixXd G = fem_.Dlam[t] * P.dphiV[q];            // 2 x locDof
            At.noalias() += (P.wv(q) * area) * (G.transpose() * G);
        }
        At *= e;
        for (int i = 0; i < locDof_; ++i)
            for (int j = 0; j < locDof_; ++j)
                trip.emplace_back(elem2dof_(t, i), elem2dof_(t, j), At(i, j));
    }

    // ---- interior faces:  -{eps gradU.n}[phi] - {eps gradphi.n}[U] + (sigma/he){eps}[U][phi] ----
    // The consistency term uses the CONSISTENT one-sided average {eps grad.n} =
    // 1/2(eps_K gradU_K.n + eps_K' gradU_K'.n), i.e. each side's normal-gradient is
    // weighted by ITS OWN eps -- NOT by {eps} applied to the whole block.  When a
    // neighbour has eps=0 (smooth side of a shock) this zeroes that side's gradient
    // term, so it is bounded by the (eps-weighted) volume energy and the operator is
    // genuinely SPSD; {eps}-weighting the whole block instead leaves an unbounded
    // gradient cross-term that no penalty can fix (would make M+gamma dt A indefinite
    // at larger dt).
    double sigma = cfg_.sigma_ip * (fem_.ord + 1) * (fem_.ord + 1);
    int L = locDof_;
    for (int e = 0; e < NE_; ++e) {
        int t1 = edge2side_(e, 0), t2 = edge2side_(e, 1);
        if (t1 == -1 || t2 == -1) continue;                    // do-nothing on boundary
        double e1v = epsK(t1), e2v = epsK(t2);
        double eps_avg = 0.5 * (e1v + e2v);
        if (eps_avg <= 0.0) continue;
        int n1 = edge_(e, 0), n2 = edge_(e, 1);
        EdgeOnElem e1 = edgeOnElem(mesh_, t1, n1, n2);         // normal n out of t1
        EdgeOnElem e2 = edgeOnElem(mesh_, t2, n1, n2);
        Vector2d nrm = e1.nout; double he = e1.he;
        double pen = sigma / he;
        MatrixXd Mloc = MatrixXd::Zero(2 * L, 2 * L);           // [side1 | side2]
        for (int q = 0; q < P.nqe; ++q) {
            RowVectorXd phi1 = P.ephi[e1.et][e1.dir][q];        // 1 x L
            RowVectorXd phi2 = P.ephi[e2.et][e2.dir][q];
            MatrixXd g1 = fem_.Dlam[t1] * P.edphi[e1.et][e1.dir][q];   // 2 x L
            MatrixXd g2 = fem_.Dlam[t2] * P.edphi[e2.et][e2.dir][q];
            RowVectorXd D1 = nrm.transpose() * g1;             // grad.n on side1
            RowVectorXd D2 = nrm.transpose() * g2;
            RowVectorXd J(2 * L), Aw(2 * L);
            J << phi1, -phi2;                                  // jump  [.] = side1 - side2
            Aw << 0.5 * e1v * D1, 0.5 * e2v * D2;              // {eps grad.n} (per-side eps)
            double whe = P.w1d(q) * he;
            Mloc.noalias() += whe *
                (-(Aw.transpose() * J) - (J.transpose() * Aw) + (eps_avg * pen) * (J.transpose() * J));
        }
        auto gdof = [&](int s, int i){ return (s == 0) ? elem2dof_(t1, i) : elem2dof_(t2, i); };
        for (int a = 0; a < 2 * L; ++a)
            for (int b = 0; b < 2 * L; ++b)
                trip.emplace_back(gdof(a / L, a % L), gdof(b / L, b % L), Mloc(a, b));
    }

    SparseMatrix<double> A(nDof_, nDof_);
    A.setFromTriplets(trip.begin(), trip.end());
    return A;
}

// MU = M U  (block-diagonal mass: M_e = area_e * Mref), per element, parallel.
void EulerDG::applyMass(const MatrixXd& U, MatrixXd& MU) const {
    MU.resize(nDof_, 4);
    parallel_ranges(NT_, [&](int lo, int hi) {
        Matrix<double, Dynamic, 4> Ue(locDof_, 4);
        for (int t = lo; t < hi; ++t) {
            for (int i = 0; i < locDof_; ++i) Ue.row(i) = U.row(elem2dof_(t, i));
            Ue = (fem_.area(t) * P_->Mref) * Ue;
            for (int i = 0; i < locDof_; ++i) MU.row(elem2dof_(t, i)) = Ue.row(i);
        }
    });
}

void EulerDG::refreshImplicit(const VectorXd& epsK, double a_dt) {
    A_ = assembleAV(epsK);                       // sparse artificial-viscosity operator
    aDtCached_ = a_dt;
    // "active" dofs = those touched by A (shock cells + their face-neighbours).
    // For every other dof, K = M (block diagonal), so K is EXACTLY block-diagonal
    // active/decoupled and the partitioned solve below is exact, not approximate.
    std::vector<char> isActive(nDof_, 0);
    for (int k = 0; k < A_.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(A_, k); it; ++it) {
            isActive[it.row()] = 1; isActive[it.col()] = 1;
        }
    activeDofs_.clear();
    std::vector<int> g2a(nDof_, -1);
    for (int i = 0; i < nDof_; ++i) if (isActive[i]) { g2a[i] = (int)activeDofs_.size(); activeDofs_.push_back(i); }
    int na = (int)activeDofs_.size();
    implicitReady_ = true;
    if (na == 0) { Aaa_.resize(0, 0); return; }   // smooth flow: K = M everywhere

    // K_aa = (M + a_dt A) restricted to active dofs.  M restricted to active = the
    // block-mass of active cells (an active cell has all its dofs active).  Also build
    // Aaa_ = A in active indexing (LTS uses it for the compact A*U2; every nonzero
    // column of A is active, so Aaa_ captures A's full action on the active dofs).
    std::vector<Triplet<double>> trip; trip.reserve((size_t)na * 8);
    std::vector<Triplet<double>> atrip; atrip.reserve((size_t)na * 8);
    for (int k = 0; k < A_.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(A_, k); it; ++it) {
            trip.emplace_back(g2a[it.row()], g2a[it.col()], a_dt * it.value());
            atrip.emplace_back(g2a[it.row()], g2a[it.col()], it.value());
        }
    Aaa_.resize(na, na);
    Aaa_.setFromTriplets(atrip.begin(), atrip.end());
    // add the M block of every active cell (all dofs of such a cell are active)
    std::vector<char> cellActive(NT_, 0);
    for (int i = 0; i < na; ++i) cellActive[activeDofs_[i] / locDof_] = 1;
    for (int t = 0; t < NT_; ++t) {
        if (!cellActive[t]) continue;
        MatrixXd Me = fem_.area(t) * P_->Mref;
        for (int i = 0; i < locDof_; ++i)
            for (int j = 0; j < locDof_; ++j)
                trip.emplace_back(g2a[elem2dof_(t, i)], g2a[elem2dof_(t, j)], Me(i, j));
    }
    Kaa_.resize(na, na);
    Kaa_.setFromTriplets(trip.begin(), trip.end());
    ldltA_.compute(Kaa_);
    implicitReady_ = (ldltA_.info() == Success);
    if (!implicitReady_) std::cerr << "EulerDG: active-block factorisation failed\n";
}

// Solve  K X = B  (K = M + gamma*dt*A) exactly via the partition: the decoupled
// dofs are a block-mass inverse (parallel); the small active sub-block is solved
// with its cached Cholesky factor (4 components on 4 threads).
MatrixXd EulerDG::implicitSolve(const MatrixXd& B) const {
    MatrixXd X = B; applyMassInverse(X);          // exact for all decoupled dofs (K=M)
    int na = (int)activeDofs_.size();
    if (na == 0) return X;
    MatrixXd Baa(na, 4);
    for (int a = 0; a < na; ++a) Baa.row(a) = B.row(activeDofs_[a]);
    MatrixXd Xaa(na, 4);
    std::array<std::thread, 4> th;
    for (int c = 0; c < 4; ++c)
        th[c] = std::thread([this, &Xaa, &Baa, c] { Xaa.col(c) = ldltA_.solve(Baa.col(c)); });
    for (auto& t : th) t.join();
    for (int a = 0; a < na; ++a) X.row(activeDofs_[a]) = Xaa.row(a);   // overwrite active dofs

    // (debug, EULERCHK) verify the partitioned solve == a full factorized solve of K=M+a_dt A
    static const bool CHK = std::getenv("EULERCHK") != nullptr;
    if (CHK) {
        SparseMatrix<double> Kfull = Mblk_ + aDtCached_ * A_;
        SimplicialLDLT<SparseMatrix<double>> full; full.compute(Kfull);
        MatrixXd Xf = full.solve(B);
        double d = (X - Xf).cwiseAbs().maxCoeff(), n = Xf.cwiseAbs().maxCoeff();
        std::cerr << "[chk] partitioned vs full solve  rel-diff=" << d / std::max(n, 1e-30)
                  << "  (active " << na << "/" << nDof_ << " dofs)\n";
    }
    return X;
}

// ===========================================================================
// Zhang-Shu (2010) positivity-preserving limiter: squeeze the polynomial toward
// its cell average until rho >= eps_rho and p >= eps_p at every volume + edge
// quadrature point.  Preserves the cell average exactly (hence conservation)
// and is inactive in smooth cells (no accuracy loss).
// ===========================================================================
// One-cell Zhang-Shu squeeze (shared by the global limiter and the LTS list limiter).
// Core Zhang-Shu squeeze on ONE cell's locDof x 4 dof block (reference-element data
// only -- no mesh/area, so it is shared verbatim by the global limiter and the
// compact LTS limiter).  In place.
void EulerDG::limitCellCore(Matrix<double, Dynamic, 4>& Ue) const {
    const Impl& P = *P_;
    int nqv = static_cast<int>(P.wv.size());
    Vector4d Ubar = Vector4d::Zero();                       // cell average
    for (int q = 0; q < nqv; ++q) {
        Vector4d Uq = Vector4d::Zero();
        for (int i = 0; i < locDof_; ++i) Uq += P.phiV[q](i) * Ue.row(i).transpose();
        Ubar += P.wv(q) * Uq;
    }
    double rb = Ubar(0);
    double pbar = pressure(Ubar);
    if (rb <= 0.0 || pbar <= 0.0) return;                   // cell mean already bad -> CFL too large
    double eps_rho = std::min(1e-13, 1e-3 * rb);
    double eps_p   = std::min(1e-13, 1e-3 * pbar);

    // collect U_h at all volume + edge quadrature points
    auto evalAll = [&](std::vector<Vector4d>& pts) {
        pts.clear();
        for (int q = 0; q < nqv; ++q) {
            Vector4d Uq = Vector4d::Zero();
            for (int i = 0; i < locDof_; ++i) Uq += P.phiV[q](i) * Ue.row(i).transpose();
            pts.push_back(Uq);
        }
        // Both edge orientations and all 3 edges: this superset covers exactly
        // the (et,dir) trace points the flux samples for THIS element on each of
        // its faces, so every flux-quadrature state is guaranteed admissible.
        // (For a symmetric 1-D Gauss rule dir=0 and dir=1 are the same point set
        // reversed, so including both is defensive -- robust to an asymmetric
        // rule -- at ~2x edge-eval cost.)
        for (int et = 0; et < 3; ++et)
            for (int dir = 0; dir < 2; ++dir)
                for (int q = 0; q < P.nqe; ++q) {
                    const RowVectorXd& ph = P.ephi[et][dir][q];
                    Vector4d Uq = Vector4d::Zero();
                    for (int i = 0; i < locDof_; ++i) Uq += ph(i) * Ue.row(i).transpose();
                    pts.push_back(Uq);
                }
    };
    std::vector<Vector4d> pts; evalAll(pts);

    // (i) density
    double rho_min = rb;
    for (const auto& Up : pts) rho_min = std::min(rho_min, Up(0));
    double th_rho = 1.0;
    if (rho_min < eps_rho) th_rho = (rb - eps_rho) / (rb - rho_min);
    th_rho = std::min(1.0, std::max(0.0, th_rho));

    // (ii) pressure (on the density-limited segment)
    double th_p = 1.0;
    for (const auto& Up : pts) {
        Vector4d U1 = Ubar + th_rho * (Up - Ubar);          // density-limited state
        if (pressure(U1) >= eps_p) continue;
        Vector4d d = U1 - Ubar;
        double A = (GAMMA - 1.0) * (d(3) * d(0) - 0.5 * (d(1) * d(1) + d(2) * d(2)));
        double B = (GAMMA - 1.0) * (Ubar(3) * d(0) + d(3) * rb - (Ubar(1) * d(1) + Ubar(2) * d(2))) - eps_p * d(0);
        double C = (GAMMA - 1.0) * (Ubar(3) * rb - 0.5 * (Ubar(1) * Ubar(1) + Ubar(2) * Ubar(2))) - eps_p * rb;
        double troot = 1.0;
        if (std::abs(A) < 1e-300) { if (std::abs(B) > 1e-300) troot = -C / B; }
        else {
            double disc = B * B - 4.0 * A * C;
            if (disc < 0) disc = 0;
            double sq = std::sqrt(disc);
            double r1 = (-B + sq) / (2.0 * A), r2 = (-B - sq) / (2.0 * A);
            troot = 1.0;
            for (double r : {r1, r2}) if (r >= -1e-12 && r <= 1.0 + 1e-12) troot = std::min(troot, std::max(0.0, r));
        }
        th_p = std::min(th_p, troot);
    }
    // theta_p is measured along the DENSITY-limited segment Ubar -> U^(1) =
    // Ubar + theta_rho(U_h - Ubar), so the two limiters COMPOSE: the final
    // squeeze in U_h coordinates is theta_rho * theta_p (NOT min).
    double theta = th_rho * th_p;
    if (theta >= 1.0) return;
    for (int i = 0; i < locDof_; ++i) {
        Vector4d c = Ue.row(i).transpose();
        Ue.row(i) = (Ubar + theta * (c - Ubar)).transpose();
    }
}

void EulerDG::limitCell(MatrixXd& U, int t) const {
    Matrix<double, Dynamic, 4> Ue(locDof_, 4);
    for (int i = 0; i < locDof_; ++i) Ue.row(i) = U.row(elem2dof_(t, i));
    limitCellCore(Ue);
    for (int i = 0; i < locDof_; ++i) U.row(elem2dof_(t, i)) = Ue.row(i);
}

void EulerDG::positivityLimit(MatrixXd& U) const {
    // each element is independent (reads/writes only its own DG block) -> parallel
    parallel_ranges(NT_, [&](int lo, int hi) { for (int t = lo; t < hi; ++t) limitCell(U, t); });
}

// LTS: limit only the given cells (the rest are untouched, e.g. the other level).
void EulerDG::positivityLimitCells(MatrixXd& U, const std::vector<int>& cells) const {
    int nc = static_cast<int>(cells.size());
    parallel_ranges(nc, [&](int lo, int hi) { for (int ci = lo; ci < hi; ++ci) limitCell(U, cells[ci]); });
}

void EulerDG::tvbLimit(MatrixXd&) const { /* optional minmod limiter -- not used (AV is primary) */ }

// ===========================================================================
// ARS(2,2,2) IMEX-RK step.  Explicit: inviscid flux load F(U)=inviscidResidual.
// Implicit: -A(eps) U (artificial viscosity), via K = M + gamma*dt*A.
//   gamma = 1 - sqrt2/2,  delta = -1/sqrt2.  Both tableaux stiffly accurate, so
//   U^{n+1} = U3 (the last implicit stage); only 2 flux evals + 2 implicit solves.
// ===========================================================================
bool EulerDG::step(double dt, double tEnd, const ExteriorStateFn& bc) {
    const double tN = tEnd - dt;
    const double g   = 1.0 - std::sqrt(2.0) / 2.0;     // gamma
    const double del = -1.0 / std::sqrt(2.0);          // delta
    const bool avOn = cfg_.use_av;

    // ---- optional phase profiler (EULERPROF=1) ----
    static const bool PROF = std::getenv("EULERPROF") != nullptr;
    static double tRefresh = 0, tResid = 0, tSolve = 0, tLimit = 0, tTotal = 0;
    auto clk = [] { return std::chrono::high_resolution_clock::now(); };
    auto secs = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    auto tStepStart = clk();

    if (avOn) {
        bool dtChanged = std::abs(g * dt - aDtCached_) > 1e-14 * std::max(1.0, std::abs(aDtCached_));
        if (!implicitReady_ || refreshCountdown_ <= 0 || dtChanged) {
            auto a = clk();
            computeViscosity(U_, epsFrozen_);
            refreshImplicit(epsFrozen_, g * dt);
            if (PROF) tRefresh += secs(a, clk());
            refreshCountdown_ = std::max(1, cfg_.av_refresh);
        }
        --refreshCountdown_;
    }

    auto Fload = [&](const MatrixXd& U, double t) { MatrixXd R; inviscidResidual(U, t, bc, R); return R; };
    auto solveK = [&](const MatrixXd& B) -> MatrixXd {
        if (!avOn) { MatrixXd X = B; applyMassInverse(X); return X; }   // K = M
        return implicitSolve(B);
    };

    MatrixXd MUn; applyMass(U_, MUn);                   // M U^n (block-diagonal, parallel)

    // Stage 1 (explicit node c=0)
    auto a1 = clk(); MatrixXd f1 = Fload(U_, tN); if (PROF) tResid += secs(a1, clk());

    // Stage 2:  (M + g dt A) U2 = M U^n + g dt f1
    auto a2 = clk(); MatrixXd U2 = solveK(MUn + (g * dt) * f1); if (PROF) tSolve += secs(a2, clk());
    auto a3 = clk(); if (cfg_.use_positivity) positivityLimit(U2); if (PROF) tLimit += secs(a3, clk());
    auto a4 = clk(); MatrixXd f2 = Fload(U2, tN + g * dt); if (PROF) tResid += secs(a4, clk());

    // Stage 3:  (M + g dt A) U3 = M U^n + dt(delta f1 + (1-delta) f2) + (1-g) dt (-A U2)
    MatrixXd rhs3 = MUn + dt * (del * f1 + (1.0 - del) * f2);
    if (avOn && activeDofs_.size())
        rhs3.noalias() -= ((1.0 - g) * dt) * (A_ * U2);          // explicit -A U2 term (A_ sparse)
    auto a5 = clk(); MatrixXd U3 = solveK(rhs3); if (PROF) tSolve += secs(a5, clk());
    auto a6 = clk(); if (cfg_.use_positivity) positivityLimit(U3); if (PROF) tLimit += secs(a6, clk());

    U_ = U3;
    ++nstep_;
    if (PROF) {
        tTotal += secs(tStepStart, clk());
        if (nstep_ % 20 == 0) {
            double phases = tRefresh + tResid + tSolve + tLimit;
            std::cerr << "[prof] step " << nstep_ << " ms/step: total=" << (1000*tTotal/nstep_)
                      << " | resid=" << (1000*tResid/nstep_) << " solve=" << (1000*tSolve/nstep_)
                      << " limit=" << (1000*tLimit/nstep_) << " refresh=" << (1000*tRefresh/nstep_)
                      << " arith/alloc=" << (1000*(tTotal-phases)/nstep_) << "\n";
        }
    }
    return true;
}

// ===========================================================================
// LOCAL TIME STEPPING (2-level, flux-register conservative subcycling)
//
// Conserved-variable DG advances cell K by  U_K += M_K^{-1} (vol - oint Hhat phi).
// The ONLY inter-cell coupling is the face flux Hhat; conservation between two
// cells sharing a face requires the time-integral of Hhat through that face to be
// equal-and-opposite for them.  Under LTS a coarse cell (one big step dt_c) and a
// fine cell (two substeps dt_c/2) would integrate that flux differently, breaking
// conservation -- the FLUX REGISTER fixes it: the coarse cell's provisional
// interface-flux contribution Phi_c (computed with a frozen fine ghost) is
// replaced by the time-integrated fine-resolved one Phi_f, by adding
// M^{-1}(Phi_f - Phi_c) to the coarse cell.  Because HLLC is conservative
// (Hhat(a,b,n) = -Hhat(b,a,-n)) both sides then see the SAME time-integrated flux
// -> exact conservation, independent of the (1st/2nd order) ghost prediction.
// ===========================================================================

// Masked inviscid residual: only cells of the advancing level get a residual; a
// face to a same-level neighbour reads the live U, a face to the OTHER level reads
// the frozen snapshot Ufrozen (the cross-level ghost).
void EulerDG::inviscidResidualMaskedCompact(const MatrixXd& Uc, const std::vector<int>& cells,
                                            double t, const ExteriorStateFn& bc, MatrixXd& Rc) const {
    const Impl& P = *P_;
    const int nc = static_cast<int>(cells.size());
    Rc.resize((size_t)nc * locDof_, 4);                     // every compact row is written below
    const int nqv = static_cast<int>(P.wv.size());
    const bool hllcFlux = cfg_.use_hllc;
    parallel_ranges(nc, [&](int lo, int hi) {
        Vector4d Fx, Fy;
        MatrixXd Ue(locDof_, 4), Unb(locDof_, 4), G(2, locDof_);
        Matrix<double, Dynamic, 4> Re(locDof_, 4);
        for (int ci = lo; ci < hi; ++ci) {
            int tt = cells[ci];
            int cls = cellClass_[tt];                           // this cell's dt-class
            Ue = Uc.middleRows((size_t)ci * locDof_, locDof_);
            Re.setZero();
            double area = fem_.area(tt);
            for (int q = 0; q < nqv; ++q) {
                Vector4d Uq = (P.phiV[q] * Ue).transpose();
                fluxes(Uq, Fx, Fy);
                G.noalias() = fem_.Dlam[tt] * P.dphiV[q];
                double wa = P.wv(q) * area;
                Re.noalias() += (wa * G.row(0).transpose()) * Fx.transpose();
                Re.noalias() += (wa * G.row(1).transpose()) * Fy.transpose();
            }
            for (int k = 0; k < 3; ++k) {
                const Impl::EE& r = P.ee[tt][k];
                bool interior = (r.nb != -1);
                if (interior) {
                    if (cellClass_[r.nb] == cls)                // same class -> compact working array
                        Unb = Uc.middleRows((size_t)g2c_[r.nb] * locDof_, locDof_);
                    else                                        // cross class -> macro-start snapshot (global)
                        for (int i = 0; i < locDof_; ++i) Unb.row(i) = ltsSnap_.row(elem2dof_(r.nb, i));
                }
                int tag = edgeTag_.size() ? edgeTag_(r.ei) : 0;
                for (int q = 0; q < P.nqe; ++q) {
                    const RowVectorXd& pa = P.ephi[r.et][r.dir][q];
                    Vector4d Um = (pa * Ue).transpose();
                    Vector4d Up;
                    if (interior) {
                        Up = (P.ephi[r.et_nb][r.dir_nb][q] * Unb).transpose();
                    } else {
                        double l1 = P.quad1d(q, 0), l2 = P.quad1d(q, 1);
                        Vector2d Pp = l1 * mesh_.node.row(r.n1).transpose() + l2 * mesh_.node.row(r.n2).transpose();
                        Up = bc ? bc(Pp.x(), Pp.y(), t, Um, r.nx, r.ny, tag) : Um;
                    }
                    Vector4d Hn = hllcFlux ? hllc(Um, Up, r.nx, r.ny) : rusanov(Um, Up, r.nx, r.ny);
                    double whe = P.w1d(q) * r.he;
                    Re.noalias() -= (whe * pa.transpose()) * Hn.transpose();
                }
            }
            Rc.middleRows((size_t)ci * locDof_, locDof_) = Re;
        }
    });
}

// Flux register: accumulate  weight * oint Hhat * phi_coarse  (the COARSE cell's
// interface face-term) into reg, over the interface edges, for the advancing level.
// Serial over the (few) interface edges -> race-free (the parallel element loop
// only writes own blocks).  The advancing cell supplies the flux (its outward
// normal); the coarse cell supplies the test functions and the reg target.  In the
// coarse phase the normal is the coarse outward normal, in the fine phase the fine
// outward normal (= -coarse) -- so the two phases' contributions cancel for a
// steady/consistent state, leaving exactly the time-resolution difference Phi_f-Phi_c.
void EulerDG::accumulateReflux(const MatrixXd& Uc, int advClass, double weight,
                               std::vector<MatrixXd>& reg) const {
    const Impl& P = *P_;
    const bool hllcFlux = cfg_.use_hllc;
    MatrixXd Uadv(locDof_, 4), Uoth(locDof_, 4);
    for (const LtsEdge& le : ltsEdges_) {
        int cc = cellClass_[le.coarse];                         // the coarse (higher) class
        bool advIsCoarse;
        if (advClass == cc)          advIsCoarse = true;        // advancing the coarse side -> Phi_c
        else if (advClass == cc - 1) advIsCoarse = false;       // advancing the fine side   -> Phi_f
        else continue;                                          // this interface not touched by advClass
        const Impl::EE& eeC = P.ee[le.coarse][le.coarseK];      // coarse cell: test fns + reg target
        int advCell = advIsCoarse ? le.coarse : le.fine;        // advancing cell -> compact stage array
        const Impl::EE& eeAdv = advIsCoarse ? eeC : P.ee[le.fine][le.fineK];
        const Impl::EE& eeOth = advIsCoarse ? P.ee[le.fine][le.fineK] : eeC;
        int othCell = advIsCoarse ? le.fine : le.coarse;
        int adv0 = g2c_[advCell] * locDof_;
        for (int i = 0; i < locDof_; ++i) {
            Uadv.row(i) = Uc.row(adv0 + i);                     // advancing class -> compact stage array
            Uoth.row(i) = ltsSnap_.row(elem2dof_(othCell, i));  // cross-class ghost = macro-start snapshot
        }
        Matrix<double, Dynamic, 4> contrib = Matrix<double, Dynamic, 4>::Zero(locDof_, 4);
        for (int q = 0; q < P.nqe; ++q) {
            Vector4d Um = (P.ephi[eeAdv.et][eeAdv.dir][q] * Uadv).transpose();
            Vector4d Up = (P.ephi[eeOth.et][eeOth.dir][q] * Uoth).transpose();
            Vector4d Hn = hllcFlux ? hllc(Um, Up, eeAdv.nx, eeAdv.ny) : rusanov(Um, Up, eeAdv.nx, eeAdv.ny);
            const RowVectorXd& phiC = P.ephi[eeC.et][eeC.dir][q];
            double whe = P.w1d(q) * eeAdv.he;
            contrib.noalias() += (weight * whe * phiC.transpose()) * Hn.transpose();
        }
        for (int i = 0; i < locDof_; ++i) reg[cc].row(elem2dof_(le.coarse, i)) += contrib.row(i);
    }
}

// One masked ARS(2,2,2) substep of class `advClass`, on COMPACT per-class arrays
// (size nc*locDof, row ci*locDof+i <-> cell cells[ci], local dof i) -- so every op is
// O(class), never O(nDof).  Implicit AV only for class 0 (eps>0 lives on the finest
// cells); coarser classes are pure explicit.  Cross-class faces read the macro-start
// snapshot ltsSnap_.  Interface fluxes fold into the (global) registers.
bool EulerDG::advanceMaskedARS(const std::vector<int>& cells, int advClass, double dt, double tN,
                               const ExteriorStateFn& bc, std::vector<MatrixXd>& reg) {
    const double g   = 1.0 - std::sqrt(2.0) / 2.0;
    const double del = -1.0 / std::sqrt(2.0);
    const bool avOn = cfg_.use_av && (advClass == 0);   // AV only on the finest class
    const int nc = static_cast<int>(cells.size());
    if (nc == 0) return true;

    for (int ci = 0; ci < nc; ++ci) g2c_[cells[ci]] = ci;   // global cell -> compact index (this class)

    MatrixXd Uc((size_t)nc * locDof_, 4);                   // gather U^n of this class
    parallel_ranges(nc, [&](int lo, int hi) {
        for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
            for (int i = 0; i < locDof_; ++i) Uc.row((size_t)ci * locDof_ + i) = U_.row(elem2dof_(t, i)); }
    });

    // Compact rows of the active (AV) dofs.  The AV operator couples a class-0 shock
    // cell to its (possibly class-1) neighbour, so some active dofs are OFF this class
    // -> mark them -1 (their RHS is 0 and their result is discarded, exactly as the
    // full-nDof path did via zeroed off-class rows; they still couple inside A).
    std::vector<int> actRow;
    if (avOn && activeDofs_.size()) {
        int na = (int)activeDofs_.size(); actRow.resize(na);
        for (int a = 0; a < na; ++a) { int gd = activeDofs_[a]; int ci = g2c_[gd / locDof_];
            actRow[a] = (ci >= 0) ? ci * locDof_ + gd % locDof_ : -1; }
    }

    auto massC = [&](const MatrixXd& X, MatrixXd& MX) {
        MX.resize((size_t)nc * locDof_, 4);
        parallel_ranges(nc, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Ue(locDof_, 4);
            for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
                Ue.noalias() = (fem_.area(t) * P_->Mref) * X.middleRows((size_t)ci * locDof_, locDof_);
                MX.middleRows((size_t)ci * locDof_, locDof_) = Ue; }
        });
    };
    auto massInvC = [&](MatrixXd& R) {
        parallel_ranges(nc, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Re(locDof_, 4);
            for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
                Re.noalias() = (MrefInv_ / fem_.area(t)) * R.middleRows((size_t)ci * locDof_, locDof_);
                R.middleRows((size_t)ci * locDof_, locDof_) = Re; }
        });
    };
    auto limitC = [&](MatrixXd& X) {
        parallel_ranges(nc, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Ue(locDof_, 4);
            for (int ci = lo; ci < hi; ++ci) {
                Ue = X.middleRows((size_t)ci * locDof_, locDof_);
                limitCellCore(Ue);
                X.middleRows((size_t)ci * locDof_, locDof_) = Ue; }
        });
    };
    const int na = avOn ? (int)activeDofs_.size() : 0;      // active (AV) dofs; 0 for explicit classes
    auto activeSolve = [&](const MatrixXd& Baa) -> MatrixXd {   // K_aa Xaa = Baa via the cached Cholesky
        MatrixXd Xaa(na, 4);
        std::array<std::thread, 4> th;
        for (int c = 0; c < 4; ++c) th[c] = std::thread([&, c] { Xaa.col(c) = ldltA_.solve(Baa.col(c)); });
        for (auto& tt : th) tt.join();
        return Xaa;
    };
    auto Fload = [&](const MatrixXd& X, double tt) { MatrixXd R; inviscidResidualMaskedCompact(X, cells, tt, bc, R); return R; };
    // gather a row of the active RHS: compact value where the dof is on this class, 0 off-class
    auto gatherActive = [&](const MatrixXd& Bc, MatrixXd& Baa) {
        for (int a = 0; a < na; ++a) { if (actRow[a] >= 0) Baa.row(a) = Bc.row(actRow[a]); else Baa.row(a).setZero(); }
    };

    MatrixXd MUn; massC(Uc, MUn);                             // M U^n on this class only
    MatrixXd f1 = Fload(Uc, tN);
    accumulateReflux(Uc, advClass, del * dt, reg);           // stage-1 weight delta*dt

    // stage 1:  (M + g dt A) U2 = M U^n + g dt f1
    MatrixXd B1 = MUn + (g * dt) * f1;
    MatrixXd U2 = B1; massInvC(U2);                          // decoupled block-mass inverse (K = M)
    MatrixXd Xaa1;
    if (na) {
        MatrixXd Baa(na, 4); gatherActive(B1, Baa);
        Xaa1 = activeSolve(Baa);                            // full active solution (incl. off-class dofs)
        for (int a = 0; a < na; ++a) if (actRow[a] >= 0) U2.row(actRow[a]) = Xaa1.row(a);
    }
    if (cfg_.use_positivity) limitC(U2);

    MatrixXd f2 = Fload(U2, tN + g * dt);
    accumulateReflux(U2, advClass, (1.0 - del) * dt, reg);   // stage-2 weight (1-delta)*dt

    // stage 2:  (M + g dt A) U3 = M U^n + dt(del f1 + (1-del) f2) - (1-g) dt A U2
    MatrixXd rhs3 = MUn + dt * (del * f1 + (1.0 - del) * f2);
    MatrixXd U3 = rhs3; massInvC(U3);
    if (na) {
        MatrixXd AU = Aaa_ * Xaa1;                          // (A U2) at ALL active dofs (active indexing)
        MatrixXd Baa(na, 4);
        for (int a = 0; a < na; ++a) {
            Eigen::RowVector4d b; if (actRow[a] >= 0) b = rhs3.row(actRow[a]); else b.setZero();
            Baa.row(a) = b - ((1.0 - g) * dt) * AU.row(a);
        }
        MatrixXd Xaa3 = activeSolve(Baa);
        for (int a = 0; a < na; ++a) if (actRow[a] >= 0) U3.row(actRow[a]) = Xaa3.row(a);
    }
    if (cfg_.use_positivity) limitC(U3);

    parallel_ranges(nc, [&](int lo, int hi) {                 // commit (scatter compact -> U_)
        for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
            for (int i = 0; i < locDof_; ++i) U_.row(elem2dof_(t, i)) = U3.row((size_t)ci * locDof_ + i); }
    });
    for (int ci = 0; ci < nc; ++ci) g2c_[cells[ci]] = -1;   // restore (keep g2c_ all -1 between advances)
    return true;
}

// Recursive subcycling: advance the hierarchy for ONE step of class `level`'s clock
// (dt), advancing class `level` once and the finer hierarchy twice (dt/2), then
// refluxing class `level` at its interface with class level-1.
void EulerDG::ltsIntegrate(int level, double dt, double tN, const ExteriorStateFn& bc,
                           std::vector<MatrixXd>& reg, int levels) {
    if (level >= 1) reg[level].setZero(nDof_, 4);             // this level's register, freshly accumulated below
    advanceMaskedARS(classCells_[level], level, dt, tN, bc, reg);   // coarse-first (provisional Phi_c)
    if (level == 0) return;                                  // finest: no finer hierarchy, no reflux
    ltsIntegrate(level - 1, 0.5 * dt, tN,            bc, reg, levels);   // first sub-step
    ltsIntegrate(level - 1, 0.5 * dt, tN + 0.5 * dt, bc, reg, levels);   // second sub-step
    // reflux class-`level` cells:  U += M^{-1} (Phi_f - Phi_c)  (reg[level] is the residual)
    const std::vector<int>& cl = classCells_[level];
    parallel_ranges(static_cast<int>(cl.size()), [&](int lo, int hi) {
        Matrix<double, Dynamic, 4> Re(locDof_, 4);
        for (int ci = lo; ci < hi; ++ci) { int t = cl[ci];
            for (int i = 0; i < locDof_; ++i) Re.row(i) = reg[level].row(elem2dof_(t, i));
            Re = (MrefInv_ / fem_.area(t)) * Re;
            for (int i = 0; i < locDof_; ++i) U_.row(elem2dof_(t, i)) += Re.row(i); }
    });
    if (cfg_.use_positivity) positivityLimitCells(U_, cl);   // refluxed cells may need a clamp
}

bool EulerDG::stepLTS(double dtMacro, double tEnd, const ExteriorStateFn& bc,
                      const std::vector<int>& cellClass, int levels) {
    cellClass_ = cellClass;
    if ((int)g2c_.size() != NT_) g2c_.assign(NT_, -1);       // compact-index scratch (kept all -1 between advances)
    ltsSnap_ = U_;                                           // cross-class ghost: state at macro start
    const double tN  = tEnd - dtMacro;
    const double g   = 1.0 - std::sqrt(2.0) / 2.0;
    const double dt0 = dtMacro / static_cast<double>(1 << (levels - 1));   // finest-class step

    // compact per-class cell lists
    classCells_.assign(levels, {});
    for (int t = 0; t < NT_; ++t) classCells_[cellClass_[t]].push_back(t);

    // class-interface edge list ("coarse" = higher class = bigger dt), deduped by t<nb
    ltsEdges_.clear();
    for (int t = 0; t < NT_; ++t)
        for (int k = 0; k < 3; ++k) {
            int nb = P_->ee[t][k].nb;
            if (nb < 0 || t > nb) continue;
            if (cellClass_[t] == cellClass_[nb]) continue;
            int coarse = cellClass_[t] > cellClass_[nb] ? t : nb;
            int fine   = cellClass_[t] > cellClass_[nb] ? nb : t;
            int ck = -1, fk = -1, ei = P_->ee[t][k].ei;
            for (int kk = 0; kk < 3; ++kk) { if (P_->ee[coarse][kk].ei == ei) ck = kk; if (P_->ee[fine][kk].ei == ei) fk = kk; }
            ltsEdges_.push_back({coarse, fine, ck, fk});
        }

    // artificial viscosity lives on class 0 (the finest) and advances with dt0.  A
    // macro step spans 2^(levels-1) fine substeps, so the shock viscosity must be
    // refreshed EVERY macro (freezing it across av_refresh macros lets a Mach-10
    // shock outrun its viscosity band and blow up).  CONFINE the AV operator to
    // class 0 (zero eps elsewhere): the implicit solve only fills class-0 cells, so
    // A coupling to another class would read uninitialised RHS -> NaN.
    static const bool LP = std::getenv("LTSPROF") != nullptr;   // micro-profile the macro step
    static double tAV = 0, tInt = 0; static int np = 0;
    auto cA = std::chrono::high_resolution_clock::now();
    if (cfg_.use_av) {
        computeViscosity(U_, epsFrozen_);
        for (int t = 0; t < NT_; ++t) if (cellClass_[t] != 0) epsFrozen_(t) = 0.0;   // zero eps off class 0
        refreshImplicit(epsFrozen_, g * dt0);
    }
    if (LP) tAV += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - cA).count();

    std::vector<MatrixXd> reg(levels);                       // reg[c] = flux register of coarse class c (c>=1)
    auto c0 = std::chrono::high_resolution_clock::now();
    ltsIntegrate(levels - 1, dtMacro, tN, bc, reg, levels);
    if (LP) {
        auto c1 = std::chrono::high_resolution_clock::now();
        tInt += std::chrono::duration<double>(c1 - c0).count(); ++np;
        if (np % 50 == 0) std::cerr << "[ltsprof] AV-refresh=" << (1000*tAV/np) << " integrate="
                                    << (1000*tInt/np) << " ms/macro (" << np << " macros)\n";
    }
    ++nstep_;
    return true;
}

} // namespace euler
