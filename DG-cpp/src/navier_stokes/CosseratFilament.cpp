#include "CosseratFilament.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace ns {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::SparseMatrix;
using Eigen::Triplet;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

void CosseratFilament::initStraight(int N_, double xs, double ys, double xe, double ye,
                                    double rhoLineIn, double EAin, double EIin)
{
    N = N_;
    EA = EAin; EI = EIin; rhoLine = rhoLineIn;

    X0.setZero(N + 1, 2);
    for (int i = 0; i <= N; ++i) {
        double t = (double)i / (double)N;
        X0(i, 0) = xs + t * (xe - xs);
        X0(i, 1) = ys + t * (ye - ys);
    }
    X = X0;
    V.setZero(N + 1, 2);
    A.setZero(N + 1, 2);

    l0.resize(N);
    for (int i = 0; i < N; ++i) {
        Vector2d e = X0.row(i + 1) - X0.row(i);
        l0(i) = e.norm();
    }
    lbar.setZero(N + 1);
    lbar(0) = 0.5 * l0(0);
    for (int i = 1; i < N; ++i) lbar(i) = 0.5 * (l0(i - 1) + l0(i));
    lbar(N) = 0.5 * l0(N - 1);

    m = rhoLine * lbar;

    Fext.setZero(2 * (N + 1));
    frozen.clear(); frozenVal.resize(0);
    frozenVel.resize(0); frozenAcc.resize(0);
}

void CosseratFilament::clampRoot(bool alsoClampNode1)
{
    frozen.clear();
    std::vector<double> vals;
    auto add = [&](int node, int dim) {
        frozen.push_back(2 * node + dim);
        vals.push_back(X0(node, dim));
    };
    add(0, 0); add(0, 1);
    if (alsoClampNode1) { add(1, 0); add(1, 1); }
    frozenVal.resize((int)vals.size());
    for (int k = 0; k < (int)vals.size(); ++k) frozenVal(k) = vals[k];
    frozenVel = VectorXd::Zero((int)vals.size());
    frozenAcc = VectorXd::Zero((int)vals.size());
}

// ---------------------------------------------------------------------------
// Internal-energy gradients
// ---------------------------------------------------------------------------

namespace {

// Stretch contribution from edge i (between nodes i, i+1):
//   V_s = (1/(2 l0)) * EA * (l - l0)^2,  l = ||X_{i+1}-X_i||
// dV/dX_i      = -EA*(l-l0)/l0 * (e/l)
// dV/dX_{i+1}  = +EA*(l-l0)/l0 * (e/l)
// Hessian (per edge):
//   d^2V / dX_a dX_b  with a,b in {i, i+1}, sign +1 on the diagonal blocks
//   and -1 off-diagonal.  The 2x2 block is
//       EA/l0 * [ (1) (e e^T)/l^2    +    (l - l0) * (I - (e e^T)/l^2) / l ]
inline void stretchGradHess(const Vector2d& Xi, const Vector2d& Xj, double l0,
                            double EA,
                            Vector2d& gi, Vector2d& gj,
                            Eigen::Matrix2d& Hii)
{
    Vector2d e = Xj - Xi;
    double l = e.norm();
    if (l < 1e-300) { gi.setZero(); gj.setZero(); Hii.setZero(); return; }
    Vector2d t = e / l;
    double s = (l - l0);
    double f = EA * s / l0;
    gi = -f * t;
    gj = +f * t;
    Eigen::Matrix2d ttT = t * t.transpose();
    Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
    // d/dXj of f*t   ->   EA/l0 * (t t^T)  +  EA*s/l0 * (I - t t^T)/l
    Hii = (EA / l0) * ttT + (EA * s / l0) * (I - ttT) / l;
    // by symmetry, the four 2x2 sub-blocks of the per-edge Hessian are
    //     H(Xi,Xi)= +Hii    H(Xj,Xj)= +Hii
    //     H(Xi,Xj)= -Hii    H(Xj,Xi)= -Hii
}

// Discrete turning angle at interior node i from edges e_a = X_i - X_{i-1},
// e_b = X_{i+1} - X_i:   theta = atan2( cross(e_a, e_b), dot(e_a, e_b) ).
// Bending energy at node i:  V_b = (EI / (2 lbar_i)) * theta^2.
//
// Closed-form gradient (Bergou et al. 2008, eq. 7):
//   dtheta/dX_{i-1} = -kb / |e_a|
//   dtheta/dX_{i+1} = +kb / |e_b|
//   dtheta/dX_i     = -(dtheta/dX_{i-1} + dtheta/dX_{i+1})
// where kb = (rotated by +pi/2 in 2D of the bisector) — i.e. in 2D this
// reduces to the unit normal to the edge multiplied by the appropriate sign.
// We compute  dtheta/dX_*  in 2D directly:
//
//   d theta / d X_{i-1} = -(1/|e_a|) * n_a
//   d theta / d X_{i+1} = +(1/|e_b|) * n_b
//   d theta / d X_i     = +(1/|e_a|) n_a - (1/|e_b|) n_b
// where n_a = perp(e_a)/|e_a|, n_b = perp(e_b)/|e_b|, perp((x,y)) = (-y, x).
inline void bendGradAtInteriorNode(const Vector2d& Xim1, const Vector2d& Xi, const Vector2d& Xip1,
                                   double EIc, double lbar,
                                   Vector2d& dim1, Vector2d& di, Vector2d& dip1,
                                   double& thetaOut)
{
    Vector2d ea = Xi   - Xim1;
    Vector2d eb = Xip1 - Xi;
    double la = ea.norm(), lb = eb.norm();
    if (la < 1e-300 || lb < 1e-300) {
        dim1.setZero(); di.setZero(); dip1.setZero();
        thetaOut = 0.0;
        return;
    }
    double cr = ea.x() * eb.y() - ea.y() * eb.x();
    double dt = ea.x() * eb.x() + ea.y() * eb.y();
    double theta = std::atan2(cr, dt);
    thetaOut = theta;
    // 2D unit normals (perp / norm)
    Vector2d na(-ea.y() / la, ea.x() / la);
    Vector2d nb(-eb.y() / lb, eb.x() / lb);
    // d theta / d X_*  in 2D.  Derivation: d(theta)/d e_a = -perp(e_a)/|e_a|^2,
    // d(theta)/d e_b = +perp(e_b)/|e_b|^2, then chain through e_a = X_i - X_{i-1}
    // (so dim1 picks up an extra minus sign) and e_b = X_{i+1} - X_i.
    Vector2d dthetadim1 = +(1.0 / la) * na;     // = +perp(e_a)/|e_a|^2
    Vector2d dthetadip1 = +(1.0 / lb) * nb;     // = +perp(e_b)/|e_b|^2
    Vector2d dthetadi   = -(dthetadim1 + dthetadip1);
    // V = (EI/(2 lbar)) theta^2  ->  dV/dX = (EI/lbar) theta * dtheta/dX
    double f = (EIc / lbar) * theta;
    dim1 = f * dthetadim1;
    di   = f * dthetadi;
    dip1 = f * dthetadip1;
}

} // anon

// ---------------------------------------------------------------------------
// Residual & tangent
// ---------------------------------------------------------------------------

void CosseratFilament::assembleResidualAndTangent(const MatrixXd& Xt, double dt,
                                                  const MatrixXd& Xn, const MatrixXd& Vn, const MatrixXd& An,
                                                  VectorXd& res, SparseMatrix<double>& K) const
{
    const int Nv = N + 1;
    const int dim = 2 * Nv;
    res.setZero(dim);

    // Newmark-beta with beta=1/4, gammaN=1/2.
    const double beta = 0.25;
    const double gammaN = 0.5;
    const double idt2 = 1.0 / (beta * dt * dt);
    const double idt  = gammaN / (beta * dt);

    // Predicted velocity / acceleration as functions of trial displacement Xt:
    //   A_p = (Xt - Xn)/(beta dt^2) - (1/(beta dt)) Vn - (1/(2 beta) - 1) An
    //   V_p = Vn + dt * ((1 - gammaN) An + gammaN A_p)
    // For the dynamic equilibrium we only need A_p (and V_p for damping).
    Eigen::MatrixXd Ap = (Xt - Xn) * idt2
                       - (1.0 / (beta * dt)) * Vn
                       - (1.0 / (2.0 * beta) - 1.0) * An;
    Eigen::MatrixXd Vp = Vn + dt * ((1.0 - gammaN) * An + gammaN * Ap);

    // ---- inertial + damping residual + diagonal stiffness M / (beta dt^2) ----
    // Residual_inertial = m * A_p + dampStr * m * V_p + drag * (V_p - V_ref)
    bool hasDrag = (dragCoef.size() == Nv) && (dragRef.rows() == Nv);
    std::vector<Triplet<double>> trips;
    trips.reserve(64 * Nv);
    for (int i = 0; i < Nv; ++i) {
        double mi = m(i);
        double ci = hasDrag ? dragCoef(i) : 0.0;
        double diag = mi * idt2 + (dampStr * mi + ci) * idt;
        for (int d = 0; d < 2; ++d) {
            res(2 * i + d) += mi * (Ap(i, d) + dampStr * Vp(i, d));
            if (hasDrag) res(2 * i + d) += ci * (Vp(i, d) - dragRef(i, d));
            trips.emplace_back(2 * i + d, 2 * i + d, diag);
        }
    }

    // ---- stretch contribution (per edge i = 0..N-1) ----
    for (int i = 0; i < N; ++i) {
        Vector2d Xi = Xt.row(i), Xj = Xt.row(i + 1);
        Vector2d gi, gj; Eigen::Matrix2d Hii;
        stretchGradHess(Xi, Xj, l0(i), EA, gi, gj, Hii);
        for (int d = 0; d < 2; ++d) {
            res(2 * i + d)       += gi(d);
            res(2 * (i + 1) + d) += gj(d);
        }
        // 4x4 stiffness on (i, i+1): [+Hii, -Hii; -Hii, +Hii]
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) {
                double v = Hii(a, b);
                trips.emplace_back(2 * i + a,       2 * i + b,       +v);
                trips.emplace_back(2 * (i+1) + a,   2 * (i+1) + b,   +v);
                trips.emplace_back(2 * i + a,       2 * (i+1) + b,   -v);
                trips.emplace_back(2 * (i+1) + a,   2 * i + b,       -v);
            }
    }

    // ---- bending contribution: gradient analytic, Hessian by finite difference ----
    // We need K = d residual / d Xt to converge Newton; for bending we get
    // d (grad V_b) / d X_* by central-differencing the analytic gradient on
    // the 6-DOF stencil per interior node.  This is cheap (O(N) 6x6 systems)
    // and dominates only the assembly cost, which is small compared with the
    // sparse linear solve.
    auto bendGradAt = [&](const MatrixXd& Xx, int i,
                          Vector2d& dim1, Vector2d& di, Vector2d& dip1) {
        double th;
        bendGradAtInteriorNode(Xx.row(i - 1), Xx.row(i), Xx.row(i + 1),
                               EI, lbar(i), dim1, di, dip1, th);
    };
    // Adaptive FD step:  scale with the typical edge length so we resolve the
    // bending gradient even when the rod is short or under heavy load.  Rule
    // of thumb  fdEps ~ 1e-4 * lbar  works robustly across stiffness regimes
    // (much smaller hits cancellation; much larger truncates the Hessian).
    double lbarMax = 0.0;
    for (int i = 0; i < (int)lbar.size(); ++i) lbarMax = std::max(lbarMax, lbar(i));
    const double fdEps = std::max(1e-8, 1e-4 * std::max(1e-6, lbarMax));
    for (int i = 1; i < N; ++i) {
        Vector2d gA, gB, gC;
        bendGradAt(Xt, i, gA, gB, gC);
        for (int d = 0; d < 2; ++d) {
            res(2 * (i - 1) + d) += gA(d);
            res(2 * i       + d) += gB(d);
            res(2 * (i + 1) + d) += gC(d);
        }
        // 6x6 numerical Hessian via central differences.
        int idx[3] = { i - 1, i, i + 1 };
        for (int p = 0; p < 3; ++p) {
            int node_p = idx[p];
            for (int dp = 0; dp < 2; ++dp) {
                MatrixXd Xpert = Xt;
                Xpert(node_p, dp) += fdEps;
                Vector2d gAp, gBp, gCp; bendGradAt(Xpert, i, gAp, gBp, gCp);
                Xpert(node_p, dp) -= 2 * fdEps;
                Vector2d gAm, gBm, gCm; bendGradAt(Xpert, i, gAm, gBm, gCm);
                // column = d g / d (node_p, dp)
                Vector2d col[3] = {
                    (gAp - gAm) / (2 * fdEps),
                    (gBp - gBm) / (2 * fdEps),
                    (gCp - gCm) / (2 * fdEps),
                };
                for (int q = 0; q < 3; ++q) {
                    int node_q = idx[q];
                    for (int dq = 0; dq < 2; ++dq) {
                        trips.emplace_back(2 * node_q + dq, 2 * node_p + dp, col[q](dq));
                    }
                }
            }
        }
    }

    // ---- external load ----
    res -= Fext;

    // Build sparse K from triplets.
    K.resize(dim, dim);
    K.setFromTriplets(trips.begin(), trips.end());
}

void CosseratFilament::applyClamps(SparseMatrix<double>& K, VectorXd& res) const
{
    if (frozen.empty()) return;
    // Zero rows + cols of the frozen DOFs and put 1 on the diagonal; set the
    // RHS to 0 (Newton update for those DOFs is forced to zero).
    std::vector<char> mask(K.rows(), 0);
    for (int idx : frozen) mask[idx] = 1;
    for (int j = 0; j < K.outerSize(); ++j) {
        for (SparseMatrix<double>::InnerIterator it(K, j); it; ++it) {
            int r = (int)it.row(), c = (int)it.col();
            if (mask[r] || mask[c]) {
                if (r == c) it.valueRef() = 1.0;
                else        it.valueRef() = 0.0;
            }
        }
    }
    for (int idx : frozen) res(idx) = 0.0;
}

// ---------------------------------------------------------------------------
// Newmark step (Newton-Raphson)
// ---------------------------------------------------------------------------

bool CosseratFilament::step(double dt)
{
    const int Nv = N + 1;
    const int dim = 2 * Nv;

    MatrixXd Xn = X, Vn = V, An = A;
    MatrixXd Xt = X;          // initial guess: previous step

    // Hold the clamped DOFs at their target positions throughout the iteration.
    auto enforceClampPos = [&](MatrixXd& Y) {
        for (int k = 0; k < (int)frozen.size(); ++k) {
            int idx = frozen[k];
            int node = idx / 2, d = idx % 2;
            Y(node, d) = frozenVal(k);
        }
    };
    enforceClampPos(Xt);

    const int maxIt = 12;
    const double tol = 1e-10;
    Eigen::SimplicialLDLT<SparseMatrix<double>> solver;

    for (int it = 0; it < maxIt; ++it) {
        VectorXd res; SparseMatrix<double> K;
        assembleResidualAndTangent(Xt, dt, Xn, Vn, An, res, K);
        applyClamps(K, res);

        double rnrm = res.norm();
        if (rnrm < tol) break;

        solver.compute(K);
        if (solver.info() != Eigen::Success) {
            std::cerr << "CosseratFilament::step: factorisation failed\n";
            return false;
        }
        VectorXd dx = solver.solve(res);   // Newton: dx = K^{-1} res; X += -dx
        if (solver.info() != Eigen::Success) {
            std::cerr << "CosseratFilament::step: solve failed\n";
            return false;
        }
        for (int i = 0; i < Nv; ++i) {
            Xt(i, 0) -= dx(2 * i);
            Xt(i, 1) -= dx(2 * i + 1);
        }
        enforceClampPos(Xt);
        if (dx.norm() < tol) break;
    }

    // Update state with the converged Newmark relations.
    const double beta = 0.25;
    const double gammaN = 0.5;
    MatrixXd An1 = (Xt - Xn) / (beta * dt * dt)
                 - Vn / (beta * dt)
                 - (1.0 / (2.0 * beta) - 1.0) * An;
    MatrixXd Vn1 = Vn + dt * ((1.0 - gammaN) * An + gammaN * An1);

    // Force clamped velocities/accelerations to their prescribed support values.
    // Stationary supports use the zero defaults; moving clamps can fill
    // frozenVel/frozenAcc before step().
    for (int kk = 0; kk < (int)frozen.size(); ++kk) {
        int idx = frozen[kk];
        int node = idx / 2, d = idx % 2;
        Vn1(node, d) = (frozenVel.size() == (int)frozen.size()) ? frozenVel(kk) : 0.0;
        An1(node, d) = (frozenAcc.size() == (int)frozen.size()) ? frozenAcc(kk) : 0.0;
    }
    X = Xt; V = Vn1; A = An1;
    return true;
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

double CosseratFilament::tipDeflection() const
{
    Vector2d d = X.row(N).transpose() - X0.row(N).transpose();
    return d.norm();
}

double CosseratFilament::potentialEnergy() const
{
    double Vs = 0.0, Vb = 0.0;
    for (int i = 0; i < N; ++i) {
        Vector2d e = X.row(i + 1) - X.row(i);
        double l = e.norm();
        Vs += 0.5 * EA * (l - l0(i)) * (l - l0(i)) / l0(i);
    }
    for (int i = 1; i < N; ++i) {
        Vector2d ea = X.row(i)     - X.row(i - 1);
        Vector2d eb = X.row(i + 1) - X.row(i);
        double cr = ea.x() * eb.y() - ea.y() * eb.x();
        double dt = ea.x() * eb.x() + ea.y() * eb.y();
        double th = std::atan2(cr, dt);
        Vb += 0.5 * EI * th * th / lbar(i);
    }
    return Vs + Vb;
}

double CosseratFilament::kineticEnergy() const
{
    double T = 0.0;
    for (int i = 0; i <= N; ++i) {
        T += 0.5 * m(i) * (V(i, 0) * V(i, 0) + V(i, 1) * V(i, 1));
    }
    return T;
}

double CosseratFilament::maxEdgeStrain() const
{
    double s = 0.0;
    for (int i = 0; i < N; ++i) {
        Vector2d e = X.row(i + 1) - X.row(i);
        double l = e.norm();
        s = std::max(s, std::abs(l - l0(i)) / l0(i));
    }
    return s;
}

} // namespace ns
