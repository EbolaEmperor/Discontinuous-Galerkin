#include "ElasticTadpole.h"

#include "FEM.h"
#include "Mesh.h"
#include "NavierStokes.h"      // MeshLocator
#include "IBCoupler.h"         // meshToRod, rodArclengthWeights

#include <algorithm>
#include <cmath>
#include <iostream>

namespace ns {

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::RowVectorXd;

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------

void ElasticTadpole::init()
{
    // Head mass / inertia (disk).
    mass = rhoHead * M_PI * rHead * rHead;
    Inertia = 0.5 * mass * rHead * rHead;

    // Build a fresh elastic tail.  The tail starts at the rear surface of
    // the head (one rHead behind the centre along the body axis) and goes
    // another Ltail in the same direction (opposite the heading).
    Vector2d tan(-std::cos(theta), -std::sin(theta));     // body-axis-backward
    Vector2d X0 = Vector2d(x, y) + rHead * tan;
    Vector2d X1 = X0 + Ltail * tan;
    tail = std::make_unique<CosseratFilament>();
    tail->initStraight(Ntail, X0.x(), X0.y(), X1.x(), X1.y(),
                       rhoTail, KS, KB);
    tail->dampStr = tailDamp;
    // Clamp the first two nodes so the wall direction is locked.
    tail->clampRoot(true);
    syncTailClamp();
}

void ElasticTadpole::syncTailClamp()
{
    if (!tail || tail->N < 1) return;

    Vector2d tan(-std::cos(theta), -std::sin(theta));     // body-axis backward
    Vector2d clamp0 = Vector2d(x, y) + rHead * tan;
    double l0seg = Ltail / std::max(1, Ntail);
    Vector2d clamp1 = clamp0 + l0seg * tan;

    Vector2d r0 = clamp0 - Vector2d(x, y);
    Vector2d r1 = clamp1 - Vector2d(x, y);
    Vector2d v0(vx - omega * r0.y(), vy + omega * r0.x());
    Vector2d v1(vx - omega * r1.y(), vy + omega * r1.x());

    if ((int)tail->frozen.size() >= 4) {
        tail->frozenVal(0) = clamp0.x();
        tail->frozenVal(1) = clamp0.y();
        tail->frozenVal(2) = clamp1.x();
        tail->frozenVal(3) = clamp1.y();
        if (tail->frozenVel.size() != (int)tail->frozen.size())
            tail->frozenVel = VectorXd::Zero((int)tail->frozen.size());
        if (tail->frozenAcc.size() != (int)tail->frozen.size())
            tail->frozenAcc = VectorXd::Zero((int)tail->frozen.size());
        tail->frozenVel(0) = v0.x();
        tail->frozenVel(1) = v0.y();
        tail->frozenVel(2) = v1.x();
        tail->frozenVel(3) = v1.y();
    }

    tail->X(0, 0) = clamp0.x(); tail->X(0, 1) = clamp0.y();
    tail->X(1, 0) = clamp1.x(); tail->X(1, 1) = clamp1.y();
    tail->V(0, 0) = v0.x();     tail->V(0, 1) = v0.y();
    tail->V(1, 0) = v1.x();     tail->V(1, 1) = v1.y();
}

void ElasticTadpole::sampleHead(MatrixXd& pts, MatrixXd& r2body) const
{
    pts.resize(nSampleHead, 2);
    r2body.resize(nSampleHead, 2);
    for (int i = 0; i < nSampleHead; ++i) {
        double a = 2.0 * M_PI * (double)i / nSampleHead;
        double dx = rHead * std::cos(a);
        double dy = rHead * std::sin(a);
        pts(i, 0) = x + dx;
        pts(i, 1) = y + dy;
        r2body(i, 0) = dx;
        r2body(i, 1) = dy;
    }
}

Vector2d ElasticTadpole::tailTip() const
{
    if (!tail || tail->X.rows() == 0) return Vector2d(x, y);
    int N = (int)tail->X.rows() - 1;
    return Vector2d(tail->X(N, 0), tail->X(N, 1));
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
namespace {
inline bool sampleVel(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                      const VectorXd& uF, const VectorXd& vF,
                      const MeshLocator& loc, double px, double py,
                      int& hint, double& u, double& v)
{
    double lam[3];
    int t = loc.locate(mesh, px, py, hint, lam);
    if (t < 0) { u = 0.0; v = 0.0; return false; }
    hint = t;
    Vector3d L(lam[0], lam[1], lam[2]);
    RowVectorXd phi = fem.computeBasisValue_all(L.transpose()).row(0);
    int locDof = fem.locDof;
    u = 0.0; v = 0.0;
    for (int i = 0; i < locDof; ++i) {
        int g = elem2dof(t, i);
        u += phi(i) * uF(g);
        v += phi(i) * vF(g);
    }
    return true;
}
} // anon

// ---------------------------------------------------------------------------
// step
// ---------------------------------------------------------------------------

bool ElasticTadpole::step(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                          const VectorXd& uField, const VectorXd& vField,
                          const MeshLocator& loc,
                          std::vector<int>& headHint,
                          std::vector<int>& tailHint,
                          double dt)
{
    if (mass <= 0.0 || !tail) init();
    const double xOld = x, yOld = y, thetaOld = theta;
    const double vxOld = vx, vyOld = vy, omegaOld = omega;

    // ====================================================================
    // (1) HEAD: collect drag + swim force, semi-implicit step.
    // ====================================================================
    MatrixXd hpts, hr;
    sampleHead(hpts, hr);
    int Mh = (int)hpts.rows();
    if ((int)headHint.size() != Mh) headHint.assign(Mh, -1);
    double dsHead = (2.0 * M_PI * rHead) / nSampleHead;

    double Fx = 0.0, Fy = 0.0, Tau = 0.0;
    int aliveH = 0;
    for (int i = 0; i < Mh; ++i) {
        double vx_pt = vx - omega * hr(i, 1);
        double vy_pt = vy + omega * hr(i, 0);
        double uF, vF;
        if (!sampleVel(fem, mesh, elem2dof, uField, vField, loc,
                       hpts(i, 0), hpts(i, 1), headHint[i], uF, vF))
            continue;
        ++aliveH;
        double Fxi = -cDragHead * (vx_pt - uF) * dsHead;
        double Fyi = -cDragHead * (vy_pt - vF) * dsHead;
        Fx += Fxi; Fy += Fyi;
        Tau += hr(i, 0) * Fyi - hr(i, 1) * Fxi;
    }
    if (aliveH == 0) return false;

    // Optional station-keeping force toward local low |u| (wake-refuge effect).
    if (kRefuge != 0.0) {
        const double eps = 0.1;
        double uE, vE, uW, vW, uN, vN, uS, vS;
        int hE = -1, hW = -1, hN = -1, hS = -1;
        bool oE = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x + eps, y, hE, uE, vE);
        bool oW = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x - eps, y, hW, uW, vW);
        bool oN = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x, y + eps, hN, uN, vN);
        bool oS = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x, y - eps, hS, uS, vS);
        if (oE && oW && oN && oS) {
            double gx = ((uE * uE + vE * vE) - (uW * uW + vW * vW)) / (2.0 * eps);
            double gy = ((uN * uN + vN * vN) - (uS * uS + vS * vS)) / (2.0 * eps);
            Fx += -kRefuge * gx;
            Fy += -kRefuge * gy;
        }
    }

    // Active upstream swim (in -x direction).
    Fx += -swimForce;

    // Semi-implicit damping update for head.
    double inv_lin = 1.0 / (1.0 + dt * dampLin);
    double inv_ang = 1.0 / (1.0 + dt * dampAng);
    vx = (vx + dt * Fx / mass) * inv_lin;
    vy = (vy + dt * Fy / mass) * inv_lin;
    omega = (omega + dt * Tau / Inertia) * inv_ang;

    // Speed cap.
    double sp = std::hypot(vx, vy);
    if (sp > maxSpeed) { double f = maxSpeed / sp; vx *= f; vy *= f; }
    double omegaCap = std::max(0.0, maxOmega);
    if (omegaCap > 0.0 && std::abs(omega) > omegaCap)
        omega = (omega > 0 ? omegaCap : -omegaCap);
    double alphaCap = std::max(0.0, maxOmegaAccel);
    if (alphaCap > 0.0) {
        double maxDelta = alphaCap * dt;
        double dOmega = omega - omegaOld;
        if (std::abs(dOmega) > maxDelta)
            omega = omegaOld + (dOmega > 0.0 ? maxDelta : -maxDelta);
    }

    // Translate head pose.
    x += dt * vx;
    y += dt * vy;
    theta += dt * omega;

    // Carry the entire tail through the head's rigid-body pose increment
    // before solving elastic/fluid deformation.  Without this, a translating
    // head pulls only the clamped nodes and injects artificial axial strain.
    if (tail) {
        double c = std::cos(theta - thetaOld);
        double s = std::sin(theta - thetaOld);
        Vector2d C0(xOld, yOld), C1(x, y);
        for (int k = 0; k <= tail->N; ++k) {
            Vector2d Xold(tail->X(k, 0), tail->X(k, 1));
            Vector2d r0 = Xold - C0;
            Vector2d bodyV0(vxOld - omegaOld * r0.y(),
                            vyOld + omegaOld * r0.x());
            Vector2d relV(tail->V(k, 0) - bodyV0.x(),
                          tail->V(k, 1) - bodyV0.y());
            Vector2d r1(c * r0.x() - s * r0.y(),
                        s * r0.x() + c * r0.y());
            Vector2d relV1(c * relV.x() - s * relV.y(),
                           s * relV.x() + c * relV.y());
            Vector2d Xnew = C1 + r1;
            Vector2d bodyV1(vx - omega * r1.y(),
                            vy + omega * r1.x());
            Vector2d Vnew = bodyV1 + relV1;
            tail->X(k, 0) = Xnew.x();
            tail->X(k, 1) = Xnew.y();
            tail->V(k, 0) = Vnew.x();
            tail->V(k, 1) = Vnew.y();
        }
    }

    // ====================================================================
    // (2) TAIL: clamp the first two nodes onto the head's rear surface,
    // step the rod against the live fluid using implicit per-node drag.
    // The rod tracks the head through its frozen-DOF endpoints; we don't
    // need to manually rigid-body-drag the whole rod -- the elastic +
    // damping forces are well-behaved as long as the clamp moves smoothly.
    // ====================================================================
    syncTailClamp();

    // Sample the fluid at every tail node, fill the implicit drag info.
    // The elastic tail is a real immersed boundary: the rod feels the local
    // fluid velocity, not an artificial blend with the head velocity.  The
    // matching reaction on the fluid is applied by ns_tadpole_elastic_main.
    int Nt = (int)tail->X.rows();
    if ((int)tailHint.size() != Nt) tailHint.assign(Nt, -1);
    tail->dragCoef.setZero(Nt);
    tail->dragRef.setZero(Nt, 2);
    Eigen::VectorXd ds = rodArclengthWeights(*tail);
    for (int k = 2; k < Nt; ++k) {
        double xk = tail->X(k, 0), yk = tail->X(k, 1);
        double uF, vF;
        if (!sampleVel(fem, mesh, elem2dof, uField, vField, loc,
                       xk, yk, tailHint[k], uF, vF))
            continue;
        tail->dragCoef(k) = cDragTail * ds(k);
        tail->dragRef(k, 0) = uF;
        tail->dragRef(k, 1) = vF;
    }
    tail->Fext.setZero();

    // Advance the elastic tail.
    MatrixXd Xpre = tail->X;
    MatrixXd Vpre = tail->V;
    if (!tail->step(dt)) return false;
    syncTailClamp();

    // Treat the tadpole tail as nearly inextensible.  The Cosserat stretch
    // energy controls Newton's solve, then this arclength projection removes
    // the residual partitioned-IB axial drift while leaving bending free.
    for (int i = 1; i < tail->N; ++i) {
        Vector2d prev(tail->X(i, 0), tail->X(i, 1));
        Vector2d cur(tail->X(i + 1, 0), tail->X(i + 1, 1));
        Vector2d e = cur - prev;
        double len = e.norm();
        if (len < 1e-14) {
            Vector2d rest(tail->X0(i + 1, 0) - tail->X0(i, 0),
                          tail->X0(i + 1, 1) - tail->X0(i, 1));
            len = rest.norm();
            e = (len > 1e-14) ? rest : Vector2d(1.0, 0.0);
        }
        Vector2d corr = prev + (tail->l0(i) / e.norm()) * e;
        tail->X(i + 1, 0) = corr.x();
        tail->X(i + 1, 1) = corr.y();
    }
    syncTailClamp();
    for (int k = 2; k <= tail->N; ++k) {
        tail->V(k, 0) = (tail->X(k, 0) - Xpre(k, 0)) / dt;
        tail->V(k, 1) = (tail->X(k, 1) - Xpre(k, 1)) / dt;
        tail->A(k, 0) = (tail->V(k, 0) - Vpre(k, 0)) / dt;
        tail->A(k, 1) = (tail->V(k, 1) - Vpre(k, 1)) / dt;
    }
    return true;
}

} // namespace ns
