#include "Tadpole.h"

#include "FEM.h"
#include "Mesh.h"
#include "NavierStokes.h"      // MeshLocator
#include "IBCoupler.h"         // meshToRod, rodArclengthWeights

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

void Tadpole::init()
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
}

void Tadpole::sampleHead(MatrixXd& pts, MatrixXd& r2body) const
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

Vector2d Tadpole::tailTip() const
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

bool Tadpole::step(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                   const VectorXd& uField, const VectorXd& vField,
                   const MeshLocator& loc,
                   std::vector<int>& headHint,
                   std::vector<int>& tailHint,
                   double dt)
{
    if (mass <= 0.0 || !tail) init();

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
    if (std::abs(omega) > 5.0) omega = (omega > 0 ? 5.0 : -5.0);

    // Translate head pose.
    x += dt * vx;
    y += dt * vy;
    theta += dt * omega;

    // ====================================================================
    // (2) TAIL: clamp the first two nodes onto the head's rear, then step
    // the rod against the live fluid using implicit per-node drag.
    // ====================================================================
    Vector2d tan(-std::cos(theta), -std::sin(theta));     // body-axis backward
    Vector2d clamp0 = Vector2d(x, y) + rHead * tan;
    // Use the original rest-edge length for node 1 so the wall angle stays
    // consistent with the rest configuration; this matches what clampRoot
    // expects (node 0 and node 1 are both frozen).
    double l0 = (Ltail / Ntail);
    Vector2d clamp1 = clamp0 + l0 * tan;
    // Update the frozen DOF target values (positions of nodes 0 and 1).
    if ((int)tail->frozen.size() >= 4) {
        tail->frozenVal(0) = clamp0.x();
        tail->frozenVal(1) = clamp0.y();
        tail->frozenVal(2) = clamp1.x();
        tail->frozenVal(3) = clamp1.y();
    }

    // Sample the fluid at every tail node, fill the implicit drag info.
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
    if (!tail->step(dt)) return false;

    return true;
}

} // namespace ns
