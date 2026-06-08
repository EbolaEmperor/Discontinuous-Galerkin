#include "Tadpole.h"

#include "FEM.h"
#include "Mesh.h"
#include "NavierStokes.h"      // MeshLocator + the basis-evaluation helpers

#include <cmath>
#include <iostream>

namespace ns {

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::RowVectorXd;

void Tadpole::initInertia()
{
    // Head: a disk of radius rHead, mass = rho * pi * r^2,
    // inertia about centre = 0.5 * m * r^2.
    double m_h = rho * M_PI * rHead * rHead;
    double I_h = 0.5 * m_h * rHead * rHead;
    // Tail: a thin rod of length L_tail and mass = rho * L_tail * (small width
    // proxy = rHead * 0.5), so the tail mass is rho * L_tail * 0.5 * rHead.
    // Inertia about head centre is m_t * (L/2)^2 + (1/12) m_t L^2 from
    // parallel axis -> (1/3) m_t L^2.
    double m_t = rho * Ltail * 0.5 * rHead;
    double I_t = (1.0 / 3.0) * m_t * Ltail * Ltail;

    mass = m_h + m_t;
    Inertia = I_h + I_t;
}

void Tadpole::sampleBody(MatrixXd& pts, MatrixXd& r2body) const
{
    int M = nSampleHead + nSampleTail;
    pts.resize(M, 2);
    r2body.resize(M, 2);
    // Head: equispaced on the perimeter.
    for (int i = 0; i < nSampleHead; ++i) {
        double a = 2.0 * M_PI * (double)i / nSampleHead;
        double dx = rHead * std::cos(a);
        double dy = rHead * std::sin(a);
        pts(i, 0) = x + dx;
        pts(i, 1) = y + dy;
        r2body(i, 0) = dx;
        r2body(i, 1) = dy;
    }
    // Tail: starts at head centre going OPPOSITE the heading direction.
    Vector2d tailDir(-std::cos(theta), -std::sin(theta));
    for (int i = 0; i < nSampleTail; ++i) {
        double s = (i + 0.5) / nSampleTail * Ltail;
        double dx = s * tailDir.x();
        double dy = s * tailDir.y();
        pts(nSampleHead + i, 0) = x + dx;
        pts(nSampleHead + i, 1) = y + dy;
        r2body(nSampleHead + i, 0) = dx;
        r2body(nSampleHead + i, 1) = dy;
    }
}

namespace {
// Sample (u, v) at world point (px, py) using the cached element hint.  Returns
// false if the point is outside the mesh; the velocity is set to zero in that
// case so the caller can decide whether to drop it.
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

bool Tadpole::step(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                   const VectorXd& uField, const VectorXd& vField,
                   const MeshLocator& loc, std::vector<int>& hint, double dt)
{
    if (mass <= 0.0 || Inertia <= 0.0) initInertia();

    MatrixXd pts, r2body;
    sampleBody(pts, r2body);
    int M = (int)pts.rows();
    if ((int)hint.size() != M) hint.assign(M, -1);

    // Per-sample arclength weight: head perimeter / nSampleHead, tail length /
    // nSampleTail.  Forces will be linear-stokes drag * weight.
    double dsHead = (2.0 * M_PI * rHead) / nSampleHead;
    double dsTail = Ltail / nSampleTail;

    // ---- 1) collect drag force + drag torque ----
    double Fx_drag = 0.0, Fy_drag = 0.0, Tau_drag = 0.0;
    int alive = 0;
    for (int i = 0; i < M; ++i) {
        // Local rigid-body velocity at this sample point
        double vx_pt = vx - omega * r2body(i, 1);
        double vy_pt = vy + omega * r2body(i, 0);
        double uF, vF;
        bool ok = sampleVel(fem, mesh, elem2dof, uField, vField, loc,
                            pts(i, 0), pts(i, 1), hint[i], uF, vF);
        if (!ok) continue;
        ++alive;
        double w = (i < nSampleHead) ? dsHead : dsTail;
        // Linear-Stokes drag on the body:  F = -c * (V_body - u_fluid) * ds
        double Fx_i = -cDrag * (vx_pt - uF) * w;
        double Fy_i = -cDrag * (vy_pt - vF) * w;
        Fx_drag += Fx_i;
        Fy_drag += Fy_i;
        Tau_drag += r2body(i, 0) * Fy_i - r2body(i, 1) * Fx_i;
    }
    if (alive == 0) return false;   // tadpole completely outside the mesh

    // ---- 2) station-keeping force toward local low |u| (wake-refuge effect) ----
    // F_refuge = -kRefuge * grad |u|^2 evaluated at the head centre.
    // We approximate grad |u|^2 by central finite differences of the speed
    // squared at the head centre.
    double Fx_refuge = 0.0, Fy_refuge = 0.0;
    {
        const double eps = 0.1;     // FD offset; ~one element size
        double uE, vE, uW, vW, uN, vN, uS, vS;
        int hE = -1, hW = -1, hN = -1, hS = -1;
        bool oE = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x + eps, y, hE, uE, vE);
        bool oW = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x - eps, y, hW, uW, vW);
        bool oN = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x, y + eps, hN, uN, vN);
        bool oS = sampleVel(fem, mesh, elem2dof, uField, vField, loc, x, y - eps, hS, uS, vS);
        if (oE && oW && oN && oS) {
            double sE = uE * uE + vE * vE;
            double sW = uW * uW + vW * vW;
            double sN = uN * uN + vN * vN;
            double sS = uS * uS + vS * vS;
            double gx = (sE - sW) / (2.0 * eps);
            double gy = (sN - sS) / (2.0 * eps);
            Fx_refuge = -kRefuge * gx;
            Fy_refuge = -kRefuge * gy;
        }
    }

    // ---- 3) total force / torque + structural damping ----
    // Active upstream swim (pointing in -x, like a fish facing the flow).
    Fx_drag += -swimForce;
    double Fx_tot = Fx_drag + Fx_refuge - dampLin * mass * vx;
    double Fy_tot = Fy_drag + Fy_refuge - dampLin * mass * vy;
    double Tau_tot = Tau_drag - dampAng * Inertia * omega;

    // ---- 4) semi-implicit time step (backward-Euler on damping for safety) ----
    //    a = (F - dampLin*m*v_new) / m  =>  v_new = (m*v + dt*F) / (m + dt*dampLin*m)
    // Equivalent to:
    //    v_new = (v + dt*F/m) / (1 + dt*dampLin)
    double inv_lin = 1.0 / (1.0 + dt * dampLin);
    double inv_ang = 1.0 / (1.0 + dt * dampAng);
    vx = (vx + dt * (Fx_drag + Fx_refuge) / mass) * inv_lin;
    vy = (vy + dt * (Fy_drag + Fy_refuge) / mass) * inv_lin;
    omega = (omega + dt * Tau_drag / Inertia) * inv_ang;
    (void)Fx_tot; (void)Fy_tot; (void)Tau_tot;

    // Speed cap.
    double sp = std::hypot(vx, vy);
    if (sp > maxSpeed) {
        double f = maxSpeed / sp;
        vx *= f; vy *= f;
    }
    if (std::abs(omega) > 5.0) omega = (omega > 0 ? 5.0 : -5.0);

    x += dt * vx;
    y += dt * vy;
    theta += dt * omega;
    return true;
}

} // namespace ns
