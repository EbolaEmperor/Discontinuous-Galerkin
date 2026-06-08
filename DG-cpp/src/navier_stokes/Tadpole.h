#ifndef TADPOLE_H
#define TADPOLE_H

// ===========================================================================
// 2-D tadpole with a RIGID circular head and an ELASTIC tail.
//
//   * Head : circular rigid body of radius r_head with 3 DOFs (xc, yc, theta).
//            Subject to per-point linear-Stokes drag against the DG flow
//            field, an upstream-pointing self-propulsion ("swim") force, and
//            mass-proportional structural damping.  Stepped semi-implicitly.
//   * Tail : a Cosserat elastic rod whose first two nodes are clamped at the
//            head's rear surface (so the tail tracks the head as a body but
//            still bends under the flow).  Stepped by Newmark-beta with
//            implicit per-node fluid drag (the same trick used for the
//            filament FSI demo, so the rod doesn't blow up under heavy drag).
//
// The clamped-end position and orientation are reset every step from the
// head's current pose, so the tail effectively follows the head while
// remaining free to flap.  The reaction force / torque the tail exerts on
// its clamped end is fed back into the head as an external load, giving a
// realistic tail-flap-drives-head-yaw coupling.
//
// One-way fluid coupling: neither the head nor the tail pushes back on the
// fluid (the tadpole is small compared with the cylinder).
// ===========================================================================

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "CosseratFilament.h"

class FEM;
class Mesh;

namespace ns {

class MeshLocator;

struct Tadpole {
    // ----- head (rigid) geometry / inertia -----
    double rHead = 0.15;         // head radius (in D = 1 units)
    double Ltail = 0.6;          // rest tail length
    double rhoHead = 1.0;        // head density (per unit area)
    int    nSampleHead = 8;      // perimeter samples on the head for drag

    // ----- tail (elastic) parameters -----
    int    Ntail   = 24;         // number of segments in the elastic tail
    double rhoTail = 1.0;        // tail line density (rho_s * h)
    double KB      = 0.02;       // bending stiffness EI
    double KS      = 1e3;        // axial stiffness EA (nearly inextensible)
    double tailDamp = 0.05;      // tail structural Rayleigh-mass damping
    double cDragTail = 6.0;      // per-unit-arclength linear-Stokes drag on the tail

    // ----- head dynamics -----
    double cDragHead = 6.0;      // per-unit-arclength linear-Stokes drag on the head
    double swimForce = 0.0;      // upstream (in -x) self-propulsion force magnitude
    double dampLin   = 0.4;      // mass-proportional translational damping
    double dampAng   = 1.0;      // angular damping
    double maxSpeed  = 1.5;      // safety speed cap

    // ----- head state (the only translational/rotational DOFs) -----
    double x = 0.0, y = 0.0;     // head centre
    double theta = 0.0;          // body axis: tail trails opposite this heading
    double vx = 0.0, vy = 0.0;   // head linear velocity
    double omega = 0.0;          // head angular velocity

    // ----- elastic tail (created by initTail()) -----
    std::unique_ptr<CosseratFilament> tail;

    // ----- cached head mass / inertia -----
    double mass = 0.0;
    double Inertia = 0.0;

    // Build the head inertia and the tail Cosserat rod (straight, in the
    // current world frame at theta).  Call once after setting parameters.
    void init();

    // Sample points on the head perimeter (for collecting head drag).
    //   pts : (nSampleHead x 2) world positions
    //   r2body : (nSampleHead x 2) offset from head centre.
    void sampleHead(Eigen::MatrixXd& pts, Eigen::MatrixXd& r2body) const;

    // World-space tail tip (last node of the elastic rod) and head centre.
    Eigen::Vector2d headCentre() const { return Eigen::Vector2d(x, y); }
    Eigen::Vector2d tailTip() const;

    // Advance head + tail by dt against the live DG flow field.  Returns
    // false on a degenerate sample.
    bool step(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
              const Eigen::VectorXd& uField, const Eigen::VectorXd& vField,
              const MeshLocator& loc,
              std::vector<int>& headHint,
              std::vector<int>& tailHint,
              double dt);
};

} // namespace ns

#endif
