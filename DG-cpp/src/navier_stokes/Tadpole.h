#ifndef TADPOLE_H
#define TADPOLE_H

// ===========================================================================
// 2-D rigid tadpole: a circular head of radius r_head with centre X_h, plus
// a straight rigid tail of length L from X_h pointing along the body axis at
// angle theta.  Three degrees of freedom in total: x, y of head centre, and
// orientation theta.
//
// The tadpole drifts under
//   1. linear Stokes-style drag at sample points along the body, against the
//      DG flow field;
//   2. a weak "station-keeping" force toward local low-speed regions, which
//      reproduces the empirical observation that fish station-keep in the
//      wake refuge of bluff bodies (Liao 2003, Karman gait).  Without this
//      term a passive rigid body in a uniform stream just gets pushed
//      downstream.
//
// Time integration is semi-implicit:  the fluid drag is collected at the
// sample points, summed and torqued onto the rigid body, and then the
// translational + angular equations of motion are stepped explicitly with a
// strong implicit damping term so the integration is unconditionally stable
// even with stiff drag coefficients.
//
// Collisions with the cylinder body and the rectangular outer walls are
// handled by an external projection routine in the driver, similar to the
// filament demo.
// ===========================================================================

#include <Eigen/Dense>
#include <vector>

class FEM;
class Mesh;

namespace ns {

class MeshLocator;

struct Tadpole {
    // Geometry / inertia
    double rHead = 0.15;         // head radius (in D = 1 units)
    double Ltail = 0.6;          // tail length
    double rho   = 1.0;          // rigid-body density (per unit area)
    int    nSampleHead = 8;      // perimeter samples on the head for drag
    int    nSampleTail = 12;     // arclength samples on the tail for drag

    // Drag / station-keeping coefficients
    double cDrag    = 6.0;       // linear drag coefficient per unit arclength
    double kRefuge  = 6.0;       // strength of  -k * grad |u|^2  (toward low-speed regions)
    double swimForce = 0.0;      // upstream-pointing self-propulsion force magnitude
                                 // (units of total force on the body); positive
                                 // values point in -x direction.  Modelling the
                                 // fish actively swimming against the flow.
    double dampLin  = 0.4;       // structural mass-proportional damping (head)
    double dampAng  = 1.0;       // structural angular damping
    double maxSpeed = 1.5;       // safety speed cap (in U_inf units)

    // State (the only DOFs)
    double x = 0.0, y = 0.0;     // head centre
    double theta = 0.0;          // orientation (tail points along (cos t, sin t))
    double vx = 0.0, vy = 0.0;   // head linear velocity
    double omega = 0.0;          // angular velocity

    // Cached mass / inertia (computed by initInertia()).
    double mass = 0.0;
    double Inertia = 0.0;

    // Compute lumped mass + moment of inertia from rho and geometry.
    void initInertia();

    // Build sample points along the body in world coordinates.  Returns
    //   pts : (M x 2) world positions
    //   r2body : (M x 2) offset from head centre (used to translate the drag
    //            force into a torque about (x, y)).
    void sampleBody(Eigen::MatrixXd& pts, Eigen::MatrixXd& r2body) const;

    // Step the rigid body forward by dt with hydrodynamic forces sampled
    // from the supplied DG flow field.  Returns false on a degenerate sample.
    bool step(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
              const Eigen::VectorXd& uField, const Eigen::VectorXd& vField,
              const MeshLocator& loc, std::vector<int>& hint, double dt);

    // World-space tail tip:  tail starts AT the head centre and ends at the
    // tip; the segment is drawn as a thick polyline by the renderer.
    inline Eigen::Vector2d headCentre() const { return Eigen::Vector2d(x, y); }
    inline Eigen::Vector2d tailTip() const {
        return Eigen::Vector2d(x + Ltail * std::cos(theta + M_PI),
                               y + Ltail * std::sin(theta + M_PI));
    }
};

} // namespace ns

#endif
