#ifndef SPOON_H
#define SPOON_H

// ===========================================================================
// 2-D rigid "spoon" blade for the stirred-soup FSI demo.
//
// The blade is an ellipse (semi-axis aRad along its long, radial axis; semi-
// axis bThk along its thin axis) submerged in the bowl and swept rigidly about
// a fixed pivot (the hand) through ONE arc-stroke, then withdrawn.  The blade
// is broadside to its own tangential motion, so it acts as a paddle: each of
// its two tips sheds a shear layer that rolls up into a vortex, and the two
// counter-rotating vortices pair into a self-propelled dipole once the blade
// stops -- the "vortices on both sides of the spoon".  (Cf. the impulsively
// started/stopped flat plate and the self-propelled vortex dipole.)
//
// Coupling to the flow is one-directional *driving*: the blade motion is
// PRESCRIBED (a known function of time), and it is imposed on the fluid by a
// direct-forcing immersed boundary,  F = (alpha/dt) * gain(t) * (V_blade - u),
// distributed at a cloud of Lagrangian markers covering the blade footprint and
// scattered through the DG basis (the consistent FE-immersed-boundary form, the
// same transfer used by the flexible-filament demo).  No back-reaction model is
// needed -- the blade trajectory is given.
//
// Stir schedule (no velocity jumps -> no numerical impulse):
//   [0, tEnter)                immerse: blade fixed at start pose, gain 0 -> 1
//   [tEnter, tEnter+tStroke)   stroke : angle sweeps 0 -> sweep with a QUADRATIC
//                                        (parabolic) angular-velocity profile
//                                        omega(tau) = omega_max * 4 tau (1 - tau),
//                                        tau = (t - tEnter)/tStroke in [0,1] -- it
//                                        starts AND ends at zero, peaking
//                                        omega_max = 1.5*sweep/tStroke at mid-stroke.
//   [.. , .. + tLift)          lift   : blade fixed at end pose, gain 1 -> 0
//   afterwards                 free   : no forcing, blade not drawn
//
// All lengths are in bowl-radius units (R = 1) and speeds in the reference
// stir-speed units used by the driver.
// ===========================================================================

#include <Eigen/Dense>
#include <vector>

class FEM;
class Mesh;

namespace ns {

class MeshLocator;

struct Spoon {
    // ---- blade geometry ----
    // Two shapes: an ellipse, or a CRESCENT (month/scoop) whose concave opening
    // faces the blade's motion direction (+tangential).  The blade rotates
    // rigidly about the pivot, so the opening always tracks the motion.
    bool   crescent = true;     // true = crescent moon (lune), false = ellipse
    // Crescent moon = a lune: the region inside an outer circle (radius R1) AND outside a smaller
    // circle (radius R2 = r2frac*R1) shifted by offset*R1 toward +tangential, so the concave notch
    // (opening) faces the motion direction and the two tips taper to sharp horns at +-aRad radially.
    // R1 is derived so the horns sit at +-aRad.  r2frac->1 / smaller offset => thinner crescent.
    double crescentR2frac = 0.86;  // R2 / R1   (closer to 1 = thinner moon)
    double crescentOffset = 0.50;  // (centre separation) / R1   (larger = deeper notch)
    double aRad  = 0.18;     // half radial extent (ellipse semi-axis / crescent horn radius)
    double bThk  = 0.035;    // half-thickness (ellipse thin semi-axis; crescent band half-thickness)
    double dmark = 0.012;    // Lagrangian marker spacing (keep ~ mesh size h)

    // ---- placement / pivot ----
    double px = 0.0, py = 0.0;   // stir pivot (the "hand"), usually the bowl centre
    double rPivot  = 0.45;       // blade-centre distance from the pivot
    double bearing0 = -0.7;      // initial angular position of the blade about the pivot (rad)

    // ---- stir schedule ----
    double tEnter  = 0.25;   // immersion ramp (gain 0->1, blade still)
    double tStroke = 1.50;   // stroke duration
    double tRamp   = 0.45;   // LEGACY (ignored by the quadratic profile; kept for config compat)
    double tLift   = 0.60;   // withdrawal ramp (gain 1->0, blade still)
    double sweep   = 1.50;   // total swept angle (rad, signed; + = CCW)

    // ---- immersed-boundary forcing ----
    double ibAlpha   = 1.0;  // direct-forcing gain:  F = (ibAlpha/dt)*gain*(V-u)*dA
    double ibForceCap = 0.0; // optional per-marker |F| cap (0 = off)

    // ---- template (filled by build()) ----
    Eigen::MatrixXd r0;      // (M x 2) marker offset from pivot at stroke start (phi=0)
    double wMark = 0.0;      // per-marker area weight (dmark^2)

    // ---- current state (filled by update()) ----
    Eigen::MatrixXd X, V;    // (M x 2) world marker positions + velocities
    double phi = 0.0;        // current swept angle
    double omega = 0.0;      // current angular velocity
    double gain = 0.0;       // current forcing envelope in [0,1]
    std::vector<int> hint;   // per-marker element-locator cache
    double cresEtaShift_ = 0.0;   // et-centroid shift used to centre the lune (set by build())

    int M() const { return (int)r0.rows(); }
    // Peak angular velocity of the parabolic stroke profile (signed).
    double peakOmega() const { return tStroke > 0.0 ? 1.5 * sweep / tStroke : 0.0; }
    double tEndStroke() const { return tEnter + tStroke; }
    double tDone()      const { return tEnter + tStroke + tLift; }
    bool submerged()    const { return gain > 1e-6; }
    // The immersed-boundary constraint is applied while the blade is in the soup:
    // through immersion + the stroke, then the blade is withdrawn (constraint off).
    // At stroke end omega -> 0, so V -> 0 and removing the constraint is smooth.
    bool active(double t) const { return t <= tEndStroke() + 1e-9; }

    // Build the marker template (ellipse grid).  Call once after setting geometry.
    void build();

    // Set phi, omega, gain and the world marker positions/velocities at time t.
    void update(double t);

    // Sample the fluid at every marker, form the direct-forcing load
    //   F_k = (ibAlpha/dt)*gain*(V_k - u_h(X_k)) * dA,
    // and scatter it (through the DG basis) into loadFx/loadFy (added to, not
    // zeroed).  Returns the net force the blade applies to the fluid (for
    // diagnostics; the reaction on the blade is its negative).
    Eigen::Vector2d applyForcing(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
                                 const Eigen::VectorXd& uField, const Eigen::VectorXd& vField,
                                 const MeshLocator& loc, double dt,
                                 Eigen::VectorXd& loadFx, Eigen::VectorXd& loadFy);

    // Crescent-moon (lune) geometry derived from aRad / crescentR2frac / crescentOffset:
    // outer circle C1(centre (0,eta1) in the er/et frame, radius R1), inner C2(centre
    // (0,eta2), radius R2), crescent = inside C1 and outside C2; horns at +-aRad.
    struct CrescentGeom { double R1, R2, eta1, eta2; };
    CrescentGeom crescentGeom() const;

    // Test whether a world-space point (x,y) lies inside the blade at its current pose.
    // Used to mask the blade interior during vorticity rendering.
    bool contains(double x, double y) const;

    // Closed-loop blade outline (nSeg+1 x 2) at the current pose, for rendering.
    Eigen::MatrixXd outline(int nSeg = 64) const;

    // Centre of the blade in world coordinates at the current pose.
    Eigen::Vector2d centre() const;
};

} // namespace ns

#endif
