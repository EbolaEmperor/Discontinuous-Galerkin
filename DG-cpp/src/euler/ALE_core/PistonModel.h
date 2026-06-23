#ifndef PISTON_MODEL_H
#define PISTON_MODEL_H

#include "Core.h"

namespace euler_ale {

struct PistonMap {
    double amp = 0.06;
    double period = 1.0;
    bool prescribed = true;
    double t0 = 0.0;
    double t1 = 1.0;
    double s0 = 0.0;
    double s1 = 0.0;

    double s(double t) const;
    double sd(double t) const;
    Vector2d refToPhys(const Vector2d& X, double t) const;
    Vector2d velocityAt(double x, double y, double t) const;
    double xLeft(double t) const;
    void setLinearMotion(double a, double b, double xa, double xb);
    double maxMeshSpeed(double t) const;
};

struct PistonState {
    double x = 0.0;
    double v = 0.0;
};

struct FSIParams {
    double mass = 1.0;
    double stiffness = 55.0;
    double damping = 0.9;
    double pExt = 1.0;
    double xMin = -0.075;
    double xMax = 0.075;
};

int pistonBoundaryTag(double x, double y, double t, const PistonMap& map);
double pistonPressureForce(const Space& sp, const MatrixXd& U, double pExt,
                           double* meanPressure = nullptr);
double structureAcceleration(const PistonState& state, double fluidForce,
                             const FSIParams& params);
PistonState advanceStructureSymplectic(PistonState state, double fluidForce,
                                       const FSIParams& params, double dt);

} // namespace euler_ale

#endif
