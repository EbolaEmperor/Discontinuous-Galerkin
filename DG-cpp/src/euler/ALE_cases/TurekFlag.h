#ifndef TUREK_FLAG_H
#define TUREK_FLAG_H

#include "Core.h"

namespace euler_ale {

struct TurekFlagGeom {
    double xa = 0.0;
    double xb = 2.5;
    double ya = 0.0;
    double yb = 0.41;
    double cx = 0.20;
    double cy = 0.20;
    double r = 0.05;
    double flagX0 = 0.25;
    double flagL = 0.35;
    double flagT = 0.02;

    double flagX1() const;
    double flagY0() const;
    double flagY1() const;
};

struct TurekFlagMap {
    TurekFlagGeom geom;
    double t0 = 0.0;
    double t1 = 1.0;
    double q0 = 0.0;
    double q1 = 0.0;
    double wakeDecay = 0.22;
    double normalDecay = 0.10;
    double wallMargin = 0.055;

    static double smooth01(double s);
    double q(double t) const;
    double qd(double t) const;
    double mode(double x) const;
    double influence(double x, double y) const;
    double flagMode(double x) const;
    double flagY(double x, double yRef, double t) const;
    Vector2d refToPhys(const Vector2d& X, double t) const;
    Vector2d velocityAt(double x, double y, double t) const;
    void setLinearMotion(double a, double b, double qa, double qb);
    double maxMeshSpeed(double t) const;
};

struct FlagModalParams {
    double mass = 0.018;
    double stiffness = 42.0;
    double damping = 0.045;
    double pExt = 1.0;
    double qMin = -0.035;
    double qMax = 0.035;
};

double signedBox(double x, double y, double x0, double x1, double y0, double y1);
double turekFluidSdf(const TurekFlagGeom& geom, double x, double y);
ModalState advanceFlagSymplectic(ModalState state, double fluidForce,
                                 const FlagModalParams& params, double dt);
Mesh makeTurekFlagMesh(const TurekFlagGeom& geom, double h, int maxIter, bool verbose);
int turekBoundaryTag(double x, double y, double t, const TurekFlagMap& map);
double flagGeneralizedForce(const Space& sp, const MatrixXd& U, const TurekFlagMap& map,
                            double time, double pExt, double* meanPressure = nullptr,
                            double* lift = nullptr, double* drag = nullptr);

int runTurek(bool quick);

} // namespace euler_ale

#endif
