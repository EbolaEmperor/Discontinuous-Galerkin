#ifndef BLAST_ROD_GEOMETRY_H
#define BLAST_ROD_GEOMETRY_H

#include "Core.h"

namespace euler_ale {

struct BlastRodGeom {
    double xa = 0.0;
    double xb = 0.975;
    double ya = 0.0;
    double yb = 0.65;
    double rodX = 0.32;
    double rodBaseY = 0.0;
    double rodW = 0.045;
    double rodL = 0.42;

    double rodLeft() const;
    double rodRight() const;
    double rodTipY() const;
};

struct BlastRodMap {
    BlastRodGeom geom;
    double t0 = 0.0;
    double t1 = 1.0;
    double q0 = 0.0;
    double q1 = 0.0;
    double normalDecay = 0.07;
    double wakeDecay = 0.20;
    double wallMargin = 0.075;

    static double smooth01(double s);
    double q(double t) const;
    double qd(double t) const;
    double rodMode(double y) const;
    double rodModeDy(double y) const;
    double influence(double x, double y) const;
    double centerX(double y, double t) const;
    Vector2d beamCenter(double y, double t) const;
    Vector2d beamTangent(double y, double t) const;
    Vector2d beamNormal(double y, double t) const;
    Vector2d beamPoint(double y, double r, double t) const;
    Vector2d beamDq(double y, double r, double t) const;
    Vector2d beamExtendedPoint(double x, double y, double t) const;
    Vector2d beamExtendedDq(double x, double y, double t) const;
    Vector2d beamDqAtPhysical(const Vector2d& p, double t) const;
    Vector2d refToPhys(const Vector2d& X, double t) const;
    Vector2d velocityAt(double x, double y, double t) const;
    void setLinearMotion(double a, double b, double qa, double qb);
    double maxMeshSpeed(double t) const;
};

struct BlastRodParams {
    double mass = 0.020;
    double stiffness = 14.0;
    double damping = 0.060;
    double pExt = 1.0;
    double qMin = -0.035;
    double qMax = 0.22;
};

double blastRodSmoothRamp(double s);

} // namespace euler_ale

#endif
