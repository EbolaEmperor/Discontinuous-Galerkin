#include "Geometry.h"

#include <algorithm>
#include <cmath>

namespace euler_ale {

double blastRodSmoothRamp(double s) {
    s = std::clamp(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
}

double BlastRodGeom::rodLeft() const { return rodX - 0.5 * rodW; }
double BlastRodGeom::rodRight() const { return rodX + 0.5 * rodW; }
double BlastRodGeom::rodRootRadius() const {
    return std::clamp(rodRootR, 0.0, 0.49 * rodW);
}
double BlastRodGeom::rodFootLeft() const { return rodLeft() - rodRootRadius(); }
double BlastRodGeom::rodFootRight() const { return rodRight() + rodRootRadius(); }
double BlastRodGeom::rodTipY() const { return rodBaseY + rodL; }
double BlastRodGeom::flowXb() const {
    double f = std::max(0.0, spongeFraction);
    return xa + (xb - xa) / (1.0 + f);
}
double BlastRodGeom::spongeStartX() const { return flowXb(); }
double BlastRodGeom::spongeWidth() const {
    return std::max(0.0, xb - spongeStartX());
}
double BlastRodGeom::spongeCoordinate(double x) const {
    double width = spongeWidth();
    if (width <= 1e-14) return 0.0;
    return std::clamp((x - spongeStartX()) / width, 0.0, 1.0);
}

double BlastRodMap::smooth01(double s) {
    return blastRodSmoothRamp(s);
}

double BlastRodMap::q(double t) const {
    double a = std::clamp((t - t0) / std::max(1e-14, t1 - t0), 0.0, 1.0);
    return (1.0 - a) * q0 + a * q1;
}

double BlastRodMap::qd(double) const {
    return (q1 - q0) / std::max(1e-14, t1 - t0);
}

double BlastRodMap::rodMode(double y) const {
    double s = std::clamp((y - geom.rodBaseY) / std::max(1e-14, geom.rodL), 0.0, 1.0);
    const double beta = 1.875104068711961;
    const double sigma = (std::cosh(beta) + std::cos(beta)) /
                         (std::sinh(beta) + std::sin(beta));
    auto raw = [&](double z) {
        double bz = beta * z;
        return std::cosh(bz) - std::cos(bz) -
               sigma * (std::sinh(bz) - std::sin(bz));
    };
    return raw(s) / raw(1.0);
}

double BlastRodMap::rodModeDy(double y) const {
    double L = std::max(1e-14, geom.rodL);
    double s = std::clamp((y - geom.rodBaseY) / L, 0.0, 1.0);
    const double beta = 1.875104068711961;
    const double sigma = (std::cosh(beta) + std::cos(beta)) /
                         (std::sinh(beta) + std::sin(beta));
    auto raw = [&](double z) {
        double bz = beta * z;
        return std::cosh(bz) - std::cos(bz) -
               sigma * (std::sinh(bz) - std::sin(bz));
    };
    double bz = beta * s;
    double dRawDs = beta * (std::sinh(bz) + std::sin(bz) -
                            sigma * (std::cosh(bz) - std::cos(bz)));
    return dRawDs / (raw(1.0) * L);
}

double BlastRodMap::influence(double x, double y) const {
    double yTip = geom.rodTipY();
    double xBand = std::max(0.0, std::abs(x - geom.rodX) - 0.5 * geom.rodW);
    double yBand = 0.0;
    if (y < geom.rodBaseY) yBand = geom.rodBaseY - y;
    else if (y > yTip) yBand = y - yTip;
    double nearRod = std::exp(-(xBand / normalDecay) * (xBand / normalDecay)) *
                     std::exp(-(yBand / normalDecay) * (yBand / normalDecay));

    double xWall = std::min(x - geom.xa, geom.xb - x);
    double yWallTop = geom.yb - y;
    double wall = smooth01(xWall / wallMargin) * smooth01(yWallTop / wallMargin);
    return nearRod * wall;
}

double BlastRodMap::centerX(double y, double t) const {
    return geom.rodX + q(t) * rodMode(y);
}

Vector2d BlastRodMap::beamCenter(double y, double t) const {
    double yc = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    return Vector2d(centerX(yc, t), yc);
}

Vector2d BlastRodMap::beamTangent(double y, double t) const {
    double slope = q(t) * rodModeDy(y);
    double inv = 1.0 / std::sqrt(1.0 + slope * slope);
    return Vector2d(slope * inv, inv);
}

Vector2d BlastRodMap::beamNormal(double y, double t) const {
    double slope = q(t) * rodModeDy(y);
    double inv = 1.0 / std::sqrt(1.0 + slope * slope);
    return Vector2d(inv, -slope * inv);
}

Vector2d BlastRodMap::beamPoint(double y, double r, double t) const {
    double yc = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    return beamCenter(yc, t) + r * beamNormal(yc, t);
}

Vector2d BlastRodMap::beamDq(double y, double r, double t) const {
    double yc = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    double mode = rodMode(yc);
    double dmode = rodModeDy(yc);
    double slope = q(t) * dmode;
    double inv = 1.0 / std::sqrt(1.0 + slope * slope);
    double inv3 = inv * inv * inv;
    Vector2d dNdq(-slope * inv3 * dmode, -inv3 * dmode);
    return Vector2d(mode, 0.0) + r * dNdq;
}

Vector2d BlastRodMap::beamExtendedPoint(double x, double y, double t) const {
    double yc = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    double r = x - geom.rodX;
    double axial = y - yc;
    return beamPoint(yc, r, t) + axial * beamTangent(yc, t);
}

Vector2d BlastRodMap::beamExtendedDq(double x, double y, double t) const {
    double yc = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    double r = x - geom.rodX;
    double axial = y - yc;
    double dmode = rodModeDy(yc);
    double slope = q(t) * dmode;
    double inv = 1.0 / std::sqrt(1.0 + slope * slope);
    double inv3 = inv * inv * inv;
    Vector2d dTdq(inv3 * dmode, -slope * inv3 * dmode);
    return beamDq(yc, r, t) + axial * dTdq;
}

Vector2d BlastRodMap::beamDqAtPhysical(const Vector2d& p, double t) const {
    double y = std::clamp(p.y(), geom.rodBaseY, geom.rodTipY());
    for (int it = 0; it < 5; ++it) {
        Vector2d c = beamCenter(y, t);
        Vector2d cy(q(t) * rodModeDy(y), 1.0);
        double den = std::max(cy.squaredNorm(), 1e-14);
        y -= (c - p).dot(cy) / den;
        y = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    }
    Vector2d n = beamNormal(y, t);
    double r = std::clamp((p - beamCenter(y, t)).dot(n),
                          -0.5 * geom.rodW, 0.5 * geom.rodW);
    return beamDq(y, r, t);
}

Vector2d BlastRodMap::refToPhys(const Vector2d& X, double t) const {
    Vector2d beamMapped = beamExtendedPoint(X.x(), X.y(), t);
    Vector2d exactDisp = beamMapped - X;
    return X + influence(X.x(), X.y()) * exactDisp;
}

Vector2d BlastRodMap::velocityAt(double x, double y, double t) const {
    double yc = std::clamp(y, geom.rodBaseY, geom.rodTipY());
    Vector2d c = beamCenter(yc, t);
    Vector2d n = beamNormal(yc, t);
    double rRaw = (x - c.x()) / std::max(n.x(), 1e-8);
    double dn = std::max(0.0, std::abs(rRaw) - 0.5 * geom.rodW);
    double nearBeam = std::exp(-(dn / normalDecay) * (dn / normalDecay));
    double xWall = std::min(x - geom.xa, geom.xb - x);
    double yWallTop = geom.yb - y;
    double wall = smooth01(xWall / wallMargin) * smooth01(yWallTop / wallMargin);
    return qd(t) * nearBeam * wall * beamExtendedDq(x, y, t);
}

void BlastRodMap::setLinearMotion(double a, double b, double qa, double qb) {
    t0 = a;
    t1 = b;
    q0 = qa;
    q1 = qb;
}

double BlastRodMap::maxMeshSpeed(double t) const { return std::abs(qd(t)); }

} // namespace euler_ale
