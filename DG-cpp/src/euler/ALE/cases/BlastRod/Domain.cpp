#include "Domain.h"

#include "BodyFittedMesh.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace euler_ale {
namespace {

std::vector<std::pair<Vector2d, Vector2d>>
blastRodDomainSegments(const BlastRodGeom& g) {
    return {
        {Vector2d(g.xa, g.ya), Vector2d(g.rodFootLeft(), g.ya)},
        {Vector2d(g.rodFootRight(), g.ya), Vector2d(g.xb, g.ya)},
        {Vector2d(g.xb, g.ya), Vector2d(g.xb, g.yb)},
        {Vector2d(g.xb, g.yb), Vector2d(g.xa, g.yb)},
        {Vector2d(g.xa, g.yb), Vector2d(g.xa, g.ya)}
    };
}

std::vector<Vector2d> roundedRodBoundaryLoop(const BlastRodGeom& g, int arcSegments) {
    std::vector<Vector2d> loop;
    double r = g.rodRootRadius();
    auto add = [&](const Vector2d& p) {
        if (loop.empty() || (loop.back() - p).norm() > 1e-12) loop.push_back(p);
    };
    if (r <= 1e-14) {
        add(Vector2d(g.rodLeft(), g.rodBaseY));
        add(Vector2d(g.rodRight(), g.rodBaseY));
        add(Vector2d(g.rodRight(), g.rodTipY()));
        add(Vector2d(g.rodLeft(), g.rodTipY()));
        return loop;
    }

    arcSegments = std::max(4, arcSegments);
    double yArc = g.rodBaseY + r;
    add(Vector2d(g.rodFootLeft(), g.rodBaseY));
    add(Vector2d(g.rodFootRight(), g.rodBaseY));
    Vector2d cRight(g.rodRight() + r, yArc);
    for (int k = 1; k <= arcSegments; ++k) {
        double a = -0.5 * M_PI - 0.5 * M_PI * static_cast<double>(k) / arcSegments;
        add(cRight + r * Vector2d(std::cos(a), std::sin(a)));
    }
    add(Vector2d(g.rodRight(), g.rodTipY()));
    add(Vector2d(g.rodLeft(), g.rodTipY()));
    Vector2d cLeft(g.rodLeft() - r, yArc);
    add(Vector2d(g.rodLeft(), yArc));
    for (int k = 1; k <= arcSegments; ++k) {
        double a = -0.5 * M_PI * static_cast<double>(k) / arcSegments;
        add(cLeft + r * Vector2d(std::cos(a), std::sin(a)));
    }
    return loop;
}

void addRoundedRodFluidBoundarySegments(SolidBodyMeshSpec& spec, const BlastRodGeom& g) {
    double hNear = spec.hNearFactor * spec.hFar;
    int nArc = std::max(6, static_cast<int>(std::ceil(0.5 * M_PI * g.rodRootRadius() /
                                                      std::max(0.55 * hNear, 1e-12))));
    std::vector<Vector2d> loop = roundedRodBoundaryLoop(g, nArc);
    if (loop.size() < 4) return;

    for (int i = 1; i + 1 < static_cast<int>(loop.size()); ++i) {
        spec.fixedSegments.push_back({loop[i], loop[i + 1]});
    }
    spec.fixedSegments.push_back({loop.back(), loop.front()});
}

} // namespace

double blastRodFluidSdf(const BlastRodGeom& g, double x, double y) {
    double dRect = -std::min({x - g.xa, g.xb - x, y - g.ya, g.yb - y});
    double dRod = signedDistancePolygon(roundedRodBoundaryLoop(g, 16), x, y);
    return std::max(dRect, -dRod);
}

Mesh makeBlastRodMesh(const BlastRodGeom& g, double h, int maxIter, bool verbose,
                      const BlastRodMeshOptions& options) {
    SolidBodyMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.hFar = h;
    spec.hNearFactor = options.hNearFactor;
    spec.gradeRadius = options.gradeRadius;
    spec.seedOffsetX = options.seedOffsetX;
    spec.seedOffsetY = options.seedOffsetY;
    spec.randomSeed = options.randomSeed;
    spec.maxIter = maxIter;
    spec.verbose = verbose;
    spec.solidDistance = [g](double x, double y) {
        return signedDistancePolygon(roundedRodBoundaryLoop(g, 16), x, y);
    };
    spec.fixedSegments = blastRodDomainSegments(g);
    addRoundedRodFluidBoundarySegments(spec, g);
    return makeSolidBodyFittedMesh(spec);
}

Mesh makeCurrentSolidBlastRodMesh(const BlastRodGeom& g, const ElasticSolid2D& solid,
                                  double h, int maxIter, bool verbose,
                                  const BlastRodMeshOptions& options) {
    std::vector<Vector2d> solidPoly = solidBoundaryLoop(solid);
    SolidBodyMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.hFar = h;
    spec.hNearFactor = options.hNearFactor;
    spec.gradeRadius = options.gradeRadius;
    spec.seedOffsetX = options.seedOffsetX;
    spec.seedOffsetY = options.seedOffsetY;
    spec.randomSeed = options.randomSeed;
    spec.maxIter = maxIter;
    spec.verbose = verbose;
    spec.solidDistance = [solidPoly](double x, double y) {
        return signedDistancePolygon(solidPoly, x, y);
    };
    spec.fixedSegments = blastRodDomainSegments(g);
    spec.solid = &solid;
    spec.addMovingSolidBoundary = true;
    return makeSolidBodyFittedMesh(spec);
}

int blastRodBoundaryTag(double x, double y, double t, const BlastRodMap& map) {
    const BlastRodGeom& g = map.geom;
    double dLeft = std::abs(x - g.xa);
    double dRight = std::abs(x - g.xb);
    double dBottom = std::abs(y - g.ya);
    double dTop = std::abs(y - g.yb);

    Vector2d p(x, y);
    double dRod = 1e9;
    const int ns = 48;
    for (int i = 0; i <= ns; ++i) {
        double yy = g.rodBaseY + g.rodL * static_cast<double>(i) / ns;
        dRod = std::min(dRod, (p - map.beamPoint(yy, -0.5 * g.rodW, t)).norm());
        dRod = std::min(dRod, (p - map.beamPoint(yy, 0.5 * g.rodW, t)).norm());
    }
    for (int i = 0; i <= 8; ++i) {
        double r = -0.5 * g.rodW + g.rodW * static_cast<double>(i) / 8.0;
        dRod = std::min(dRod, (p - map.beamPoint(g.rodTipY(), r, t)).norm());
    }

    double best = dLeft;
    int tag = TAG_EXACT;
    if (dRight < best) {
        best = dRight;
        tag = TAG_OUTFLOW;
    }
    if (dBottom < best) {
        best = dBottom;
        tag = TAG_SLIP_WALL;
    }
    if (dTop < best) {
        best = dTop;
        tag = TAG_SLIP_WALL;
    }
    if (dRod < best) tag = TAG_MOVING_WALL;
    return tag;
}

int blastRodBoundaryTag(double x, double y, double t, const SolidALEMap& map,
                        const BlastRodGeom& g) {
    double dLeft = std::abs(x - g.xa);
    double dRight = std::abs(x - g.xb);
    double dBottom = std::abs(y - g.ya);
    double dTop = std::abs(y - g.yb);
    double dRod = map.distanceToBoundary(x, y, t);

    double best = dLeft;
    int tag = TAG_EXACT;
    if (dRight < best) {
        best = dRight;
        tag = TAG_OUTFLOW;
    }
    if (dBottom < best) {
        best = dBottom;
        tag = TAG_SLIP_WALL;
    }
    if (dTop < best) {
        best = dTop;
        tag = TAG_SLIP_WALL;
    }
    if (dRod < best) tag = TAG_MOVING_WALL;
    return tag;
}

} // namespace euler_ale
