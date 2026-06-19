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
        {Vector2d(g.xa, g.ya), Vector2d(g.rodLeft(), g.ya)},
        {Vector2d(g.rodRight(), g.ya), Vector2d(g.xb, g.ya)},
        {Vector2d(g.xb, g.ya), Vector2d(g.xb, g.yb)},
        {Vector2d(g.xb, g.yb), Vector2d(g.xa, g.yb)},
        {Vector2d(g.xa, g.yb), Vector2d(g.xa, g.ya)}
    };
}

} // namespace

double blastRodFluidSdf(const BlastRodGeom& g, double x, double y) {
    double dRect = -std::min({x - g.xa, g.xb - x, y - g.ya, g.yb - y});
    double dRod = signedBox(x, y, g.rodLeft(), g.rodRight(), g.rodBaseY, g.rodTipY());
    return std::max(dRect, -dRod);
}

Mesh makeBlastRodMesh(const BlastRodGeom& g, double h, int maxIter, bool verbose) {
    SolidBodyMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.hFar = h;
    spec.maxIter = maxIter;
    spec.verbose = verbose;
    spec.solidDistance = [g](double x, double y) {
        return signedBox(x, y, g.rodLeft(), g.rodRight(), g.rodBaseY, g.rodTipY());
    };
    spec.fixedSegments = blastRodDomainSegments(g);
    spec.fixedSegments.push_back({Vector2d(g.rodLeft(), g.rodBaseY),
                                  Vector2d(g.rodLeft(), g.rodTipY())});
    spec.fixedSegments.push_back({Vector2d(g.rodLeft(), g.rodTipY()),
                                  Vector2d(g.rodRight(), g.rodTipY())});
    spec.fixedSegments.push_back({Vector2d(g.rodRight(), g.rodTipY()),
                                  Vector2d(g.rodRight(), g.rodBaseY)});
    return makeSolidBodyFittedMesh(spec);
}

Mesh makeCurrentSolidBlastRodMesh(const BlastRodGeom& g, const ElasticSolid2D& solid,
                                  double h, int maxIter, bool verbose) {
    std::vector<Vector2d> solidPoly = solidBoundaryLoop(solid);
    SolidBodyMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.hFar = h;
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
