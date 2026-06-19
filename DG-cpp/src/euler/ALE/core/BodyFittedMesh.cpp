#include "BodyFittedMesh.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace euler_ale {
namespace {

void addFixedPointUnique(DistanceMeshSpec& spec, const Vector2d& p, double tol) {
    for (const auto& q : spec.fixedPoints) {
        if ((p - q).norm() < tol) return;
    }
    spec.fixedPoints.push_back(p);
}

double smoothRamp(double s) {
    s = std::clamp(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
}

} // namespace

double signedBox(double x, double y, double x0, double x1, double y0, double y1) {
    double dx = std::max({x0 - x, 0.0, x - x1});
    double dy = std::max({y0 - y, 0.0, y - y1});
    double outside = std::hypot(dx, dy);
    if (dx > 0.0 || dy > 0.0) return outside;
    return -std::min({x - x0, x1 - x, y - y0, y1 - y});
}

double segmentProjection(const Vector2d& p, const Vector2d& a, const Vector2d& b,
                         double& dist2) {
    Vector2d ab = b - a;
    double den = ab.squaredNorm();
    double s = 0.0;
    if (den > 1e-30) s = std::clamp((p - a).dot(ab) / den, 0.0, 1.0);
    Vector2d q = a + s * ab;
    dist2 = (p - q).squaredNorm();
    return s;
}

double signedDistancePolygon(const std::vector<Vector2d>& poly, double x, double y) {
    if (poly.size() < 3) return 1e300;
    Vector2d p(x, y);
    double d2min = 1e300;
    bool inside = false;
    for (int i = 0, j = static_cast<int>(poly.size()) - 1;
         i < static_cast<int>(poly.size()); j = i++) {
        double d2 = 0.0;
        segmentProjection(p, poly[j], poly[i], d2);
        d2min = std::min(d2min, d2);
        double denom = poly[j].y() - poly[i].y();
        bool cross = ((poly[i].y() > y) != (poly[j].y() > y)) &&
                     (x < (poly[j].x() - poly[i].x()) * (y - poly[i].y()) /
                              denom +
                          poly[i].x());
        if (cross) inside = !inside;
    }
    double d = std::sqrt(d2min);
    return inside ? -d : d;
}

std::vector<Vector2d> solidBoundaryLoop(const ElasticSolid2D& solid) {
    const MatrixXd& x = solid.currentNodes();
    std::vector<SolidBoundarySegment> bottom, right, top, left;
    for (const auto& s : solid.boundarySegments()) {
        if (s.side == SOLID_BOTTOM) bottom.push_back(s);
        else if (s.side == SOLID_RIGHT) right.push_back(s);
        else if (s.side == SOLID_TOP) top.push_back(s);
        else if (s.side == SOLID_LEFT) left.push_back(s);
    }
    auto midx = [&](const SolidBoundarySegment& s) { return 0.5 * (x(s.a, 0) + x(s.b, 0)); };
    auto midy = [&](const SolidBoundarySegment& s) { return 0.5 * (x(s.a, 1) + x(s.b, 1)); };
    std::sort(bottom.begin(), bottom.end(), [&](const auto& a, const auto& b) { return midx(a) < midx(b); });
    std::sort(right.begin(), right.end(), [&](const auto& a, const auto& b) { return midy(a) < midy(b); });
    std::sort(top.begin(), top.end(), [&](const auto& a, const auto& b) { return midx(a) > midx(b); });
    std::sort(left.begin(), left.end(), [&](const auto& a, const auto& b) { return midy(a) > midy(b); });

    std::vector<Vector2d> loop;
    auto appendPoint = [&](const Vector2d& p) {
        if (loop.empty() || (loop.back() - p).norm() > 1e-12) loop.push_back(p);
    };
    auto appendSeg = [&](const SolidBoundarySegment& s, bool forward) {
        Vector2d a = x.row(forward ? s.a : s.b).transpose();
        Vector2d b = x.row(forward ? s.b : s.a).transpose();
        appendPoint(a);
        appendPoint(b);
    };
    for (const auto& s : bottom) appendSeg(s, true);
    for (const auto& s : right) appendSeg(s, true);
    for (const auto& s : top) appendSeg(s, false);
    for (const auto& s : left) appendSeg(s, false);
    if (!loop.empty() && (loop.front() - loop.back()).norm() < 1e-12) loop.pop_back();
    return loop;
}

void addGradedFixedSegment(DistanceMeshSpec& spec, const Vector2d& a, const Vector2d& b,
                           const std::function<double(double, double)>& sizeFn,
                           double hMin) {
    double L = (b - a).norm();
    double tol = 0.15 * hMin;
    if (L < 1e-14) {
        addFixedPointUnique(spec, a, tol);
        return;
    }

    Vector2d dir = (b - a) / L;
    double s = 0.0;
    while (s < L) {
        Vector2d p = a + s * dir;
        addFixedPointUnique(spec, p, tol);
        double hLocal = std::max(hMin, sizeFn(p.x(), p.y()));
        s += hLocal;
        if (L - s < 0.45 * hLocal) break;
    }
    addFixedPointUnique(spec, b, tol);
}

Mesh makeSolidBodyFittedMesh(const SolidBodyMeshSpec& cfg) {
    double hNear = cfg.hNearFactor * cfg.hFar;
    auto solidDistance = cfg.solidDistance;
    auto sizeField = [solidDistance, hNear, hFar = cfg.hFar,
                      gradeRadius = cfg.gradeRadius](double x, double y) {
        double dSolid = std::max(0.0, solidDistance(x, y));
        double a = smoothRamp(dSolid / gradeRadius);
        return hNear + (hFar - hNear) * a;
    };

    DistanceMeshSpec spec;
    spec.xa = cfg.xa;
    spec.xb = cfg.xb;
    spec.ya = cfg.ya;
    spec.yb = cfg.yb;
    spec.h0 = cfg.hFar;
    spec.seedH = hNear;
    spec.signedDistance = [cfg, solidDistance](double x, double y) {
        double dRect = -std::min({x - cfg.xa, cfg.xb - x, y - cfg.ya, cfg.yb - y});
        double dSolid = solidDistance(x, y);
        return std::max(dRect, -dSolid);
    };
    spec.targetSize = sizeField;

    for (const auto& segment : cfg.fixedSegments) {
        addGradedFixedSegment(spec, segment.first, segment.second, sizeField, hNear);
    }
    if (cfg.addMovingSolidBoundary && cfg.solid) {
        const MatrixXd& sx = cfg.solid->currentNodes();
        for (const auto& seg : cfg.solid->movingBoundarySegments()) {
            Vector2d a = sx.row(seg.a).transpose();
            Vector2d b = sx.row(seg.b).transpose();
            addGradedFixedSegment(spec, a, b, sizeField, hNear);
        }
    }

    if (cfg.verbose) {
        std::cout << "  body-fitted mesh size field: h_near=" << hNear
                  << " h_far=" << cfg.hFar
                  << " grade_radius=" << cfg.gradeRadius << "\n";
    }

    Mesh mesh;
    generateDistanceMesh(mesh, spec, cfg.maxIter, cfg.verbose);
    return mesh;
}

} // namespace euler_ale
