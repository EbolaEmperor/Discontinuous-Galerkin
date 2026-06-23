#ifndef EULER_GEO_UTILS_H
#define EULER_GEO_UTILS_H

#include "Mesh.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

namespace euler {

using namespace Eigen;

template <class Worker>
void parallelRanges(int n, const Worker& worker) {
    unsigned hw = std::thread::hardware_concurrency();
    unsigned nt = std::max(1u, std::min<unsigned>(hw ? hw : 1u, (unsigned)std::max(1, n)));
    if (nt <= 1 || n < 64) { worker(0, n); return; }
    int chunk = (n + (int)nt - 1) / (int)nt;
    std::vector<std::thread> ths;
    for (unsigned t = 0; t < nt; ++t) {
        int lo = (int)t * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        ths.emplace_back([lo, hi, &worker] { worker(lo, hi); });
    }
    for (auto& th : ths) th.join();
}

inline double smooth01(double s) {
    s = std::clamp(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
}

inline double orient2d(const Vector2d& a, const Vector2d& b, const Vector2d& c) {
    return (b.x() - a.x()) * (c.y() - a.y()) -
           (b.y() - a.y()) * (c.x() - a.x());
}

inline double signedBox(double x, double y, double x0, double x1, double y0, double y1) {
    double dx = std::max({x0 - x, 0.0, x - x1});
    double dy = std::max({y0 - y, 0.0, y - y1});
    double outside = std::hypot(dx, dy);
    if (dx > 0.0 || dy > 0.0) return outside;
    return -std::min({x - x0, x1 - x, y - y0, y1 - y});
}

inline double segmentProjection(const Vector2d& p, const Vector2d& a, const Vector2d& b,
                                double& dist2) {
    Vector2d ab = b - a;
    double den = ab.squaredNorm();
    double s = 0.0;
    if (den > 1e-30) s = std::clamp((p - a).dot(ab) / den, 0.0, 1.0);
    Vector2d q = a + s * ab;
    dist2 = (p - q).squaredNorm();
    return s;
}

inline double triArea(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0));
    Vector2d p1 = m.node.row(m.elem(t, 1));
    Vector2d p2 = m.node.row(m.elem(t, 2));
    return 0.5 * std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
                          (p2.x() - p0.x()) * (p1.y() - p0.y()));
}

inline double longestEdge(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0));
    Vector2d p1 = m.node.row(m.elem(t, 1));
    Vector2d p2 = m.node.row(m.elem(t, 2));
    return std::max({(p1 - p0).norm(), (p2 - p1).norm(), (p0 - p2).norm()});
}

inline double hCFL(const Mesh& m, int t) {
    return 2.0 * triArea(m, t) / std::max(longestEdge(m, t), 1e-300);
}

} // namespace euler

#endif
