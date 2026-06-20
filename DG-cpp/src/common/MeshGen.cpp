#include "MeshGen.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef DG_USE_CGAL
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#endif

namespace {

// 2-D orientation: > 0 if (a,b,c) is counter-clockwise.
inline double orient2d(const Vector2d& a, const Vector2d& b, const Vector2d& c) {
    return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
}

// In-circle predicate for a CCW triangle (a,b,c): > 0 if d is strictly inside
// the circumcircle.  Returns the determinant value (caller compares with a tol).
inline double inCircleDet(const Vector2d& a, const Vector2d& b,
                          const Vector2d& c, const Vector2d& d) {
    double ax = a.x() - d.x(), ay = a.y() - d.y();
    double bx = b.x() - d.x(), by = b.y() - d.y();
    double cx = c.x() - d.x(), cy = c.y() - d.y();
    double a2 = ax * ax + ay * ay;
    double b2 = bx * bx + by * by;
    double c2 = cx * cx + cy * cy;
    return a2 * (bx * cy - cx * by) - b2 * (ax * cy - cx * ay) + c2 * (ax * by - bx * ay);
}

inline uint64_t edgeKey(int a, int b) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

inline uint64_t undirectedEdgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return edgeKey(a, b);
}

inline int edgeKeyA(uint64_t key) {
    return static_cast<int>(key >> 32);
}

inline int edgeKeyB(uint64_t key) {
    return static_cast<int>(key & 0xffffffffu);
}

inline uint64_t pointBucketKey(long long ix, long long iy) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(ix)) << 32) ^
           static_cast<uint32_t>(iy);
}

uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

double signedUnitFromHash(uint64_t x) {
    const double inv = 1.0 / static_cast<double>(uint64_t{1} << 53);
    return 2.0 * static_cast<double>(x >> 11) * inv - 1.0;
}

uint64_t triangleKey(int a, int b, int c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 42) ^
           (static_cast<uint64_t>(b) << 21) ^
           static_cast<uint64_t>(c);
}

int compactNearDuplicatePoints(std::vector<Vector2d>& points, int fixedCount, double tol) {
    if (points.empty() || !(tol > 0.0)) return fixedCount;
    const double inv = 1.0 / tol;
    const double tol2 = tol * tol;
    std::vector<Vector2d> kept;
    kept.reserve(points.size());
    std::unordered_map<uint64_t, std::vector<int>> buckets;
    buckets.reserve(points.size() * 2);

    auto cell = [&](const Vector2d& p) {
        return std::pair<long long, long long>{
            static_cast<long long>(std::floor(p.x() * inv)),
            static_cast<long long>(std::floor(p.y() * inv))
        };
    };

    auto addKept = [&](const Vector2d& p) {
        int idx = static_cast<int>(kept.size());
        kept.push_back(p);
        auto [ix, iy] = cell(p);
        buckets[pointBucketKey(ix, iy)].push_back(idx);
    };

    int newFixedCount = 0;
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        const Vector2d& p = points[i];
        auto [ix, iy] = cell(p);
        bool duplicate = false;
        for (long long dx = -1; dx <= 1 && !duplicate; ++dx) {
            for (long long dy = -1; dy <= 1 && !duplicate; ++dy) {
                auto it = buckets.find(pointBucketKey(ix + dx, iy + dy));
                if (it == buckets.end()) continue;
                for (int j : it->second) {
                    if ((p - kept[j]).squaredNorm() <= tol2) {
                        duplicate = true;
                        break;
                    }
                }
            }
        }
        if (duplicate) continue;
        if (i < fixedCount) ++newFixedCount;
        addKept(p);
    }

    points.swap(kept);
    return newFixedCount;
}

double wallSecondsSince(const std::chrono::steady_clock::time_point& start) {
    using seconds = std::chrono::duration<double>;
    return std::chrono::duration_cast<seconds>(std::chrono::steady_clock::now() - start).count();
}

struct DelaunayOptions {
    int maxLiveTriangles = 0;
    double jitter = 0.0;
    uint64_t jitterSeed = 0;
};

struct DelaunayResult {
    bool ok = true;
    std::string reason;
};

struct EdgeAccum {
    int a = -1;
    int b = -1;
    int count = 0;
};

#ifdef DG_USE_CGAL
DelaunayResult cgalDelaunayTriangulateChecked(const std::vector<Vector2d>& ptsIn,
                                              std::vector<Vector3i>& tris,
                                              const DelaunayOptions& opt) {
    tris.clear();
    int N = static_cast<int>(ptsIn.size());
    if (N < 3) return {};

    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using VertexBase = CGAL::Triangulation_vertex_base_with_info_2<int, Kernel>;
    using FaceBase = CGAL::Triangulation_face_base_2<Kernel>;
    using Tds = CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>;
    using Delaunay = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
    using Point = Kernel::Point_2;

    double xmin = ptsIn[0].x(), xmax = xmin, ymin = ptsIn[0].y(), ymax = ymin;
    for (const auto& p : ptsIn) {
        xmin = std::min(xmin, p.x()); xmax = std::max(xmax, p.x());
        ymin = std::min(ymin, p.y()); ymax = std::max(ymax, p.y());
    }
    double dmax = std::max(xmax - xmin, ymax - ymin);
    if (dmax <= 0.0) dmax = 1.0;
    double areaTol = 1e-28 * dmax * dmax;

    std::vector<std::pair<Point, int>> input;
    input.reserve(N);
    for (int i = 0; i < N; ++i) {
        double x = ptsIn[i].x();
        double y = ptsIn[i].y();
        if (opt.jitter > 0.0) {
            uint64_t h = splitmix64(opt.jitterSeed ^ static_cast<uint64_t>(i));
            x += opt.jitter * signedUnitFromHash(h);
            y += opt.jitter * signedUnitFromHash(splitmix64(h));
        }
        input.emplace_back(Point(x, y), i);
    }

    try {
        Delaunay dt;
        dt.insert(input.begin(), input.end());
        tris.reserve(static_cast<size_t>(dt.number_of_faces()));
        for (auto fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
            int a = fit->vertex(0)->info();
            int b = fit->vertex(1)->info();
            int c = fit->vertex(2)->info();
            if (a == b || b == c || c == a) continue;
            if (std::abs(orient2d(ptsIn[a], ptsIn[b], ptsIn[c])) <= areaTol) continue;
            if (orient2d(ptsIn[a], ptsIn[b], ptsIn[c]) < 0.0) std::swap(b, c);
            tris.emplace_back(a, b, c);
        }
    } catch (const std::exception& e) {
        tris.clear();
        return {false, std::string("CGAL Delaunay exception: ") + e.what()};
    }

    if (opt.maxLiveTriangles > 0 && static_cast<int>(tris.size()) > opt.maxLiveTriangles) {
        std::ostringstream oss;
        oss << "CGAL triangle count guard: raw=" << tris.size()
            << " limit=" << opt.maxLiveTriangles;
        tris.clear();
        return {false, oss.str()};
    }
    return {};
}
#endif

DelaunayResult delaunayTriangulateChecked(const std::vector<Vector2d>& ptsIn,
                                          std::vector<Vector3i>& tris,
                                          const DelaunayOptions& opt) {
    tris.clear();
    int N = static_cast<int>(ptsIn.size());
    if (N < 3) return {};

#ifdef DG_USE_CGAL
    DelaunayResult cgalResult = cgalDelaunayTriangulateChecked(ptsIn, tris, opt);
    if (cgalResult.ok) return cgalResult;
#endif

    // Bounding box -> super-triangle that comfortably contains all points.
    double xmin = ptsIn[0].x(), xmax = xmin, ymin = ptsIn[0].y(), ymax = ymin;
    for (const auto& p : ptsIn) {
        xmin = std::min(xmin, p.x()); xmax = std::max(xmax, p.x());
        ymin = std::min(ymin, p.y()); ymax = std::max(ymax, p.y());
    }
    double dmax = std::max(xmax - xmin, ymax - ymin);
    if (dmax <= 0) dmax = 1.0;
    double midx = 0.5 * (xmin + xmax), midy = 0.5 * (ymin + ymax);

    std::vector<Vector2d> pts = ptsIn;
    if (opt.jitter > 0.0) {
        for (int i = 0; i < N; ++i) {
            uint64_t h = splitmix64(opt.jitterSeed ^ static_cast<uint64_t>(i));
            pts[i].x() += opt.jitter * signedUnitFromHash(h);
            pts[i].y() += opt.jitter * signedUnitFromHash(splitmix64(h));
        }
    }
    pts.emplace_back(midx - 20.0 * dmax, midy - dmax);
    pts.emplace_back(midx,               midy + 20.0 * dmax);
    pts.emplace_back(midx + 20.0 * dmax, midy - dmax);
    int s0 = N, s1 = N + 1, s2 = N + 2;

    // Scale-aware tolerance for the in-circle test.  The deterministic jitter
    // handles exact cocircular ties; this tolerance is kept small so we do not
    // accidentally leave holes in the Bowyer-Watson cavity.
    const double tol = 1e-14 * dmax * dmax * dmax * dmax;
    const double areaTol = 1e-28 * dmax * dmax;

    auto pushCCW = [&](std::vector<Vector3i>& out, int a, int b, int c) {
        if (std::abs(orient2d(pts[a], pts[b], pts[c])) <= areaTol) return;
        if (orient2d(pts[a], pts[b], pts[c]) < 0) std::swap(b, c);
        out.emplace_back(a, b, c);
    };

    auto fail = [&](const std::string& reason) {
        tris.clear();
        return DelaunayResult{false, reason};
    };

    auto deduplicateTriangles = [&](std::vector<Vector3i>& T) {
        std::vector<Vector3i> unique;
        unique.reserve(T.size());
        std::unordered_set<uint64_t> seen;
        seen.reserve(T.size() * 2 + 1);
        for (auto t : T) {
            if (std::abs(orient2d(pts[t[0]], pts[t[1]], pts[t[2]])) <= areaTol) continue;
            if (orient2d(pts[t[0]], pts[t[1]], pts[t[2]]) < 0) std::swap(t[1], t[2]);
            uint64_t key = triangleKey(t[0], t[1], t[2]);
            if (!seen.insert(key).second) continue;
            unique.push_back(t);
        }
        T.swap(unique);
    };

    std::vector<Vector3i> T;
    pushCCW(T, s0, s1, s2);

    for (int ip = 0; ip < N; ++ip) {
        const Vector2d& p = pts[ip];

        // Edges of all "bad" triangles (those whose circumcircle holds p).
        // A valid Bowyer-Watson cavity boundary consists of edges used exactly
        // once by the bad-triangle set.  Counting undirected edges is more
        // conservative than checking only for a reverse directed edge: if a
        // prior degeneracy creates duplicate/overlapping triangles, non-manifold
        // edges are rejected instead of being promoted to new boundary edges.
        std::unordered_map<uint64_t, EdgeAccum> edgeCount;
        edgeCount.reserve(96);
        std::vector<Vector3i> good;
        good.reserve(T.size());
        bool nonManifoldCavity = false;

        auto addBadEdge = [&](int a, int b) {
            EdgeAccum& e = edgeCount[undirectedEdgeKey(a, b)];
            if (e.count == 0) {
                e.a = a;
                e.b = b;
            }
            ++e.count;
            if (e.count > 2) nonManifoldCavity = true;
        };

        for (const auto& t : T) {
            double det = inCircleDet(pts[t[0]], pts[t[1]], pts[t[2]], p);
            if (det > tol) {
                addBadEdge(t[0], t[1]);
                addBadEdge(t[1], t[2]);
                addBadEdge(t[2], t[0]);
            } else {
                good.push_back(t);
            }
        }
        if (nonManifoldCavity) {
            std::ostringstream oss;
            oss << "non-manifold cavity at point " << ip;
            return fail(oss.str());
        }

        // Cavity boundary = undirected edges used by exactly one bad triangle.
        good.reserve(good.size() + edgeCount.size());
        for (const auto& kv : edgeCount) {
            const EdgeAccum& e = kv.second;
            if (e.count == 1) pushCCW(good, e.a, e.b, ip);
        }
        T.swap(good);
        if (opt.maxLiveTriangles > 0 && static_cast<int>(T.size()) > opt.maxLiveTriangles) {
            deduplicateTriangles(T);
            if (static_cast<int>(T.size()) > opt.maxLiveTriangles) {
                std::ostringstream oss;
                oss << "triangle growth guard at point " << ip
                    << ": live=" << T.size()
                    << " limit=" << opt.maxLiveTriangles;
                return fail(oss.str());
            }
        }
    }

    // Drop triangles that still reference a super-triangle vertex.
    deduplicateTriangles(T);
    if (opt.maxLiveTriangles > 0 && static_cast<int>(T.size()) > opt.maxLiveTriangles) {
        std::ostringstream oss;
        oss << "triangle growth guard after dedup: live=" << T.size()
            << " limit=" << opt.maxLiveTriangles;
        return fail(oss.str());
    }
    tris.reserve(T.size());
    for (const auto& t : T) {
        if (t[0] >= N || t[1] >= N || t[2] >= N) continue;
        tris.push_back(t);
    }
    return {};
}

} // namespace

double CylinderGeom::sdist(double x, double y) const {
    double dRect = -std::min(std::min(y - ya, yb - y), std::min(x - xa, xb - x));
    double dCirc = std::hypot(x - cx, y - cy) - r;
    return std::max(dRect, -dCirc);
}

double BowlGeom::sdist(double x, double y) const {
    return std::hypot(x - cx, y - cy) - R;
}

void delaunayTriangulate(const std::vector<Vector2d>& ptsIn,
                         std::vector<Vector3i>& tris) {
    DelaunayOptions opt;
    opt.maxLiveTriangles = std::max(128, 12 * static_cast<int>(ptsIn.size()));
    DelaunayResult result = delaunayTriangulateChecked(ptsIn, tris, opt);
    if (!result.ok) tris.clear();
}

int generateDistanceMesh(Mesh& mesh, const DistanceMeshSpec& spec,
                         int maxIter, bool verbose) {
    if (!spec.signedDistance) throw std::runtime_error("generateDistanceMesh: missing signedDistance");
    if (!spec.targetSize) throw std::runtime_error("generateDistanceMesh: missing targetSize");
    if (!(spec.h0 > 0.0)) throw std::runtime_error("generateDistanceMesh: h0 must be positive");

    const double h0 = spec.h0;
    const double seedH = (spec.seedH > 0.0) ? spec.seedH : h0;
    const double hTol = std::min(h0, seedH);
    const double geps = 1e-3 * hTol;
    const double deps = 1e-7 * hTol;
    const double Fscale = 1.2;
    const double deltat = 0.2;
    const double ttol = 0.10;
    const double dptol = 5e-3;

    std::vector<Vector2d> p = spec.fixedPoints;
    int nFix = static_cast<int>(p.size());

    std::mt19937 rng(spec.randomSeed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double dy = seedH * std::sqrt(3.0) / 2.0;
    const double r0max = 1.0 / (seedH * seedH);
    double yOffset = std::fmod(std::max(0.0, spec.seedOffsetY), dy);
    double xOffset = std::fmod(std::max(0.0, spec.seedOffsetX), seedH);
    int row = 0;
    for (double y = spec.ya + yOffset; y <= spec.yb + 1e-9; y += dy, ++row) {
        double xoff = (row % 2) ? 0.5 * seedH : 0.0;
        for (double x = spec.xa + xOffset + xoff; x <= spec.xb + 1e-9; x += seedH) {
            if (spec.signedDistance(x, y) > -geps) continue;
            double s = std::max(spec.targetSize(x, y), 1e-12);
            double r0 = 1.0 / (s * s);
            if (U(rng) > std::min(1.0, r0 / r0max)) continue;
            bool tooClose = false;
            for (int j = 0; j < nFix; ++j) {
                if ((p[j] - Vector2d(x, y)).norm() < 0.6 * spec.targetSize(x, y)) {
                    tooClose = true;
                    break;
                }
            }
            if (!tooClose) p.emplace_back(x, y);
        }
    }

    const double duplicateTol = 1e-4 * hTol;
    int nBeforeCompact = static_cast<int>(p.size());
    nFix = compactNearDuplicatePoints(p, nFix, duplicateTol);
    int Np = static_cast<int>(p.size());
    if (verbose) {
        std::cout << "  distance mesh seed points=" << Np;
        if (Np != nBeforeCompact) std::cout << " (dedup " << nBeforeCompact - Np << ")";
        std::cout << " fixed=" << nFix << "\n" << std::flush;
    }
    auto meshStart = std::chrono::steady_clock::now();
    std::vector<Vector2d> pold(Np, Vector2d(1e10, 1e10));
    std::vector<Vector3i> tris;

    auto bigMove = [&]() {
        if (static_cast<int>(pold.size()) != Np) return 1e300;
        double m = 0.0;
        for (int i = 0; i < Np; ++i) m = std::max(m, (p[i] - pold[i]).norm());
        return m;
    };

    int retriangulations = 0;
    auto retriangulate = [&](bool logThis) {
        auto start = std::chrono::steady_clock::now();
        const std::vector<Vector2d> basePoints = p;
        const int baseFix = nFix;
        const int baseCount = static_cast<int>(basePoints.size());
        if (logThis) {
            std::cout << "  Delaunay start: points=" << baseCount
                      << " fixed=" << baseFix
                      << " call=" << (retriangulations + 1) << "\n" << std::flush;
        }

        std::vector<Vector3i> allTris;
        std::vector<Vector2d> acceptedPoints;
        int acceptedFix = baseFix;
        int acceptedDedup = 0;
        int acceptedAttempt = -1;
        double acceptedTol = duplicateTol;
        double acceptedJitter = 0.0;
        std::string lastReason;

        const int maxAttempts = 6;
        for (int attempt = 0; attempt < maxAttempts; ++attempt) {
            std::vector<Vector2d> trialPoints = basePoints;
            double tolThis = duplicateTol * std::pow(10.0, std::min(attempt, 3));
            int trialFix = compactNearDuplicatePoints(trialPoints, baseFix, tolThis);
            int trialCount = static_cast<int>(trialPoints.size());
            int maxLiveTris = std::max(256, 12 * (trialCount + 3));

            DelaunayOptions opt;
            opt.maxLiveTriangles = maxLiveTris;
            opt.jitter = (attempt == 0)
                             ? 0.0
                             : std::pow(10.0, attempt - 13) *
                                   std::max(spec.xb - spec.xa, spec.yb - spec.ya);
            opt.jitterSeed = 0x5eed1234ull + 7919ull * static_cast<uint64_t>(attempt);

            std::vector<Vector3i> trialTris;
            DelaunayResult result = delaunayTriangulateChecked(trialPoints, trialTris, opt);
            if (result.ok && static_cast<int>(trialTris.size()) <= maxLiveTris) {
                acceptedPoints = std::move(trialPoints);
                acceptedFix = trialFix;
                acceptedDedup = baseCount - trialCount;
                acceptedAttempt = attempt;
                acceptedTol = tolThis;
                acceptedJitter = opt.jitter;
                allTris = std::move(trialTris);
                break;
            }

            lastReason = result.ok ? "raw triangle count exceeded guard" : result.reason;
            if (logThis) {
                std::cout << "  Delaunay retry: attempt=" << attempt
                          << " points=" << trialCount
                          << " dedup=" << (baseCount - trialCount)
                          << " tol=" << tolThis
                          << " jitter=" << opt.jitter
                          << " reason=" << lastReason << "\n" << std::flush;
            }
        }

        if (acceptedAttempt < 0) {
            std::ostringstream oss;
            oss << "generateDistanceMesh: Delaunay failed after retries";
            if (!lastReason.empty()) oss << " (" << lastReason << ")";
            throw std::runtime_error(oss.str());
        }

        p = std::move(acceptedPoints);
        nFix = acceptedFix;
        Np = static_cast<int>(p.size());
        if (static_cast<int>(pold.size()) != Np) pold.assign(Np, Vector2d(1e10, 1e10));

        tris.clear();
        for (const auto& t : allTris) {
            Vector2d cmid = (p[t[0]] + p[t[1]] + p[t[2]]) / 3.0;
            if (spec.signedDistance(cmid.x(), cmid.y()) < -geps) tris.push_back(t);
        }
        ++retriangulations;
        if (logThis) {
            std::cout << "  Delaunay done: raw_tris=" << allTris.size()
                      << " kept=" << tris.size()
                      << " points=" << Np
                      << " fixed=" << nFix
                      << " dedup=" << acceptedDedup
                      << " attempt=" << acceptedAttempt
                      << " tol=" << acceptedTol
                      << " jitter=" << acceptedJitter
                      << " wall=" << wallSecondsSince(start) << "s\n" << std::flush;
        }
    };

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        if (bigMove() > ttol * hTol) {
            pold = p;
            retriangulate(verbose);
        }

        std::set<std::pair<int, int>> barSet;
        for (const auto& t : tris) {
            int a = t[0], b = t[1], c = t[2];
            barSet.insert({std::min(a, b), std::max(a, b)});
            barSet.insert({std::min(b, c), std::max(b, c)});
            barSet.insert({std::min(c, a), std::max(c, a)});
        }
        std::vector<std::pair<int, int>> bars(barSet.begin(), barSet.end());

        int nB = static_cast<int>(bars.size());
        std::vector<double> L(nB), L0(nB);
        double sumL2 = 0.0, sumH2 = 0.0;
        for (int i = 0; i < nB; ++i) {
            Vector2d d = p[bars[i].second] - p[bars[i].first];
            L[i] = std::max(d.norm(), 1e-14);
            Vector2d mid = 0.5 * (p[bars[i].first] + p[bars[i].second]);
            double hb = std::max(spec.targetSize(mid.x(), mid.y()), 1e-12);
            L0[i] = hb;
            sumL2 += L[i] * L[i];
            sumH2 += hb * hb;
        }
        double scale = Fscale * std::sqrt(sumL2 / std::max(sumH2, 1e-300));

        std::vector<Vector2d> F(Np, Vector2d::Zero());
        for (int i = 0; i < nB; ++i) {
            Vector2d barvec = p[bars[i].second] - p[bars[i].first];
            double l0 = L0[i] * scale;
            double Fbar = std::max(l0 - L[i], 0.0);
            Vector2d fv = (Fbar / L[i]) * barvec;
            F[bars[i].first] -= fv;
            F[bars[i].second] += fv;
        }

        double maxInteriorMove = 0.0;
        for (int i = nFix; i < Np; ++i) {
            Vector2d step = deltat * F[i];
            p[i] += step;
            double d = spec.signedDistance(p[i].x(), p[i].y());
            if (d > 0.0) {
                if (spec.projectInside) {
                    p[i] = spec.projectInside(p[i]);
                } else {
                    double dgx = (spec.signedDistance(p[i].x() + deps, p[i].y()) - d) / deps;
                    double dgy = (spec.signedDistance(p[i].x(), p[i].y() + deps) - d) / deps;
                    p[i] -= Vector2d(d * dgx, d * dgy);
                }
            } else {
                maxInteriorMove = std::max(maxInteriorMove, step.norm());
            }
        }

        if (verbose && (iter < 3 || iter % 10 == 0)) {
            std::cout << "  DistMesh iter=" << iter
                      << " points=" << Np
                      << " tris=" << tris.size()
                      << " bars=" << nB
                      << " max_move=" << maxInteriorMove
                      << " retriangulations=" << retriangulations
                      << " wall=" << wallSecondsSince(meshStart) << "s\n" << std::flush;
        }
        if (maxInteriorMove < dptol * hTol && iter > 5) { ++iter; break; }
    }

    retriangulate(verbose);

    std::vector<int> remap(p.size(), -1);
    std::vector<Vector2d> usedNodes;
    for (auto& t : tris) {
        for (int k = 0; k < 3; ++k) {
            if (remap[t[k]] < 0) {
                remap[t[k]] = static_cast<int>(usedNodes.size());
                usedNodes.push_back(p[t[k]]);
            }
        }
    }

    mesh.node.resize(static_cast<int>(usedNodes.size()), 2);
    for (int i = 0; i < static_cast<int>(usedNodes.size()); ++i) {
        mesh.node(i, 0) = usedNodes[i].x();
        mesh.node(i, 1) = usedNodes[i].y();
    }
    mesh.elem.resize(static_cast<int>(tris.size()), 3);
    for (int i = 0; i < static_cast<int>(tris.size()); ++i) {
        int a = remap[tris[i][0]], b = remap[tris[i][1]], c = remap[tris[i][2]];
        Vector2d pa = usedNodes[a], pb = usedNodes[b], pc = usedNodes[c];
        if (orient2d(pa, pb, pc) < 0) std::swap(b, c);
        mesh.elem(i, 0) = a;
        mesh.elem(i, 1) = b;
        mesh.elem(i, 2) = c;
    }

    if (verbose) {
        double minA, meanA, minAr, maxAr;
        meshQuality(mesh, minA, meanA, minAr, maxAr);
        std::cout << "  distance mesh: " << mesh.node.rows() << " nodes, " << mesh.elem.rows()
                  << " triangles (" << iter << " DistMesh iters)\n"
                  << "  quality: min angle " << minA << " deg, mean angle " << meanA
                  << " deg; area in [" << minAr << ", " << maxAr << "]\n";
    }
    return iter;
}

// ---------------------------------------------------------------------------
// DistMesh size function: ~ h on the cylinder, growing to h*farRatio away.
// ---------------------------------------------------------------------------
namespace {
struct SizeFn {
    CylinderGeom g;
    double h, farRatio, grade;
    double operator()(double x, double y) const {
        double dOut = std::max(0.0, std::hypot(x - g.cx, y - g.cy) - g.r); // dist outside cylinder
        return std::min(h * farRatio, h * (1.0 + grade / h * dOut));
    }
};
} // namespace

int generateCylinderMesh(Mesh& mesh, const CylinderGeom& g, double h,
                         double sizeFarRatio, double gradeRate,
                         int maxIter, bool verbose) {
    SizeFn fh{g, h, sizeFarRatio, gradeRate};
    DistanceMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.h0 = h;
    spec.signedDistance = [g](double x, double y) { return g.sdist(x, y); };
    spec.targetSize = [fh](double x, double y) { return fh(x, y); };
    spec.fixedPoints.emplace_back(g.xa, g.ya);
    spec.fixedPoints.emplace_back(g.xb, g.ya);
    spec.fixedPoints.emplace_back(g.xb, g.yb);
    spec.fixedPoints.emplace_back(g.xa, g.yb);
    int nCyl = std::max(16, static_cast<int>(std::round(2.0 * M_PI * g.r / h)));
    for (int k = 0; k < nCyl; ++k) {
        double th = 2.0 * M_PI * k / nCyl;
        spec.fixedPoints.emplace_back(g.cx + g.r * std::cos(th), g.cy + g.r * std::sin(th));
    }
    return generateDistanceMesh(mesh, spec, maxIter, verbose);
}

VectorXi classifyEdges(const Mesh& mesh, const MatrixXi& edge,
                       const MatrixXi& edge2side, const CylinderGeom& g) {
    int NE = edge.rows();
    VectorXi tag = VectorXi::Constant(NE, BD_INTERIOR);
    for (int e = 0; e < NE; ++e) {
        if (edge2side(e, 0) != -1 && edge2side(e, 1) != -1) continue; // interior
        Vector2d p1 = mesh.node.row(edge(e, 0));
        Vector2d p2 = mesh.node.row(edge(e, 1));
        Vector2d m = 0.5 * (p1 + p2);
        double dL = std::abs(m.x() - g.xa);
        double dR = std::abs(m.x() - g.xb);
        double dB = std::abs(m.y() - g.ya);
        double dT = std::abs(m.y() - g.yb);
        double dC = std::abs(std::hypot(m.x() - g.cx, m.y() - g.cy) - g.r);
        double best = dL; int which = BD_INFLOW;
        if (dR < best) { best = dR; which = BD_OUTFLOW; }
        if (dB < best) { best = dB; which = BD_WALL; }
        if (dT < best) { best = dT; which = BD_WALL; }
        if (dC < best) { best = dC; which = BD_CYL; }
        tag(e) = which;
    }
    return tag;
}

// ---------------------------------------------------------------------------
// Closed circular bowl: DistMesh on a disk SDF (no interior hole).
// ---------------------------------------------------------------------------
namespace {
struct BowlSizeFn {
    BowlGeom g;
    double h, farRatio, grade, bandLo, bandHi;
    double operator()(double x, double y) const {
        double r = std::hypot(x - g.cx, y - g.cy);
        // distance from the fine annular band [bandLo, bandHi]
        double d = (r < bandLo) ? (bandLo - r) : (r > bandHi ? r - bandHi : 0.0);
        return std::min(h * farRatio, h * (1.0 + grade / h * d));
    }
};
} // namespace

int generateBowlMesh(Mesh& mesh, const BowlGeom& g, double h,
                     double farRatio, double gradeRate, double bandLo, double bandHi,
                     int maxIter, bool verbose) {
    BowlSizeFn fh{g, h, farRatio, gradeRate, bandLo, bandHi};
    DistanceMeshSpec spec;
    spec.xa = g.cx - g.R;
    spec.xb = g.cx + g.R;
    spec.ya = g.cy - g.R;
    spec.yb = g.cy + g.R;
    spec.h0 = h;
    spec.signedDistance = [g](double x, double y) { return g.sdist(x, y); };
    spec.targetSize = [fh](double x, double y) { return fh(x, y); };
    spec.projectInside = [g](const Vector2d& p) {
        double r = std::hypot(p.x() - g.cx, p.y() - g.cy);
        if (r <= g.R) return p;
        double s = g.R / std::max(r, 1e-14);
        return Vector2d(g.cx + (p.x() - g.cx) * s,
                        g.cy + (p.y() - g.cy) * s);
    };

    double hRim = fh(g.cx + g.R, g.cy);
    int nRim = std::max(24, static_cast<int>(std::round(2.0 * M_PI * g.R / hRim)));
    for (int k = 0; k < nRim; ++k) {
        double th = 2.0 * M_PI * k / nRim;
        spec.fixedPoints.emplace_back(g.cx + g.R * std::cos(th), g.cy + g.R * std::sin(th));
    }
    return generateDistanceMesh(mesh, spec, maxIter, verbose);
}

VectorXi classifyBowlEdges(const Mesh& mesh, const MatrixXi& edge,
                           const MatrixXi& edge2side, const BowlGeom& g) {
    (void)g;
    int NE = edge.rows();
    VectorXi tag = VectorXi::Constant(NE, BD_BOWL_INTERIOR);
    for (int e = 0; e < NE; ++e)
        if (edge2side(e, 0) == -1 || edge2side(e, 1) == -1)
            tag(e) = BD_BOWL_WALL;
    return tag;
}

void meshQuality(const Mesh& mesh, double& minAngleDeg, double& meanAngleDeg,
                 double& minArea, double& maxArea) {
    int NT = mesh.elem.rows();
    minAngleDeg = 180.0; meanAngleDeg = 0.0;
    minArea = 1e300; maxArea = 0.0;
    long cnt = 0;
    for (int t = 0; t < NT; ++t) {
        Vector2d a = mesh.node.row(mesh.elem(t, 0));
        Vector2d b = mesh.node.row(mesh.elem(t, 1));
        Vector2d c = mesh.node.row(mesh.elem(t, 2));
        double area = 0.5 * std::abs(orient2d(a, b, c));
        minArea = std::min(minArea, area);
        maxArea = std::max(maxArea, area);
        Vector2d v[3] = {a, b, c};
        for (int k = 0; k < 3; ++k) {
            Vector2d e1 = v[(k + 1) % 3] - v[k];
            Vector2d e2 = v[(k + 2) % 3] - v[k];
            double ang = std::acos(std::max(-1.0, std::min(1.0,
                            e1.dot(e2) / (e1.norm() * e2.norm())))) * 180.0 / M_PI;
            minAngleDeg = std::min(minAngleDeg, ang);
            meanAngleDeg += ang;
            ++cnt;
        }
    }
    if (cnt) meanAngleDeg /= cnt;
}
