#include "BodyFittedMesh.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <sstream>
#include <unordered_map>

#ifdef DG_USE_CGAL
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/number_utils.h>
#endif

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

std::vector<Vector2d> sampleGradedSegment(const Vector2d& a, const Vector2d& b,
                                          const std::function<double(double, double)>& sizeFn,
                                          double hMin) {
    std::vector<Vector2d> pts;
    double L = (b - a).norm();
    if (L < 1e-14) {
        pts.push_back(a);
        return pts;
    }

    Vector2d dir = (b - a) / L;
    double s = 0.0;
    while (s < L) {
        Vector2d p = a + s * dir;
        pts.push_back(p);
        double hLocal = std::max(hMin, sizeFn(p.x(), p.y()));
        s += hLocal;
        if (L - s < 0.45 * hLocal) break;
    }
    pts.push_back(b);
    return pts;
}

void addGradedFixedClosedLoop(DistanceMeshSpec& spec, const std::vector<Vector2d>& loop,
                              const std::function<double(double, double)>& sizeFn,
                              double hMin) {
    if (loop.size() < 3) return;
    double tol = 0.15 * hMin;
    double total = 0.0;
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        total += (loop[(i + 1) % loop.size()] - loop[i]).norm();
    }
    if (!(total > 0.0)) return;

    double arclen = 0.0;
    double nextSample = 0.0;
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        const Vector2d& a = loop[i];
        const Vector2d& b = loop[(i + 1) % loop.size()];
        double L = (b - a).norm();
        if (L <= 1e-14) continue;
        while (nextSample <= arclen + L + 1e-12 && nextSample < total - 1e-12) {
            double s = std::clamp((nextSample - arclen) / L, 0.0, 1.0);
            Vector2d p = a + s * (b - a);
            addFixedPointUnique(spec, p, tol);
            double hLocal = std::max(hMin, sizeFn(p.x(), p.y()));
            nextSample += hLocal;
        }
        arclen += L;
    }

    const double sharpCos = std::cos(35.0 * M_PI / 180.0);
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        const Vector2d& prev = loop[(i + static_cast<int>(loop.size()) - 1) %
                                    loop.size()];
        const Vector2d& p = loop[i];
        const Vector2d& next = loop[(i + 1) % loop.size()];
        Vector2d u = p - prev;
        Vector2d v = next - p;
        double lu = u.norm();
        double lv = v.norm();
        if (lu <= 1e-14 || lv <= 1e-14) continue;
        double c = (u / lu).dot(v / lv);
        if (c < sharpCos) addFixedPointUnique(spec, p, tol);
    }
}

#ifdef DG_USE_CGAL
using CgalKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using CgalVertexBase = CGAL::Triangulation_vertex_base_2<CgalKernel>;
using CgalFaceBase = CGAL::Delaunay_mesh_face_base_2<CgalKernel>;
using CgalTds = CGAL::Triangulation_data_structure_2<CgalVertexBase, CgalFaceBase>;
using CgalCdt = CGAL::Constrained_Delaunay_triangulation_2<
    CgalKernel, CgalTds, CGAL::Exact_intersections_tag>;
using CgalPoint = CgalKernel::Point_2;

template <class CDT>
class LocalSizeCriteria : public CGAL::Delaunay_mesh_criteria_2<CDT> {
public:
    using Base = CGAL::Delaunay_mesh_criteria_2<CDT>;
    using GeomTraits = typename CDT::Geom_traits;
    using FaceHandle = typename CDT::Face_handle;

    struct Quality : public std::pair<double, double> {
        Quality() : std::pair<double, double>(1.0, 0.0) {}
        Quality(double sine, double size) : std::pair<double, double>(sine, size) {}

        double sine() const { return this->first; }
        double size() const { return this->second; }

        bool operator<(const Quality& q) const {
            if (size() > 1.0) {
                if (q.size() > 1.0) return size() > q.size();
                return true;
            }
            if (q.size() > 1.0) return false;
            return sine() < q.sine();
        }
    };

    LocalSizeCriteria(double minAngleDeg, double maxEdgeFactor,
                      std::function<double(double, double)> targetSize,
                      const GeomTraits& traits = GeomTraits())
        : Base(std::pow(std::sin(minAngleDeg * M_PI / 180.0), 2.0), traits),
          maxEdgeFactor_(maxEdgeFactor),
          targetSize_(std::move(targetSize)),
          traits_(traits) {}

    class IsBad {
    public:
        using Point2 = typename CDT::Point;

        IsBad(double aspectBound, double maxEdgeFactor,
              std::function<double(double, double)> targetSize,
              const GeomTraits& traits)
            : aspectBound_(aspectBound),
              maxEdgeFactor_(maxEdgeFactor),
              targetSize_(std::move(targetSize)),
              traits_(traits) {}

        CGAL::Mesh_2::Face_badness operator()(const Quality& q) const {
            if (q.size() > 1.0) return CGAL::Mesh_2::IMPERATIVELY_BAD;
            return (q.sine() < aspectBound_) ? CGAL::Mesh_2::BAD
                                             : CGAL::Mesh_2::NOT_BAD;
        }

        CGAL::Mesh_2::Face_badness operator()(const FaceHandle& fh,
                                              Quality& q) const {
            auto squaredDistance = traits_.compute_squared_distance_2_object();
            auto area2 = traits_.compute_area_2_object();

            const Point2& pa = fh->vertex(0)->point();
            const Point2& pb = fh->vertex(1)->point();
            const Point2& pc = fh->vertex(2)->point();

            double a = CGAL::to_double(squaredDistance(pb, pc));
            double b = CGAL::to_double(squaredDistance(pc, pa));
            double c = CGAL::to_double(squaredDistance(pa, pb));

            double maxSq = a;
            double secondMaxSq = b;
            if (maxSq < secondMaxSq) std::swap(maxSq, secondMaxSq);
            if (c > maxSq) {
                secondMaxSq = maxSq;
                maxSq = c;
            } else if (c > secondMaxSq) {
                secondMaxSq = c;
            }

            double cx = (CGAL::to_double(pa.x()) + CGAL::to_double(pb.x()) +
                         CGAL::to_double(pc.x())) / 3.0;
            double cy = (CGAL::to_double(pa.y()) + CGAL::to_double(pb.y()) +
                         CGAL::to_double(pc.y())) / 3.0;
            double h = std::max(1e-14, maxEdgeFactor_ * targetSize_(cx, cy));
            double sizeQuality = maxSq / (h * h);
            if (sizeQuality > 1.0) {
                q = Quality(1.0, sizeQuality);
                return CGAL::Mesh_2::IMPERATIVELY_BAD;
            }

            double twiceArea = 2.0 * CGAL::to_double(area2(pa, pb, pc));
            double sineQuality = (twiceArea * twiceArea) /
                                 std::max(maxSq * secondMaxSq, 1e-300);
            q = Quality(sineQuality, sizeQuality);
            return operator()(q);
        }

    private:
        double aspectBound_;
        double maxEdgeFactor_;
        std::function<double(double, double)> targetSize_;
        const GeomTraits& traits_;
    };

    IsBad is_bad_object() const {
        return IsBad(this->bound(), maxEdgeFactor_, targetSize_, traits_);
    }

private:
    double maxEdgeFactor_;
    std::function<double(double, double)> targetSize_;
    GeomTraits traits_;
};

Vector2d fromCgalPoint(const CgalPoint& p) {
    return Vector2d(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
}

void insertConstraintPolyline(CgalCdt& cdt, const std::vector<Vector2d>& polyline) {
    if (polyline.size() < 2) return;
    for (int i = 0; i + 1 < static_cast<int>(polyline.size()); ++i) {
        if ((polyline[i + 1] - polyline[i]).norm() < 1e-14) continue;
        cdt.insert_constraint(CgalPoint(polyline[i].x(), polyline[i].y()),
                              CgalPoint(polyline[i + 1].x(), polyline[i + 1].y()));
    }
}

Mesh makeCgalQualityMesh(const SolidBodyMeshSpec& cfg,
                         const std::function<double(double, double)>& signedDistance,
                         const std::function<double(double, double)>& sizeField,
                         double hNear) {
    const double seedH = std::max(hNear, 1e-14);
    const double geps = 1e-3 * std::min(cfg.hFar, seedH);
    const double snapTol = 0.08 * hNear;
    std::vector<std::vector<Vector2d>> constraintPolylines;
    std::vector<Vector2d> boundarySamples;

    auto snapBoundaryPoint = [&](const Vector2d& p) {
        for (const auto& q : boundarySamples) {
            if ((p - q).norm() <= snapTol) return q;
        }
        boundarySamples.push_back(p);
        return p;
    };
    auto addConstraint = [&](const Vector2d& a, const Vector2d& b) {
        std::vector<Vector2d> raw = sampleGradedSegment(a, b, sizeField, hNear);
        std::vector<Vector2d> snapped;
        snapped.reserve(raw.size());
        for (const auto& p : raw) {
            Vector2d q = snapBoundaryPoint(p);
            if (snapped.empty() || (snapped.back() - q).norm() > 1e-12) {
                snapped.push_back(q);
            }
        }
        if (snapped.size() >= 2) constraintPolylines.push_back(std::move(snapped));
    };

    for (const auto& segment : cfg.fixedSegments) {
        addConstraint(segment.first, segment.second);
    }
    if (cfg.addMovingSolidBoundary && cfg.solid) {
        const MatrixXd& sx = cfg.solid->currentNodes();
        for (const auto& seg : cfg.solid->movingBoundarySegments()) {
            Vector2d a = sx.row(seg.a).transpose();
            Vector2d b = sx.row(seg.b).transpose();
            addConstraint(a, b);
        }
    }
    if (constraintPolylines.empty()) {
        throw std::runtime_error("CGAL quality mesher requires constrained boundary segments");
    }

    CgalCdt cdt;
    for (const auto& polyline : constraintPolylines) insertConstraintPolyline(cdt, polyline);

    int seedCount = 0;
    if (cfg.qualityUseInteriorSeeds) {
        std::mt19937 rng(cfg.randomSeed);
        std::uniform_real_distribution<double> U(0.0, 1.0);
        const double dy = seedH * std::sqrt(3.0) / 2.0;
        const double r0max = 1.0 / (seedH * seedH);
        double yOffset = std::fmod(std::max(0.0, cfg.seedOffsetY), dy);
        double xOffset = std::fmod(std::max(0.0, cfg.seedOffsetX), seedH);
        int row = 0;
        for (double y = cfg.ya + yOffset; y <= cfg.yb + 1e-9; y += dy, ++row) {
            double xoff = (row % 2) ? 0.5 * seedH : 0.0;
            for (double x = cfg.xa + xOffset + xoff; x <= cfg.xb + 1e-9; x += seedH) {
                if (signedDistance(x, y) > -geps) continue;
                double h = std::max(sizeField(x, y), 1e-12);
                double r0 = 1.0 / (h * h);
                if (U(rng) > std::min(1.0, r0 / r0max)) continue;

                bool tooClose = false;
                Vector2d p(x, y);
                for (const auto& q : boundarySamples) {
                    if ((p - q).norm() < 0.45 * h) {
                        tooClose = true;
                        break;
                    }
                }
                if (tooClose) continue;
                cdt.insert(CgalPoint(x, y));
                ++seedCount;
            }
        }
    }

    using Criteria = LocalSizeCriteria<CgalCdt>;
    Criteria criteria(std::max(1.0, cfg.qualityMinAngleDeg),
                      std::max(1.0, cfg.qualityMaxEdgeFactor), sizeField);
    CGAL::Delaunay_mesher_2<CgalCdt, Criteria> mesher(cdt, criteria);
    mesher.init(false);

    auto markSignedDistanceDomain = [&]() {
        std::vector<CgalCdt::Face_handle> badFaces;
        auto isBad = criteria.is_bad_object();
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            Vector2d a = fromCgalPoint(fit->vertex(0)->point());
            Vector2d b = fromCgalPoint(fit->vertex(1)->point());
            Vector2d c = fromCgalPoint(fit->vertex(2)->point());
            Vector2d mid = (a + b + c) / 3.0;
            bool inside = signedDistance(mid.x(), mid.y()) < -geps;
            fit->set_in_domain(inside);
            if (!inside) continue;
            typename Criteria::Quality q;
            if (isBad(fit, q) != CGAL::Mesh_2::NOT_BAD) badFaces.push_back(fit);
        }
        mesher.set_bad_faces(badFaces.begin(), badFaces.end());
        return static_cast<int>(badFaces.size());
    };
    int badFaceCount = markSignedDistanceDomain();

    const int initialVertices = static_cast<int>(cdt.number_of_vertices());
    const int maxSteps = (cfg.qualityMaxRefineSteps > 0)
                             ? cfg.qualityMaxRefineSteps
                             : std::max(80000, 14 * initialVertices);
    int inserted = 0;
    while (!mesher.is_refinement_done()) {
        if (inserted >= maxSteps) {
            std::ostringstream oss;
            oss << "CGAL quality mesher exceeded refine-step guard"
                << " inserted=" << inserted
                << " vertices=" << cdt.number_of_vertices()
                << " initial_vertices=" << initialVertices;
            throw std::runtime_error(oss.str());
        }
        if (!mesher.step_by_step_refine_mesh()) break;
        ++inserted;
        if ((inserted % 256) == 0) badFaceCount = markSignedDistanceDomain();
    }
    badFaceCount = markSignedDistanceDomain();
    if (!mesher.is_refinement_done()) {
        throw std::runtime_error("CGAL quality mesher stopped before satisfying criteria");
    }

    Mesh mesh;
    std::unordered_map<CgalCdt::Vertex_handle, int, CGAL::Handle_hash_function> nodeId;
    std::vector<Vector2d> nodes;
    std::vector<Vector3i> elems;

    auto addVertex = [&](CgalCdt::Vertex_handle vh) {
        auto it = nodeId.find(vh);
        if (it != nodeId.end()) return it->second;
        int id = static_cast<int>(nodes.size());
        nodes.push_back(fromCgalPoint(vh->point()));
        nodeId.emplace(vh, id);
        return id;
    };

    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        if (!fit->is_in_domain()) continue;
        Vector2d a = fromCgalPoint(fit->vertex(0)->point());
        Vector2d b = fromCgalPoint(fit->vertex(1)->point());
        Vector2d c = fromCgalPoint(fit->vertex(2)->point());
        Vector2d mid = (a + b + c) / 3.0;
        if (signedDistance(mid.x(), mid.y()) >= -geps) continue;
        double area2 = (b.x() - a.x()) * (c.y() - a.y()) -
                       (b.y() - a.y()) * (c.x() - a.x());
        if (std::abs(area2) < 1e-28) continue;

        int ia = addVertex(fit->vertex(0));
        int ib = addVertex(fit->vertex(1));
        int ic = addVertex(fit->vertex(2));
        if (area2 < 0.0) std::swap(ib, ic);
        elems.emplace_back(ia, ib, ic);
    }

    if (elems.empty()) throw std::runtime_error("CGAL quality mesher produced no fluid cells");

    mesh.node.resize(static_cast<int>(nodes.size()), 2);
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        mesh.node(i, 0) = nodes[i].x();
        mesh.node(i, 1) = nodes[i].y();
    }
    mesh.elem.resize(static_cast<int>(elems.size()), 3);
    for (int i = 0; i < static_cast<int>(elems.size()); ++i) {
        mesh.elem(i, 0) = elems[i].x();
        mesh.elem(i, 1) = elems[i].y();
        mesh.elem(i, 2) = elems[i].z();
    }

    if (cfg.verbose) {
        double minA, meanA, minAr, maxAr;
        meshQuality(mesh, minA, meanA, minAr, maxAr);
        double minEdge = 1e300;
        Vector2d minAreaCentroid = Vector2d::Zero();
        Vector2d minAreaTri[3] = {Vector2d::Zero(), Vector2d::Zero(), Vector2d::Zero()};
        for (int e = 0; e < mesh.elem.rows(); ++e) {
            Vector2d p[3] = {
                mesh.node.row(mesh.elem(e, 0)).transpose(),
                mesh.node.row(mesh.elem(e, 1)).transpose(),
                mesh.node.row(mesh.elem(e, 2)).transpose()
            };
            for (int k = 0; k < 3; ++k) {
                minEdge = std::min(minEdge, (p[(k + 1) % 3] - p[k]).norm());
            }
            double area = 0.5 * std::abs((p[1].x() - p[0].x()) * (p[2].y() - p[0].y()) -
                                         (p[1].y() - p[0].y()) * (p[2].x() - p[0].x()));
            if (std::abs(area - minAr) <= std::max(1e-30, 1e-12 * minAr)) {
                minAreaCentroid = (p[0] + p[1] + p[2]) / 3.0;
                minAreaTri[0] = p[0];
                minAreaTri[1] = p[1];
                minAreaTri[2] = p[2];
            }
        }
        std::cout << "  CGAL quality mesh: boundary_samples=" << boundarySamples.size()
                  << " seeds=" << seedCount
                  << " refine_insertions=" << inserted
                  << " remaining_bad_faces=" << badFaceCount
                  << " nodes=" << mesh.node.rows()
                  << " triangles=" << mesh.elem.rows() << "\n"
                  << "  quality: min angle " << minA << " deg, mean angle " << meanA
                  << " deg; min_edge=" << minEdge
                  << " min_area_centroid=(" << minAreaCentroid.x() << ","
                  << minAreaCentroid.y() << ")"
                  << " min_area_tri=[(" << minAreaTri[0].x() << "," << minAreaTri[0].y()
                  << "),(" << minAreaTri[1].x() << "," << minAreaTri[1].y()
                  << "),(" << minAreaTri[2].x() << "," << minAreaTri[2].y() << ")]"
                  << " area in [" << minAr << ", " << maxAr << "]\n";
    }

    return mesh;
}
#endif

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
    auto coordinateSortedFallback = [&]() {
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
        auto appendSeg = [&](const SolidBoundarySegment& s, int mode) {
            Vector2d a = x.row(s.a).transpose();
            Vector2d b = x.row(s.b).transpose();
            if ((mode == SOLID_BOTTOM && a.x() > b.x()) ||
                (mode == SOLID_RIGHT && a.y() > b.y()) ||
                (mode == SOLID_TOP && a.x() < b.x()) ||
                (mode == SOLID_LEFT && a.y() < b.y())) {
                std::swap(a, b);
            }
            appendPoint(a);
            appendPoint(b);
        };
        for (const auto& s : bottom) appendSeg(s, SOLID_BOTTOM);
        for (const auto& s : right) appendSeg(s, SOLID_RIGHT);
        for (const auto& s : top) appendSeg(s, SOLID_TOP);
        for (const auto& s : left) appendSeg(s, SOLID_LEFT);
        if (!loop.empty() && (loop.front() - loop.back()).norm() < 1e-12) loop.pop_back();
        return loop;
    };

    std::unordered_map<int, std::vector<int>> adj;
    adj.reserve(solid.boundarySegments().size());
    for (const auto& s : solid.boundarySegments()) {
        if (s.a < 0 || s.b < 0 || s.a >= x.rows() || s.b >= x.rows() || s.a == s.b)
            continue;
        adj[s.a].push_back(s.b);
        adj[s.b].push_back(s.a);
    }
    if (adj.size() < 3) return coordinateSortedFallback();

    int start = -1;
    for (const auto& kv : adj) {
        int id = kv.first;
        if (kv.second.size() != 2) return coordinateSortedFallback();
        if (start < 0 ||
            x(id, 1) < x(start, 1) - 1e-12 ||
            (std::abs(x(id, 1) - x(start, 1)) <= 1e-12 && x(id, 0) < x(start, 0))) {
            start = id;
        }
    }
    if (start < 0) return coordinateSortedFallback();

    auto chooseStartNeighbor = [&]() {
        const std::vector<int>& nb = adj[start];
        int best = nb[0];
        for (int cand : nb) {
            if (x(cand, 1) < x(best, 1) - 1e-12 ||
                (std::abs(x(cand, 1) - x(best, 1)) <= 1e-12 &&
                 x(cand, 0) > x(best, 0))) {
                best = cand;
            }
        }
        return best;
    };

    std::vector<int> ids;
    ids.reserve(adj.size());
    int prev = -1;
    int curr = start;
    int next = chooseStartNeighbor();
    for (int guard = 0; guard <= static_cast<int>(adj.size()) + 2; ++guard) {
        ids.push_back(curr);
        prev = curr;
        curr = next;
        if (curr == start) break;
        auto it = adj.find(curr);
        if (it == adj.end() || it->second.size() != 2) return coordinateSortedFallback();
        next = (it->second[0] == prev) ? it->second[1] : it->second[0];
        if (next == prev) return coordinateSortedFallback();
    }
    if (curr != start || ids.size() != adj.size()) return coordinateSortedFallback();

    std::vector<Vector2d> loop;
    loop.reserve(ids.size());
    double signedArea2 = 0.0;
    for (int id : ids) loop.push_back(x.row(id).transpose());
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        const Vector2d& a = loop[i];
        const Vector2d& b = loop[(i + 1) % loop.size()];
        signedArea2 += a.x() * b.y() - a.y() * b.x();
    }
    if (std::abs(signedArea2) < 1e-18) return coordinateSortedFallback();
    if (signedArea2 < 0.0) std::reverse(loop.begin(), loop.end());
    return loop;
}

void addGradedFixedSegment(DistanceMeshSpec& spec, const Vector2d& a, const Vector2d& b,
                           const std::function<double(double, double)>& sizeFn,
                           double hMin) {
    double tol = 0.15 * hMin;
    std::vector<Vector2d> points = sampleGradedSegment(a, b, sizeFn, hMin);
    for (const auto& p : points) addFixedPointUnique(spec, p, tol);
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
    spec.seedOffsetX = cfg.seedOffsetX;
    spec.seedOffsetY = cfg.seedOffsetY;
    spec.randomSeed = cfg.randomSeed;
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
        addGradedFixedClosedLoop(spec, solidBoundaryLoop(*cfg.solid), sizeField, hNear);
    }

    if (cfg.verbose) {
        std::cout << "  body-fitted mesh size field: h_near=" << hNear
                  << " h_far=" << cfg.hFar
                  << " grade_radius=" << cfg.gradeRadius << "\n";
    }

#ifdef DG_USE_CGAL
    if (cfg.useQualityMesher) {
        return makeCgalQualityMesh(cfg, spec.signedDistance, sizeField, hNear);
    }
#endif

    Mesh mesh;
    generateDistanceMesh(mesh, spec, cfg.maxIter, cfg.verbose);
    return mesh;
}

} // namespace euler_ale
