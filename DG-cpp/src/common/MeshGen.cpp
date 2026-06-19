#include "MeshGen.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

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
    tris.clear();
    int N = static_cast<int>(ptsIn.size());
    if (N < 3) return;

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
    pts.emplace_back(midx - 20.0 * dmax, midy - dmax);
    pts.emplace_back(midx,               midy + 20.0 * dmax);
    pts.emplace_back(midx + 20.0 * dmax, midy - dmax);
    int s0 = N, s1 = N + 1, s2 = N + 2;

    // Scale-aware tolerance for the in-circle test (degeneracy guard).
    const double tol = 1e-12 * dmax * dmax * dmax * dmax;

    auto pushCCW = [&](std::vector<Vector3i>& out, int a, int b, int c) {
        if (orient2d(pts[a], pts[b], pts[c]) < 0) std::swap(b, c);
        out.emplace_back(a, b, c);
    };

    std::vector<Vector3i> T;
    pushCCW(T, s0, s1, s2);

    for (int ip = 0; ip < N; ++ip) {
        const Vector2d& p = pts[ip];

        // Directed edges of all "bad" triangles (those whose circumcircle holds p).
        std::map<std::pair<int, int>, int> edgeCount;
        std::vector<Vector3i> good;
        good.reserve(T.size());

        for (const auto& t : T) {
            double det = inCircleDet(pts[t[0]], pts[t[1]], pts[t[2]], p);
            if (det > tol) {
                edgeCount[{t[0], t[1]}]++;
                edgeCount[{t[1], t[2]}]++;
                edgeCount[{t[2], t[0]}]++;
            } else {
                good.push_back(t);
            }
        }

        // Cavity boundary = directed edges whose reverse is absent among bad tris.
        for (const auto& kv : edgeCount) {
            int u = kv.first.first, v = kv.first.second;
            if (edgeCount.find({v, u}) == edgeCount.end()) {
                pushCCW(good, u, v, ip);
            }
        }
        T.swap(good);
    }

    // Drop triangles that still reference a super-triangle vertex.
    tris.reserve(T.size());
    for (const auto& t : T) {
        if (t[0] >= N || t[1] >= N || t[2] >= N) continue;
        tris.push_back(t);
    }
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
    const int nFix = static_cast<int>(p.size());

    std::mt19937 rng(12345u);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double dy = seedH * std::sqrt(3.0) / 2.0;
    const double r0max = 1.0 / (seedH * seedH);
    int row = 0;
    for (double y = spec.ya; y <= spec.yb + 1e-9; y += dy, ++row) {
        double xoff = (row % 2) ? 0.5 * seedH : 0.0;
        for (double x = spec.xa + xoff; x <= spec.xb + 1e-9; x += seedH) {
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

    int Np = static_cast<int>(p.size());
    std::vector<Vector2d> pold(Np, Vector2d(1e10, 1e10));
    std::vector<Vector3i> tris;

    auto bigMove = [&]() {
        double m = 0.0;
        for (int i = 0; i < Np; ++i) m = std::max(m, (p[i] - pold[i]).norm());
        return m;
    };

    auto retriangulate = [&]() {
        std::vector<Vector3i> allTris;
        delaunayTriangulate(p, allTris);
        tris.clear();
        for (const auto& t : allTris) {
            Vector2d cmid = (p[t[0]] + p[t[1]] + p[t[2]]) / 3.0;
            if (spec.signedDistance(cmid.x(), cmid.y()) < -geps) tris.push_back(t);
        }
    };

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        if (bigMove() > ttol * hTol) {
            pold = p;
            retriangulate();
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

        if (maxInteriorMove < dptol * hTol && iter > 5) { ++iter; break; }
    }

    retriangulate();

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
