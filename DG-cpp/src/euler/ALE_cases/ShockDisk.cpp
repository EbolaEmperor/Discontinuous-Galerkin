#include "BodyFittedMesh.h"
#include "SolidALEMap.h"
#include "Checkpoint.h"
#include "Core.h"
#include "DGState.h"
#include "Movie.h"
#include "NeoHookeanSolid.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace euler_ale {
namespace {

using namespace Eigen;

constexpr double PI = 3.141592653589793238462643383279502884;

struct DiskGeom {
    double xa = 0.0;
    double xb = 1.20;
    double ya = 0.0;
    double yb = 0.80;
    double cx = 0.40;
    double cy = 0.40;
    double radius = 0.125;
    double pinRadius = 0.017;
    double spongeFraction = 0.20;

    double spongeStartX() const {
        return xa + (1.0 - spongeFraction) * (xb - xa);
    }

    double spongeWidth() const {
        return std::max(0.0, xb - spongeStartX());
    }

    double spongeCoordinate(double x) const {
        double width = spongeWidth();
        if (width <= 1e-14) return 0.0;
        return std::clamp((x - spongeStartX()) / width, 0.0, 1.0);
    }
};

double smooth01(double s) {
    s = std::clamp(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
}

std::vector<Vector2d> diskBoundarySamples(const DiskGeom& g, int n) {
    std::vector<Vector2d> loop;
    n = std::max(32, n);
    loop.reserve(n);
    for (int i = 0; i < n; ++i) {
        double th = 2.0 * PI * static_cast<double>(i) / n;
        loop.emplace_back(g.cx + g.radius * std::cos(th),
                          g.cy + g.radius * std::sin(th));
    }
    return loop;
}

double closedLoopPerimeter(const std::vector<Vector2d>& loop) {
    double perim = 0.0;
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        perim += (loop[(i + 1) % loop.size()] - loop[i]).norm();
    }
    return perim;
}

std::vector<Vector2d> resampleClosedLoopByArcLength(const std::vector<Vector2d>& loop,
                                                    int n) {
    std::vector<Vector2d> out;
    if (loop.empty() || n <= 0) return out;
    if (loop.size() == 1) {
        out.assign(n, loop.front());
        return out;
    }

    std::vector<double> segLen(loop.size(), 0.0);
    double perim = 0.0;
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        segLen[i] = (loop[(i + 1) % loop.size()] - loop[i]).norm();
        perim += segLen[i];
    }
    if (!(perim > 1e-14)) {
        out.assign(n, loop.front());
        return out;
    }

    out.reserve(n);
    int seg = 0;
    double accum = 0.0;
    for (int k = 0; k < n; ++k) {
        double target = perim * static_cast<double>(k) / static_cast<double>(n);
        while (seg + 1 < static_cast<int>(loop.size()) &&
               accum + segLen[seg] < target) {
            accum += segLen[seg++];
        }
        double s = (segLen[seg] > 1e-14) ? (target - accum) / segLen[seg] : 0.0;
        out.push_back(((1.0 - s) * loop[seg] + s * loop[(seg + 1) % loop.size()]).eval());
    }
    return out;
}

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

int edgeKeyA(uint64_t key) {
    return static_cast<int>(key >> 32);
}

int edgeKeyB(uint64_t key) {
    return static_cast<int>(key & 0xffffffffu);
}

double polygonSignedArea2(const std::vector<int>& ids, const MatrixXd& node) {
    double a2 = 0.0;
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
        const Vector2d a = node.row(ids[i]).transpose();
        const Vector2d b = node.row(ids[(i + 1) % ids.size()]).transpose();
        a2 += a.x() * b.y() - a.y() * b.x();
    }
    return a2;
}

void rotateLoopToLowerLeft(std::vector<int>& ids, const MatrixXd& node) {
    if (ids.empty()) return;
    int best = 0;
    for (int i = 1; i < static_cast<int>(ids.size()); ++i) {
        int a = ids[i];
        int b = ids[best];
        if (node(a, 1) < node(b, 1) - 1e-12 ||
            (std::abs(node(a, 1) - node(b, 1)) <= 1e-12 &&
             node(a, 0) < node(b, 0))) {
            best = i;
        }
    }
    std::rotate(ids.begin(), ids.begin() + best, ids.end());
}

std::vector<int> orderedClosedLoopIds(const std::vector<std::pair<int, int>>& edges,
                                      const MatrixXd& node) {
    std::unordered_map<int, std::vector<int>> adj;
    adj.reserve(edges.size());
    for (const auto& e : edges) {
        if (e.first < 0 || e.second < 0 || e.first >= node.rows() ||
            e.second >= node.rows() || e.first == e.second) {
            continue;
        }
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    if (adj.size() < 3) return {};

    int start = -1;
    for (const auto& kv : adj) {
        int id = kv.first;
        if (kv.second.size() != 2) return {};
        if (start < 0 || node(id, 1) < node(start, 1) - 1e-12 ||
            (std::abs(node(id, 1) - node(start, 1)) <= 1e-12 &&
             node(id, 0) < node(start, 0))) {
            start = id;
        }
    }
    if (start < 0) return {};

    const std::vector<int>& nb0 = adj[start];
    int next = nb0[0];
    if (node(nb0[1], 1) < node(next, 1) - 1e-12 ||
        (std::abs(node(nb0[1], 1) - node(next, 1)) <= 1e-12 &&
         node(nb0[1], 0) > node(next, 0))) {
        next = nb0[1];
    }

    std::vector<int> ids;
    ids.reserve(adj.size());
    int prev = -1;
    int curr = start;
    for (int guard = 0; guard <= static_cast<int>(adj.size()) + 2; ++guard) {
        ids.push_back(curr);
        int oldPrev = prev;
        prev = curr;
        curr = next;
        if (curr == start) break;
        auto it = adj.find(curr);
        if (it == adj.end() || it->second.size() != 2) return {};
        next = (it->second[0] == prev) ? it->second[1] : it->second[0];
        if (next == oldPrev) return {};
    }
    if (curr != start || ids.size() != adj.size()) return {};
    if (std::abs(polygonSignedArea2(ids, node)) < 1e-18) return {};
    if (polygonSignedArea2(ids, node) < 0.0) std::reverse(ids.begin(), ids.end());
    rotateLoopToLowerLeft(ids, node);
    return ids;
}

std::vector<int> orderedSolidBoundaryIds(const ElasticSolid2D& solid) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(solid.boundarySegments().size());
    for (const auto& s : solid.boundarySegments()) edges.push_back({s.a, s.b});
    return orderedClosedLoopIds(edges, solid.currentNodes());
}

std::vector<int> orderedInnerFluidBoundaryIds(const Mesh& mesh, const DiskGeom& g) {
    std::unordered_map<uint64_t, int> edgeCount;
    edgeCount.reserve(static_cast<size_t>(mesh.elem.rows()) * 3);
    for (int e = 0; e < mesh.elem.rows(); ++e) {
        for (int k = 0; k < 3; ++k) {
            int a = mesh.elem(e, k);
            int b = mesh.elem(e, (k + 1) % 3);
            ++edgeCount[edgeKey(a, b)];
        }
    }

    const double outerTol = 1e-9;
    auto onOuter = [&](const Vector2d& p) {
        return std::abs(p.x() - g.xa) <= outerTol ||
               std::abs(p.x() - g.xb) <= outerTol ||
               std::abs(p.y() - g.ya) <= outerTol ||
               std::abs(p.y() - g.yb) <= outerTol;
    };

    std::vector<std::pair<int, int>> innerEdges;
    for (const auto& kv : edgeCount) {
        if (kv.second != 1) continue;
        int a = edgeKeyA(kv.first);
        int b = edgeKeyB(kv.first);
        Vector2d mid = 0.5 * (mesh.node.row(a) + mesh.node.row(b)).transpose();
        if (onOuter(mid)) continue;
        innerEdges.push_back({a, b});
    }
    return orderedClosedLoopIds(innerEdges, mesh.node);
}

Vector2d sampleLoopAtFraction(const std::vector<Vector2d>& loop, double fraction) {
    if (loop.empty()) return Vector2d::Zero();
    if (loop.size() == 1) return loop.front();
    std::vector<double> segLen(loop.size(), 0.0);
    double perim = 0.0;
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        segLen[i] = (loop[(i + 1) % loop.size()] - loop[i]).norm();
        perim += segLen[i];
    }
    if (!(perim > 1e-14)) return loop.front();

    double target = std::clamp(fraction, 0.0, 1.0) * perim;
    double accum = 0.0;
    for (int i = 0; i < static_cast<int>(loop.size()); ++i) {
        if (accum + segLen[i] >= target || i + 1 == static_cast<int>(loop.size())) {
            double s = (segLen[i] > 1e-14) ? (target - accum) / segLen[i] : 0.0;
            return ((1.0 - s) * loop[i] + s * loop[(i + 1) % loop.size()]).eval();
        }
        accum += segLen[i];
    }
    return loop.front();
}

double maxAbsMatrixDiff(const MatrixXd& a, const MatrixXd& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return 1e300;
    double m = 0.0;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            m = std::max(m, std::abs(a(i, j) - b(i, j)));
    return m;
}

bool sameElementMatrix(const MatrixXi& a, const MatrixXi& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            if (a(i, j) != b(i, j)) return false;
    return true;
}

void buildSpaceOnMesh(const Mesh& mesh, int ord, Space& sp) {
    sp.mesh = mesh;
    sp.fem = std::make_unique<FEM>(ord, sp.mesh, false);
    sp.fem->getDOF(sp.mesh, sp.e2d, sp.nDof);
    sp.mesh.getEdge2Side(sp.edge, sp.e2s);
    sp.minH = 1e300;
    for (int k = 0; k < sp.mesh.elem.rows(); ++k) {
        sp.minH = std::min(sp.minH, hCFL(sp.mesh, k));
    }
    sp.tag = VectorXi::Zero(sp.edge.rows());
    sp.volumeVelocityHints.clear();
    sp.edgeVelocityHints.clear();
}

bool recoverAleReferenceNodesFromFluidBase(const Mesh& fluidBase,
                                           const ElasticSolid2D& solid,
                                           const DiskGeom& g,
                                           MatrixXd& recovered,
                                           double* maxBoundaryShift) {
    std::vector<int> solidIds = orderedSolidBoundaryIds(solid);
    std::vector<int> fluidIds = orderedInnerFluidBoundaryIds(fluidBase, g);
    if (solidIds.size() < 8 || fluidIds.size() < 8) return false;

    std::vector<Vector2d> solidLoop;
    solidLoop.reserve(solidIds.size());
    for (int id : solidIds) solidLoop.push_back(solid.currentNodes().row(id).transpose());
    std::vector<Vector2d> fluidLoop;
    fluidLoop.reserve(fluidIds.size());
    for (int id : fluidIds) fluidLoop.push_back(fluidBase.node.row(id).transpose());

    std::vector<double> solidArc(solidLoop.size(), 0.0);
    double solidPerim = 0.0;
    for (int i = 0; i < static_cast<int>(solidLoop.size()); ++i) {
        solidArc[i] = solidPerim;
        solidPerim += (solidLoop[(i + 1) % solidLoop.size()] - solidLoop[i]).norm();
    }
    if (!(solidPerim > 1e-14)) return false;

    recovered = solid.currentNodes();
    double maxShift = 0.0;
    for (int k = 0; k < static_cast<int>(solidIds.size()); ++k) {
        double f = solidArc[k] / solidPerim;
        Vector2d p = sampleLoopAtFraction(fluidLoop, f);
        int id = solidIds[k];
        maxShift = std::max(maxShift, (p - recovered.row(id).transpose()).norm());
        recovered.row(id) = p.transpose();
    }
    if (maxBoundaryShift) *maxBoundaryShift = maxShift;
    return true;
}

double diskSolidDistance(const DiskGeom& g, double x, double y) {
    return std::hypot(x - g.cx, y - g.cy) - g.radius;
}

void addOuterSegments(SolidBodyMeshSpec& spec, const DiskGeom& g) {
    spec.fixedSegments.push_back({Vector2d(g.xa, g.ya), Vector2d(g.xb, g.ya)});
    spec.fixedSegments.push_back({Vector2d(g.xb, g.ya), Vector2d(g.xb, g.yb)});
    spec.fixedSegments.push_back({Vector2d(g.xb, g.yb), Vector2d(g.xa, g.yb)});
    spec.fixedSegments.push_back({Vector2d(g.xa, g.yb), Vector2d(g.xa, g.ya)});
}

Mesh makeDiskSolidReferenceMesh(const DiskGeom& g, double h, int maxIter, bool verbose) {
    DistanceMeshSpec spec;
    spec.xa = g.cx - g.radius - 0.018;
    spec.xb = g.cx + g.radius + 0.018;
    spec.ya = g.cy - g.radius - 0.018;
    spec.yb = g.cy + g.radius + 0.018;
    spec.h0 = h;
    spec.seedH = 0.85 * h;
    spec.randomSeed = 20260623u;
    spec.signedDistance = [g](double x, double y) {
        return diskSolidDistance(g, x, y);
    };
    spec.targetSize = [h](double, double) { return h; };
    spec.fixedPoints = diskBoundarySamples(g, 96);
    spec.fixedPoints.emplace_back(g.cx, g.cy);
    for (int k = 0; k < 8; ++k) {
        double th = 2.0 * PI * static_cast<double>(k) / 8.0;
        spec.fixedPoints.emplace_back(g.cx + g.pinRadius * std::cos(th),
                                      g.cy + g.pinRadius * std::sin(th));
    }

    Mesh mesh;
    generateDistanceMesh(mesh, spec, maxIter, verbose);
    return mesh;
}

Mesh makeDiskSolidMeshFromCurrentBoundary(const DiskGeom& g,
                                          const ElasticSolid2D& solid,
                                          double h, int maxIter, bool verbose) {
    std::vector<Vector2d> loop = solidBoundaryLoop(solid);
    if (loop.size() < 8) return makeDiskSolidReferenceMesh(g, h, maxIter, verbose);

    double perim = closedLoopPerimeter(loop);
    if (!(perim > 1e-10) || !std::isfinite(perim)) {
        return makeDiskSolidReferenceMesh(g, h, maxIter, verbose);
    }

    int nBoundary = std::clamp(static_cast<int>(std::ceil(perim / std::max(0.85 * h, 1e-12))),
                               80, 192);
    std::vector<Vector2d> boundary = resampleClosedLoopByArcLength(loop, nBoundary);

    Vector2d lo = boundary.front();
    Vector2d hi = boundary.front();
    for (const auto& p : boundary) {
        if (!std::isfinite(p.x()) || !std::isfinite(p.y())) {
            throw std::runtime_error("shock_disk solid remesh: non-finite current boundary");
        }
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    Vector2d span = hi - lo;
    double maxSpan = std::max(span.x(), span.y());
    if (!(maxSpan > 0.25 * g.radius) || maxSpan > 4.0 * g.radius) {
        throw std::runtime_error("shock_disk solid remesh: current boundary outside guard");
    }

    DistanceMeshSpec spec;
    double pad = std::max(2.0 * h, 0.01);
    spec.xa = lo.x() - pad;
    spec.xb = hi.x() + pad;
    spec.ya = lo.y() - pad;
    spec.yb = hi.y() + pad;
    spec.h0 = h;
    spec.seedH = 0.85 * h;
    spec.randomSeed = 20260623u;
    spec.signedDistance = [boundary](double x, double y) {
        return signedDistancePolygon(boundary, x, y);
    };
    spec.targetSize = [h](double, double) { return h; };
    spec.fixedPoints = boundary;
    spec.fixedPoints.emplace_back(g.cx, g.cy);
    for (int k = 0; k < 8; ++k) {
        double th = 2.0 * PI * static_cast<double>(k) / 8.0;
        spec.fixedPoints.emplace_back(g.cx + g.pinRadius * std::cos(th),
                                      g.cy + g.pinRadius * std::sin(th));
    }

    Mesh mesh;
    generateDistanceMesh(mesh, spec, std::min(maxIter, 8), verbose);
    return mesh;
}

Mesh makeDiskFluidMesh(const DiskGeom& g, const ElasticSolid2D& solid,
                       double h, int maxIter, bool verbose) {
    std::vector<Vector2d> loop = solidBoundaryLoop(solid);
    SolidBodyMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.hFar = h;
    spec.hNearFactor = 0.22;
    spec.gradeRadius = 0.38;
    spec.randomSeed = 20260624u;
    spec.maxIter = maxIter;
    spec.verbose = verbose;
    spec.solidDistance = [loop](double x, double y) {
        return signedDistancePolygon(loop, x, y);
    };
    spec.solid = &solid;
    spec.addMovingSolidBoundary = true;
    addOuterSegments(spec, g);
    return makeSolidBodyFittedMesh(spec);
}

int diskBoundaryTag(double x, double y, double time,
                    const SolidALEMap& map, const DiskGeom& g) {
    double dLeft = std::abs(x - g.xa);
    double dRight = std::abs(x - g.xb);
    double dBottom = std::abs(y - g.ya);
    double dTop = std::abs(y - g.yb);
    double dSolid = map.distanceToBoundary(x, y, time);

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
    if (dSolid < best) tag = TAG_MOVING_WALL;
    return tag;
}

Vector2d movingCentroid(const ElasticSolid2D& solid) {
    const MatrixXd& x = solid.currentNodes();
    const VectorXi& fixed = solid.fixedMask();
    Vector2d c = Vector2d::Zero();
    int n = 0;
    for (int i = 0; i < x.rows(); ++i) {
        if (fixed(i)) continue;
        c += x.row(i).transpose();
        ++n;
    }
    if (n == 0) return Vector2d::Zero();
    return c / static_cast<double>(n);
}

double maxDisplacement(const ElasticSolid2D& solid) {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    double best = 0.0;
    for (int i = 0; i < x.rows(); ++i) {
        best = std::max(best, (x.row(i) - X.row(i)).norm());
    }
    return best;
}

double extremeReferenceDisplacement(const ElasticSolid2D& solid,
                                    int referenceAxis, bool maximumSide,
                                    int displacementComponent) {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    int best = 0;
    for (int i = 1; i < X.rows(); ++i) {
        if (maximumSide) {
            if (X(i, referenceAxis) > X(best, referenceAxis)) best = i;
        } else {
            if (X(i, referenceAxis) < X(best, referenceAxis)) best = i;
        }
    }
    return x(best, displacementComponent) - X(best, displacementComponent);
}

double maxFixedDisplacement(const ElasticSolid2D& solid) {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    const VectorXi& fixed = solid.fixedMask();
    double best = 0.0;
    for (int i = 0; i < x.rows(); ++i) {
        if (!fixed(i)) continue;
        best = std::max(best, (x.row(i) - X.row(i)).norm());
    }
    return best;
}

void copyFrameStill(const std::string& framePath, const std::string& stillPath) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::copy_file(framePath, stillPath, fs::copy_options::overwrite_existing, ec);
    if (ec) {
        std::cerr << "Warning: cannot write still " << stillPath << ": "
                  << ec.message() << "\n";
    }
}

std::string frameName(const std::string& folder, int frame) {
    char fn[512];
    std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", folder.c_str(), frame);
    return std::string(fn);
}

struct SolidTexture {
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgb;

    bool valid() const {
        return width > 0 && height > 0 &&
               rgb.size() == static_cast<size_t>(width) * height * 3;
    }
};

std::string shellQuote(const std::string& path) {
    std::string quoted = "'";
    for (char c : path) {
        if (c == '\'') quoted += "'\\''";
        else quoted += c;
    }
    quoted += "'";
    return quoted;
}

std::string ffmpegExecutable() {
    namespace fs = std::filesystem;
    const std::array<std::string, 3> candidates = {
        "/opt/homebrew/bin/ffmpeg",
        "/usr/local/bin/ffmpeg",
        "ffmpeg"
    };
    for (const auto& candidate : candidates) {
        if (candidate.find('/') == std::string::npos || fs::exists(candidate)) {
            return candidate;
        }
    }
    return "ffmpeg";
}

SolidTexture loadSolidTexture() {
    namespace fs = std::filesystem;

    SolidTexture texture;
    fs::create_directories("out");

    fs::path sourcePng = fs::path(__FILE__).parent_path() / "milk_dragon_badge.png";
    if (!fs::exists(sourcePng)) {
        sourcePng = fs::path("DG-cpp/src/euler/ALE_cases/milk_dragon_badge.png");
    }
    fs::path texturePpm = fs::path("out") / "milk_dragon_badge.ppm";

    std::error_code ecSrc;
    std::error_code ecPpm;
    bool needsConvert = fs::exists(sourcePng) && !fs::exists(texturePpm);
    if (fs::exists(sourcePng) && fs::exists(texturePpm)) {
        auto sourceTime = fs::last_write_time(sourcePng, ecSrc);
        auto ppmTime = fs::last_write_time(texturePpm, ecPpm);
        if (!ecSrc && !ecPpm && sourceTime > ppmTime) needsConvert = true;
    }

    if (needsConvert) {
        std::string command = shellQuote(ffmpegExecutable()) +
            " -y -v error -i " + shellQuote(sourcePng.string()) +
            " -frames:v 1 -update 1 " + shellQuote(texturePpm.string());
        runOutputCommand(command);
    }

    if (!readPPM(texturePpm.string(), texture.width, texture.height, texture.rgb)) {
        std::cerr << "Warning: cannot load shock disk texture " << texturePpm
                  << "; density frames will omit the milk-dragon disk texture\n";
        texture = SolidTexture();
    } else {
        std::cout << "  shock disk texture=" << texturePpm.string()
                  << " size=" << texture.width << "x" << texture.height << "\n";
    }
    return texture;
}

double cross2(const Vector2d& a, const Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

bool barycentricCurrentTriangle(const Vector2d& p, const Vector2d& a,
                                const Vector2d& b, const Vector2d& c,
                                Vector3d& lam) {
    Vector2d ab = b - a;
    Vector2d ac = c - a;
    Vector2d ap = p - a;
    double den = cross2(ab, ac);
    if (std::abs(den) < 1e-20) return false;
    lam(1) = cross2(ap, ac) / den;
    lam(2) = cross2(ab, ap) / den;
    lam(0) = 1.0 - lam(1) - lam(2);
    return lam.minCoeff() >= -1e-8 && lam.maxCoeff() <= 1.0 + 1e-8;
}

bool nearReferenceDiskBoundary(const Vector2d& p, const DiskGeom& g) {
    double r = (p - Vector2d(g.cx, g.cy)).norm();
    double tol = std::max(0.003, 0.06 * g.radius);
    return std::abs(r - g.radius) <= tol;
}

Vector2d p2ReferenceMidpoint(const Vector2d& a, const Vector2d& b,
                             const DiskGeom& g) {
    Vector2d mid = 0.5 * (a + b);
    if (nearReferenceDiskBoundary(a, g) && nearReferenceDiskBoundary(b, g)) {
        Vector2d radial = mid - Vector2d(g.cx, g.cy);
        double len = radial.norm();
        if (len > 1e-14) {
            return Vector2d(g.cx, g.cy) + (g.radius / len) * radial;
        }
    }
    return mid;
}

std::array<unsigned char, 3> sampleTextureBilinear(const SolidTexture& texture,
                                                   double u, double v) {
    u = std::clamp(u, 0.0, 1.0);
    v = std::clamp(v, 0.0, 1.0);
    double x = u * static_cast<double>(texture.width - 1);
    double y = v * static_cast<double>(texture.height - 1);
    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));
    int x1 = std::min(x0 + 1, texture.width - 1);
    int y1 = std::min(y0 + 1, texture.height - 1);
    double fx = x - x0;
    double fy = y - y0;

    auto texel = [&](int px, int py, int c) {
        size_t k = (static_cast<size_t>(py) * texture.width + px) * 3 + c;
        return static_cast<double>(texture.rgb[k]);
    };

    std::array<unsigned char, 3> rgb{};
    for (int c = 0; c < 3; ++c) {
        double c00 = texel(x0, y0, c);
        double c10 = texel(x1, y0, c);
        double c01 = texel(x0, y1, c);
        double c11 = texel(x1, y1, c);
        double value =
            (1.0 - fx) * (1.0 - fy) * c00 +
            fx * (1.0 - fy) * c10 +
            (1.0 - fx) * fy * c01 +
            fx * fy * c11;
        rgb[c] = static_cast<unsigned char>(
            std::clamp(std::lround(value), 0l, 255l));
    }
    return rgb;
}

void overlaySolidTextureP2(std::vector<unsigned char>& image, int width, int height,
                           const ElasticSolid2D& solid, const SolidTexture& texture,
                           const DiskGeom& geom,
                           double xa, double xb, double ya, double yb) {
    if (!texture.valid() || image.size() != static_cast<size_t>(width) * height * 3) {
        return;
    }

    const MatrixXd& x = solid.currentNodes();
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXi& elem = solid.elements();
    const double textureInset = 0.028;
    const double textureScale = 1.0 - 2.0 * textureInset;

    auto pixX = [&](const Vector2d& p) {
        return (p.x() - xa) / (xb - xa) * static_cast<double>(width);
    };
    auto pixY = [&](const Vector2d& p) {
        return (yb - p.y()) / (yb - ya) * static_cast<double>(height);
    };
    auto pixelToWorld = [&](int px, int py) {
        double wx = xa + (static_cast<double>(px) + 0.5) / width * (xb - xa);
        double wy = yb - (static_cast<double>(py) + 0.5) / height * (yb - ya);
        return Vector2d(wx, wy);
    };

    for (int e = 0; e < elem.rows(); ++e) {
        int i0 = elem(e, 0);
        int i1 = elem(e, 1);
        int i2 = elem(e, 2);
        Vector2d p0 = x.row(i0).transpose();
        Vector2d p1 = x.row(i1).transpose();
        Vector2d p2 = x.row(i2).transpose();
        if (std::abs(cross2(p1 - p0, p2 - p0)) < 1e-18) continue;

        double minPx = std::min({pixX(p0), pixX(p1), pixX(p2)});
        double maxPx = std::max({pixX(p0), pixX(p1), pixX(p2)});
        double minPy = std::min({pixY(p0), pixY(p1), pixY(p2)});
        double maxPy = std::max({pixY(p0), pixY(p1), pixY(p2)});
        int xBegin = std::clamp(static_cast<int>(std::floor(minPx)) - 1, 0, width - 1);
        int xEnd = std::clamp(static_cast<int>(std::ceil(maxPx)) + 1, 0, width - 1);
        int yBegin = std::clamp(static_cast<int>(std::floor(minPy)) - 1, 0, height - 1);
        int yEnd = std::clamp(static_cast<int>(std::ceil(maxPy)) + 1, 0, height - 1);

        Vector2d X0 = X.row(i0).transpose();
        Vector2d X1 = X.row(i1).transpose();
        Vector2d X2 = X.row(i2).transpose();
        Vector2d X01 = p2ReferenceMidpoint(X0, X1, geom);
        Vector2d X12 = p2ReferenceMidpoint(X1, X2, geom);
        Vector2d X20 = p2ReferenceMidpoint(X2, X0, geom);

        for (int py = yBegin; py <= yEnd; ++py) {
            for (int px = xBegin; px <= xEnd; ++px) {
                Vector3d lam = Vector3d::Zero();
                if (!barycentricCurrentTriangle(pixelToWorld(px, py), p0, p1, p2, lam)) {
                    continue;
                }

                double l0 = lam(0);
                double l1 = lam(1);
                double l2 = lam(2);
                double n0 = l0 * (2.0 * l0 - 1.0);
                double n1 = l1 * (2.0 * l1 - 1.0);
                double n2 = l2 * (2.0 * l2 - 1.0);
                double n01 = 4.0 * l0 * l1;
                double n12 = 4.0 * l1 * l2;
                double n20 = 4.0 * l2 * l0;
                Vector2d ref = n0 * X0 + n1 * X1 + n2 * X2 +
                               n01 * X01 + n12 * X12 + n20 * X20;

                double u = 0.5 + (ref.x() - geom.cx) / (2.0 * geom.radius);
                double v = 0.5 - (ref.y() - geom.cy) / (2.0 * geom.radius);
                u = textureInset + textureScale * u;
                v = textureInset + textureScale * v;
                std::array<unsigned char, 3> rgb = sampleTextureBilinear(texture, u, v);
                size_t k = (static_cast<size_t>(py) * width + px) * 3;
                image[k] = rgb[0];
                image[k + 1] = rgb[1];
                image[k + 2] = rgb[2];
            }
        }
    }
}

} // namespace

int runShockDisk(bool quick, bool freshStart = false) {
    namespace fs = std::filesystem;

    const int ord = 1;
    const int nFrames = quick ? 28 : 1800;
    const double tEnd = quick ? 0.22 : 7.20;
    const double hFluid = quick ? 0.052 : 0.021;
    const double hSolid = quick ? 0.021 : 0.010;
    const int fluidIter = quick ? 20 : 30;
    const int solidIter = quick ? 25 : 35;
    const double cfl = quick ? 0.16 : 0.13;
    const double rhoFloor = 0.08;
    const double pFloor = 0.06;
    const double speedMax = 4.5;
    const double pExt = 1.0;
    const double remeshFluidMinAngle = quick ? 4.0 : 5.0;
    const double remeshSolidMinAngle = quick ? 6.0 : 8.0;
    const int forceSmoothPasses = quick ? 2 : 4;
    const double forceSmoothBlend = 0.55;

    DiskGeom geom;
    SolidMaterial material;
    material.density = 8.0;
    material.thickness = 1.0;
    material.young = 12.0;
    material.poisson = 0.33;
    material.damping = 0.35;
    NeoHookeanSolidModel solidModel(NeoHookeanMaterial::fromSolidMaterial(material));

    std::cout << "Shock-loaded centre-pinned elastic disk: Euler ALE + NeoHookeanSolidModel\n";
    std::cout << "  domain=[" << geom.xa << "," << geom.xb << "]x["
              << geom.ya << "," << geom.yb << "] aspect="
              << (geom.xb - geom.xa) / (geom.yb - geom.ya)
              << " sponge=[" << geom.spongeStartX() << "," << geom.xb << "]\n";
    std::cout << "  generating solid FEM mesh h=" << hSolid << "...\n";
    Mesh solidMesh = makeDiskSolidReferenceMesh(geom, hSolid, solidIter, true);
    ElasticSolid2D solid;
    solid.resetReferenceMesh(solidMesh, material);
    solid.setFixedNodesInDisk(Vector2d(geom.cx, geom.cy), geom.pinRadius, true);
    solid.setAllBoundarySegmentsMoving();
    SolidMeshQuality solidQ = solid.meshQuality();
    std::cout << "  solid nodes=" << solid.numNodes()
              << " elems=" << solid.numElements()
              << " min_angle=" << solidQ.minAngleDeg
              << " fixed_center_radius=" << geom.pinRadius
              << " mass=" << solid.totalMass()
              << " force_smooth_passes=" << forceSmoothPasses
              << " force_smooth_blend=" << forceSmoothBlend << "\n";

    SolidALEMap map;
    map.setSolid(&solid);
    map.setDomain(geom.xa, geom.xb, geom.ya, geom.yb);
    map.setInfluence(0.16, 0.12);
    map.setCurrent(0.0, solid.currentNodes(), solid.velocities());

    std::cout << "  generating body-fitted fluid mesh h=" << hFluid << "...\n";
    Mesh base = makeDiskFluidMesh(geom, solid, hFluid, fluidIter, true);
    double minAng = 0.0;
    double meanAng = 0.0;
    double minArea = 0.0;
    double maxArea = 0.0;
    meshQuality(base, minAng, meanAng, minArea, maxArea);
    std::cout << "  fluid nodes=" << base.node.rows()
              << " elems=" << base.elem.rows()
              << " min_angle=" << minAng
              << " mean_angle=" << meanAng << "\n";

    RefMapFn refMap = [&](const Vector2d& X, double time) {
        return map.refToPhys(X, time);
    };
    MaxMeshSpeedFn maxSpeed = [&](double time) {
        return map.maxMeshSpeed(time);
    };
    MeshVelocityFn meshVel = [&](double x, double y, double time) {
        return map.velocityAt(x, y, time);
    };
    Tagger tagger = [&](double x, double y, double time) {
        return diskBoundaryTag(x, y, time, map, geom);
    };
    auto driverPrim = [](double time) {
        double period = 1.3;
        double onDuration = 0.8;
        double phase = std::fmod(time, period);
        double active = (phase < onDuration) ? 1.0 : 0.0;
        double rampUp = smooth01(std::min(phase, 0.030) / 0.030);
        double rampDown = smooth01(std::clamp((onDuration - phase) / 0.030, 0.0, 1.0));
        double envelope = active * std::min(rampUp, rampDown);
        int pulseIndex = static_cast<int>(std::floor(time / period));
        double intensity = 1.0 + 0.5 * std::max(0, pulseIndex - 2);
        double rho = 1.0 + 0.80 * intensity * envelope;
        double u = 0.98 * intensity * envelope;
        double p = 1.0 + 2.35 * intensity * envelope;
        return Vector4d(rho, u, 0.0, p);
    };
    ALEBCFn bc = [&](double, double, double time, const Vector4d& Um,
                     double nx, double ny, int tag, double wn) {
        if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL) {
            return movingWallGhost(Um, nx, ny, (tag == TAG_MOVING_WALL) ? wn : 0.0);
        }
        if (tag == TAG_EXACT) {
            Vector4d pr = driverPrim(time);
            return euler::primToCons(pr(0), pr(1), pr(2), pr(3));
        }
        if (tag == TAG_OUTFLOW) {
            return characteristicPressureOutletGhost(Um, nx, ny, wn,
                                                     Vector4d(1.0, 0.0, 0.0, pExt));
        }
        return Um;
    };

    ALEAdaptiveForest forest(base, ord, 4);
    Space sp;
    rebuildSpace(forest, ord, refMap, 0.0, tagger, sp);
    MatrixXd U = euler::projectInitial(*sp.fem, sp.mesh, sp.e2d,
        [&](double x, double) {
            if (x < 0.16) return Vector4d(1.74, 0.76, 0.0, 3.22);
            return Vector4d(1.0, 0.0, 0.0, 1.0);
        });
    applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
    const Vector4d spongeReference = euler::primToCons(1.0, 0.0, 0.0, pExt);
    const double spongeSigmaMax =
        3.0 * (0.98 + std::sqrt(euler::GAMMA * pExt)) /
        std::max(geom.spongeWidth(), 1e-12);
    auto applyRightSponge = [&](MatrixXd& state, const Space& space, double dtStep) {
        if (geom.spongeWidth() <= 1e-14 || dtStep <= 0.0) return;
        MatrixXd dofLam = space.fem->lagrangeNodes();
        for (int elem = 0; elem < space.mesh.elem.rows(); ++elem) {
            Vector2d p0 = space.mesh.node.row(space.mesh.elem(elem, 0));
            Vector2d p1 = space.mesh.node.row(space.mesh.elem(elem, 1));
            Vector2d p2 = space.mesh.node.row(space.mesh.elem(elem, 2));
            for (int i = 0; i < space.fem->locDof; ++i) {
                Vector3d lam = dofLam.row(i).transpose();
                Vector2d p = lam(0) * p0 + lam(1) * p1 + lam(2) * p2;
                double s = geom.spongeCoordinate(p.x());
                if (s <= 0.0) continue;
                double ramp = s * s * (3.0 - 2.0 * s);
                double alpha = std::exp(-spongeSigmaMax * ramp * dtStep);
                int dof = space.e2d(elem, i);
                Vector4d Ui = state.row(dof).transpose();
                state.row(dof) = (spongeReference + alpha * (Ui - spongeReference)).transpose();
            }
        }
    };

    double t = 0.0;
    double nextFrame = 0.0;
    double frameDt = tEnd / std::max(1, nFrames);
    int step = 0;
    int frame = 0;
    int remeshCount = 0;

    const std::string checkpointPrefix = "shock_disk";
    std::string prefix = quick ? "shock_disk_quick" : "shock_disk";
    std::string dir = "out/" + prefix + "_frames";
    std::string dirMesh = "out/" + prefix + "_mesh_frames";
    std::string dirSch = "out/" + prefix + "_schlieren_frames";
    std::string csvPath = "out/" + prefix + "_diagnostics.csv";
    fs::create_directories("out");
    fs::create_directories(dir);
    fs::create_directories(dirMesh);
    fs::create_directories(dirSch);
    SolidTexture solidTexture = loadSolidTexture();

    std::vector<CheckpointMilestone> checkpointPlan = checkpointSchedule(quick);
    std::vector<int> checkpointDone(checkpointPlan.size(), 0);
    bool resumed = false;

    if (!freshStart) {
        std::optional<RunCheckpoint> resumeCP =
            loadLatestCheckpoint(checkpointPrefix, quick, ord, nFrames, tEnd,
                                 hFluid, solid.numNodes(), true);
        if (resumeCP.has_value()) {
            const RunCheckpoint& cp = *resumeCP;
            if (cp.solidReferenceMesh.node.rows() > 0) {
                solid.resetReferenceMesh(cp.solidReferenceMesh, material);
                solid.setFixedNodesInDisk(Vector2d(geom.cx, geom.cy), geom.pinRadius, true);
                solid.setAllBoundarySegmentsMoving();
                map.setSolid(&solid);
            }
            solid.setState(cp.solidNodes, cp.solidVelocities);
            MatrixXd aleReferenceNodes = cp.aleReferenceNodes;
            if (aleReferenceNodes.rows() == solid.numNodes() &&
                aleReferenceNodes.cols() == 2) {
                double materialDiff =
                    maxAbsMatrixDiff(aleReferenceNodes, solid.referenceNodes());
                if (cp.remeshCount > 0 && materialDiff < 1e-12) {
                    MatrixXd recoveredAleReference;
                    double maxBoundaryShift = 0.0;
                    if (recoverAleReferenceNodesFromFluidBase(cp.referenceMesh, solid, geom,
                                                              recoveredAleReference,
                                                              &maxBoundaryShift)) {
                        aleReferenceNodes = recoveredAleReference;
                        std::cout << "  recovered legacy ALE reference boundary from "
                                  << "fluid checkpoint mesh; max_boundary_shift="
                                  << maxBoundaryShift << "\n";
                    } else {
                        std::cerr << "Warning: legacy checkpoint has material ALE "
                                  << "reference nodes and boundary recovery failed\n";
                    }
                }
                map.setReferenceNodes(aleReferenceNodes);
            }
            map.setCurrent(cp.time, solid.currentNodes(), solid.velocities());
            base = cp.referenceMesh;
            forest = ALEAdaptiveForest(base, ord, 4);
            rebuildSpace(forest, ord, refMap, cp.time, tagger, sp);
            if (cp.U.rows() == sp.nDof && cp.U.cols() == 4) {
                U = cp.U;
                if (cp.currentMesh.node.rows() > 0 && cp.currentMesh.elem.rows() > 0) {
                    Space checkpointSp;
                    buildSpaceOnMesh(cp.currentMesh, ord, checkpointSp);
                    if (cp.U.rows() == checkpointSp.nDof &&
                        cp.currentMesh.node.rows() == sp.mesh.node.rows() &&
                        cp.currentMesh.node.cols() == sp.mesh.node.cols()) {
                        double meshNodeDiff =
                            maxAbsMatrixDiff(cp.currentMesh.node, sp.mesh.node);
                        bool sameElem = sameElementMatrix(cp.currentMesh.elem, sp.mesh.elem);
                        if (meshNodeDiff > 1e-11 || !sameElem) {
                            U = interpolateDGToSpace(checkpointSp, cp.U, sp);
                            std::cout << "  checkpoint flow state interpolated from "
                                      << "saved physical mesh; node_diff="
                                      << meshNodeDiff
                                      << " same_elem=" << sameElem << "\n";
                        }
                    } else if (cp.U.rows() == checkpointSp.nDof) {
                        U = interpolateDGToSpace(checkpointSp, cp.U, sp);
                        std::cout << "  checkpoint flow state interpolated from "
                                  << "saved physical mesh with resized topology\n";
                    }
                }
                applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
                t = cp.time;
                nextFrame = cp.nextFrame;
                step = cp.step;
                frame = cp.frame;
                remeshCount = cp.remeshCount;
                if (cp.milestoneDone.size() == checkpointDone.size()) {
                    checkpointDone = cp.milestoneDone;
                } else {
                    double frac = (tEnd > 0.0) ? cp.time / tEnd : 1.0;
                    for (int i = 0; i < static_cast<int>(checkpointPlan.size()); ++i)
                        checkpointDone[i] =
                            (frac + 1e-12 >= checkpointPlan[i].fraction) ? 1 : 0;
                }
                pruneFramesFrom(dir, frame);
                pruneFramesFrom(dirMesh, frame);
                pruneFramesFrom(dirSch, frame);
                trimDiagnosticsToTime(csvPath, t);
                resumed = true;
                std::cout << "  resumed from checkpoint at t=" << std::fixed
                          << std::setprecision(5) << t
                          << " step=" << step
                          << " frame=" << frame
                          << " remesh=" << remeshCount << "\n";
            } else {
                std::cerr << "Warning: checkpoint DG state incompatible; cold start\n";
            }
        }
    } else {
        std::cout << "  fresh start; checkpoint auto-resume disabled\n";
        if (!quick) {
            pruneOldCheckpoints(checkpointPrefix, quick, 0);
        }
    }

    if (!resumed) {
        clearFrameDirectory(dir);
        clearFrameDirectory(dirMesh);
        clearFrameDirectory(dirSch);
    }

    std::ofstream diag;
    if (resumed) {
        diag.open(csvPath, std::ios::app);
    } else {
        diag.open(csvPath);
        diag << "time,centroid_x,centroid_y,left_dx,right_dx,top_dy,bottom_dy,"
                "max_displacement,pin_max_displacement,max_speed,"
                "drag,lift,mean_pressure,fluid_triangles,solid_nodes,solid_elements,"
                "solid_min_angle,solid_inverted,rho_min,rho_max,min_h,strain_energy,"
                "kinetic_energy\n";
    }

    const int W = quick ? 960 : 1440;
    const int H = quick ? 640 : 960;
    const int ssaa = quick ? 1 : 2;
    const double viewXa = geom.xa;
    const double viewXb = geom.xb;
    const double viewYa = geom.ya;
    const double viewYb = geom.yb;

    auto writeFrame = [&](int idx, double time) {
        rebuildSpace(forest, ord, refMap, time, tagger, sp);
        int hiW = W * ssaa;
        int hiH = H * ssaa;
        std::vector<unsigned char> hi =
            euler::renderScalarPPMImage(*sp.fem, sp.mesh, sp.e2d, U.col(0),
                                        hiW, hiH, viewXa, viewXb, viewYa, viewYb,
                                        0.55, 3.25, euler::CM_INFERNO);
        std::vector<unsigned char> hiFlow = hi;
        overlaySolidTextureP2(hi, hiW, hiH, solid, solidTexture, geom,
                              viewXa, viewXb, viewYa, viewYb);
        int outW = hiW;
        int outH = hiH;
        std::vector<unsigned char> img = hi;
        if (ssaa > 1) img = downsampleImage(hi, hiW, hiH, ssaa, outW, outH);
        writePPM(frameName(dir, idx), outW, outH, img);

        int meshW = hiW;
        int meshH = hiH;
        std::vector<unsigned char> meshImg = hiFlow;
        if (ssaa > 1) meshImg = downsampleImage(hiFlow, hiW, hiH, ssaa, meshW, meshH);
        overlayMesh(meshImg, meshW, meshH, sp.mesh, viewXa, viewXb, viewYa, viewYb);
        overlaySolidMesh(meshImg, meshW, meshH, solid, viewXa, viewXb, viewYa, viewYb, true);
        writePPM(frameName(dirMesh, idx), meshW, meshH, meshImg);

        std::string schPath = frameName(dirSch, idx);
        euler::writeSchlierenPPM(schPath, *sp.fem, sp.mesh, sp.e2d, U.col(0),
                                 outW, outH, viewXa, viewXb, viewYa, viewYb, 9.0);
        overlaySolidMesh(schPath, solid, viewXa, viewXb, viewYa, viewYb, true);
    };

    if (!resumed) {
        writeFrame(frame++, t);
        nextFrame += frameDt;
    }

    while (t < tEnd - 1e-14) {
        map.setCurrent(t, solid.currentNodes(), solid.velocities());
        rebuildSpace(forest, ord, refMap, t, tagger, sp);

        double pMean = pExt;
        double drag = 0.0;
        double lift = 0.0;
        solid.clearExternalForces();
        loadSolidFromFluidPressure(sp, U, solid, pExt, &pMean, &drag, &lift);
        solid.smoothMovingBoundaryForces(forceSmoothPasses, forceSmoothBlend);

        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        double dtSolid = solidModel.stableTimeStep(solid, 0.30);
        dt = std::min(dt, dtSolid);
        dt = std::min(dt, quick ? 0.00115 : 0.00085);
        if (!(dt > 1e-15)) {
            std::cerr << "DT-DEBUG t=" << t << " step=" << step
                      << " estDt=" << estimateDt(sp, U, maxSpeed, t, ord, cfl)
                      << " solidDt=" << dtSolid
                      << " minH=" << sp.minH
                      << " elems=" << sp.mesh.elem.rows()
                      << " Urows=" << U.rows() << " nDof=" << sp.nDof
                      << " rho_min=" << U.col(0).minCoeff()
                      << " rho_max=" << U.col(0).maxCoeff() << std::endl;
            break;
        }

        MatrixXd solidX0 = solid.currentNodes();
        MatrixXd solidV0 = solid.velocities();
        solidModel.advanceExplicit(solid, dt);
        MatrixXd solidX1 = solid.currentNodes();
        MatrixXd solidV1 = solid.velocities();
        map.setMotion(t, t + dt, solidX0, solidV0, solidX1, solidV1);
        U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
        t += dt;
        ++step;

        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        applyRightSponge(U, sp, dt);
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);

        bool solidNeedsRemesh = false;
        bool fluidNeedsRemesh = false;
        SolidMeshQuality solidQ2 = solid.currentMeshQuality();
        if (solidQ2.invertedElements > 0 || solidQ2.minAngleDeg < remeshSolidMinAngle) {
            solidNeedsRemesh = true;
        }
        if (!solidNeedsRemesh) {
            double fluidMinAng = 0.0, fluidMeanAng = 0.0, fluidMinA = 0.0, fluidMaxA = 0.0;
            meshQuality(sp.mesh, fluidMinAng, fluidMeanAng, fluidMinA, fluidMaxA);
            if (fluidMinAng < remeshFluidMinAngle) fluidNeedsRemesh = true;
        }

        if ((solidNeedsRemesh || fluidNeedsRemesh) && t < tEnd - 1e-14) {
            double fluidMinAngBefore = 0.0, fluidMeanAngBefore = 0.0, fluidMinABefore = 0.0, fluidMaxABefore = 0.0;
            meshQuality(sp.mesh, fluidMinAngBefore, fluidMeanAngBefore, fluidMinABefore, fluidMaxABefore);
            std::cout << "  coupled remesh at t=" << std::fixed << std::setprecision(4) << t
                      << " solid_need=" << solidNeedsRemesh
                      << " fluid_need=" << fluidNeedsRemesh
                      << " solid_min_angle=" << solidQ2.minAngleDeg
                      << " solid_inverted=" << solidQ2.invertedElements
                      << " fluid_min_angle=" << fluidMinAngBefore << "\n" << std::flush;
            if (solidNeedsRemesh) {
                Mesh newSolidMesh = makeDiskSolidMeshFromCurrentBoundary(geom, solid, hSolid,
                                                                         solidIter, false);
                double newSolidMinAng = 0.0, newSolidMeanAng = 0.0;
                double newSolidMinArea = 0.0, newSolidMaxArea = 0.0;
                meshQuality(newSolidMesh, newSolidMinAng, newSolidMeanAng,
                            newSolidMinArea, newSolidMaxArea);
                std::cout << "    new solid mesh: nodes=" << newSolidMesh.node.rows()
                          << " elems=" << newSolidMesh.elem.rows()
                          << " min_angle=" << newSolidMinAng
                          << " area=[" << newSolidMinArea << "," << newSolidMaxArea
                          << "]\n" << std::flush;
                if (!(newSolidMinAng > 6.0) || !(newSolidMinArea > 0.0)) {
                    throw std::runtime_error("shock_disk solid remesh: generated mesh quality guard failed");
                }
                solid.remeshToCurrentMesh(newSolidMesh, false);
                solid.setFixedNodesInDisk(Vector2d(geom.cx, geom.cy), geom.pinRadius, true);
                solid.setAllBoundarySegmentsMoving();
                map.setSolid(&solid);
            }
            map.setReferenceNodes(solid.currentNodes());
            map.setCurrent(t, solid.currentNodes(), solid.velocities());

            Space oldSp = std::move(sp);
            MatrixXd oldU = U;
            Mesh newFluidBase = makeDiskFluidMesh(geom, solid, hFluid, fluidIter, false);
            base = newFluidBase;
            forest = ALEAdaptiveForest(base, ord, 4);
            rebuildSpace(forest, ord, refMap, t, tagger, sp);
            U = interpolateDGToSpace(oldSp, oldU, sp);
            applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
            ++remeshCount;
            double fluidMinAngAfter = 0.0, fluidMeanAngAfter = 0.0, fluidMinAAfter = 0.0, fluidMaxAAfter = 0.0;
            meshQuality(sp.mesh, fluidMinAngAfter, fluidMeanAngAfter, fluidMinAAfter, fluidMaxAAfter);
            std::cout << "    new fluid mesh: nodes=" << sp.mesh.node.rows()
                      << " elems=" << sp.mesh.elem.rows()
                      << " fluid_min_angle=" << fluidMinAngAfter
                      << " solid_min_angle=" << solid.currentMeshQuality().minAngleDeg
                      << " solid_dt=" << solidModel.stableTimeStep(solid, 0.30) << "\n" << std::flush;
        }

        double rmin = U.col(0).minCoeff();
        double rmax = U.col(0).maxCoeff();
        SolidMeshQuality sq = solid.currentMeshQuality();
        Vector2d c = movingCentroid(solid);
        diag << std::setprecision(12) << t << "," << c.x() << "," << c.y()
             << "," << extremeReferenceDisplacement(solid, 0, false, 0)
             << "," << extremeReferenceDisplacement(solid, 0, true, 0)
             << "," << extremeReferenceDisplacement(solid, 1, true, 1)
             << "," << extremeReferenceDisplacement(solid, 1, false, 1)
             << "," << maxDisplacement(solid)
             << "," << maxFixedDisplacement(solid)
             << "," << solid.maxNodeSpeed()
             << "," << drag << "," << lift << "," << pMean
             << "," << sp.mesh.elem.rows()
             << "," << solid.numNodes()
             << "," << solid.numElements()
             << "," << sq.minAngleDeg
             << "," << sq.invertedElements
             << "," << rmin << "," << rmax
             << "," << sp.minH
             << "," << solidModel.strainEnergy(solid)
             << "," << solidModel.kineticEnergy(solid) << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step
                      << " leftDx=" << extremeReferenceDisplacement(solid, 0, false, 0)
                      << " rightDx=" << extremeReferenceDisplacement(solid, 0, true, 0)
                      << " maxDisp=" << maxDisplacement(solid)
                      << " pinDisp=" << maxFixedDisplacement(solid)
                      << " speed=" << solid.maxNodeSpeed()
                      << " F=(" << drag << "," << lift << ")"
                      << " pMean=" << pMean
                      << " fluidTris=" << sp.mesh.elem.rows()
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " solidMinAng=" << sq.minAngleDeg
                      << " inv=" << sq.invertedElements
                      << " remesh=" << remeshCount << "\n";
            nextFrame += frameDt;
        }

        double frac = (tEnd > 0.0) ? t / tEnd : 1.0;
        for (int i = 0; i < static_cast<int>(checkpointPlan.size()); ++i) {
            if (checkpointDone[i]) continue;
            if (frac + 1e-12 < checkpointPlan[i].fraction) continue;
            checkpointDone[i] = 1;
            diag.flush();
            RunCheckpoint cp;
            cp.quick = quick;
            cp.ord = ord;
            cp.nFrames = nFrames;
            cp.tEnd = tEnd;
            cp.h = hFluid;
            cp.time = t;
            cp.nextFrame = nextFrame;
            cp.step = step;
            cp.frame = frame;
            cp.remeshCount = remeshCount;
            cp.milestoneDone = checkpointDone;
            cp.referenceMesh = base;
            cp.currentMesh = sp.mesh;
            cp.solidReferenceMesh.node = solid.referenceNodes();
            cp.solidReferenceMesh.elem = solid.elements();
            cp.U = U;
            cp.solidNodes = solid.currentNodes();
            cp.solidVelocities = solid.velocities();
            cp.aleReferenceNodes = map.referenceNodes();
            fs::path cpPath = checkpointPath(checkpointPrefix, quick,
                                             checkpointPlan[i].label);
            if (writeCheckpointAtomic(cpPath, cp)) {
                if (!quick) pruneOldCheckpoints(checkpointPrefix, quick, 3);
                std::cout << "  checkpoint " << checkpointPlan[i].label
                          << "% written at t=" << std::fixed << std::setprecision(5)
                          << t << " -> " << cpPath.string() << "\n";
            }
        }
    }

    std::string still = "out/" + prefix + ".ppm";
    std::string stillPng = "out/" + prefix + ".png";
    std::string stillMesh = "out/" + prefix + "_mesh.ppm";
    std::string stillMeshPng = "out/" + prefix + "_mesh.png";
    std::string stillSch = "out/" + prefix + "_schlieren.ppm";
    std::string stillSchPng = "out/" + prefix + "_schlieren.png";
    copyFrameStill(frameName(dir, frame - 1), still);
    copyFrameStill(frameName(dirMesh, frame - 1), stillMesh);
    copyFrameStill(frameName(dirSch, frame - 1), stillSch);

    int fps = quick ? 30 : 60;
    std::string video = "out/" + prefix + ".mp4";
    std::string videoMesh = "out/" + prefix + "_mesh.mp4";
    std::string videoSch = "out/" + prefix + "_schlieren.mp4";
    std::string ffmpeg = shellQuote(ffmpegExecutable());
    runOutputCommand(ffmpeg + " -y -i " + shellQuote(still) +
                     " -frames:v 1 -update 1 " + shellQuote(stillPng));
    runOutputCommand(ffmpeg + " -y -i " + shellQuote(stillMesh) +
                     " -frames:v 1 -update 1 " + shellQuote(stillMeshPng));
    runOutputCommand(ffmpeg + " -y -i " + shellQuote(stillSch) +
                     " -frames:v 1 -update 1 " + shellQuote(stillSchPng));
    runOutputCommand(ffmpeg + " -y -framerate " + std::to_string(fps) +
                     " -i " + shellQuote(dir + "/frame_%05d.ppm") +
                     " -c:v libx264 -pix_fmt yuv420p -crf 16 " + shellQuote(video));
    runOutputCommand(ffmpeg + " -y -framerate " + std::to_string(fps) +
                     " -i " + shellQuote(dirMesh + "/frame_%05d.ppm") +
                     " -c:v libx264 -pix_fmt yuv420p -crf 17 " + shellQuote(videoMesh));
    runOutputCommand(ffmpeg + " -y -framerate " + std::to_string(fps) +
                     " -i " + shellQuote(dirSch + "/frame_%05d.ppm") +
                     " -c:v libx264 -pix_fmt yuv420p -crf 16 " + shellQuote(videoSch));

    std::cout << "Done. frames=" << frame
              << " density=" << video
              << " mesh=" << videoMesh
              << " schlieren=" << videoSch
              << " diagnostics=" << csvPath << "\n";
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    bool fresh = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
        if (std::string(argv[i]) == "--fresh") fresh = true;
    }
    return euler_ale::runShockDisk(quick, fresh);
}
