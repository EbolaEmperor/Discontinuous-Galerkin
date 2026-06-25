#include "SolidALEMap.h"

#include "Core.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

namespace euler_ale {
namespace {

double clamp01(double s) { return std::clamp(s, 0.0, 1.0); }

double smooth01(double s) {
    s = clamp01(s);
    return s * s * (3.0 - 2.0 * s);
}

void drawColorLine(std::vector<unsigned char>& img, int W, int H,
                   int x0, int y0, int x1, int y1,
                   unsigned char r, unsigned char g, unsigned char b, double alpha) {
    int dx = std::abs(x1 - x0), dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1, err = dx + dy;
    while (true) {
        if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) {
            size_t k = (static_cast<size_t>(y0) * W + x0) * 3;
            img[k]     = static_cast<unsigned char>(alpha * r + (1.0 - alpha) * img[k]);
            img[k + 1] = static_cast<unsigned char>(alpha * g + (1.0 - alpha) * img[k + 1]);
            img[k + 2] = static_cast<unsigned char>(alpha * b + (1.0 - alpha) * img[k + 2]);
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

double segmentProjection(const Vector2d& p, const Vector2d& a, const Vector2d& b, double& dist2) {
    Vector2d ab = b - a;
    double den = ab.squaredNorm();
    double s = 0.0;
    if (den > 1e-30) s = clamp01((p - a).dot(ab) / den);
    Vector2d q = a + s * ab;
    dist2 = (p - q).squaredNorm();
    return s;
}

} // namespace

struct SolidALEMap::SegmentTree {
    struct Node {
        Vector2d lo = Vector2d::Zero();
        Vector2d hi = Vector2d::Zero();
        int left = -1;
        int right = -1;
        int begin = 0;
        int end = 0;
    };

    static constexpr int leafSize = 8;
    std::vector<int> ids;
    std::vector<Vector2d> segLo;
    std::vector<Vector2d> segHi;
    std::vector<Vector2d> center;
    std::vector<Node> nodes;

    void clear() {
        ids.clear();
        segLo.clear();
        segHi.clear();
        center.clear();
        nodes.clear();
    }

    bool empty() const {
        return nodes.empty();
    }

    double rootDistance2(const Vector2d& p) const {
        return empty() ? 1e300 : boxDistance2(p, nodes[0]);
    }

    void build(const std::vector<Vector2d>& lo, const std::vector<Vector2d>& hi,
               const std::vector<Vector2d>& c) {
        clear();
        segLo = lo;
        segHi = hi;
        center = c;
        ids.resize(segLo.size());
        std::iota(ids.begin(), ids.end(), 0);
        if (!ids.empty()) buildNode(0, static_cast<int>(ids.size()));
    }

    int buildNode(int begin, int end) {
        Node node;
        node.begin = begin;
        node.end = end;
        node.lo = Vector2d(1e300, 1e300);
        node.hi = Vector2d(-1e300, -1e300);
        for (int i = begin; i < end; ++i) {
            int id = ids[i];
            node.lo = node.lo.cwiseMin(segLo[id]);
            node.hi = node.hi.cwiseMax(segHi[id]);
        }

        int idx = static_cast<int>(nodes.size());
        nodes.push_back(node);
        if (end - begin <= leafSize) return idx;

        Vector2d ext = node.hi - node.lo;
        int axis = (ext.x() >= ext.y()) ? 0 : 1;
        int mid = begin + (end - begin) / 2;
        std::nth_element(ids.begin() + begin, ids.begin() + mid, ids.begin() + end,
                         [&](int a, int b) { return center[a](axis) < center[b](axis); });
        nodes[idx].left = buildNode(begin, mid);
        nodes[idx].right = buildNode(mid, end);
        return idx;
    }

    static double boxDistance2(const Vector2d& p, const Node& n) {
        double dx = 0.0;
        if (p.x() < n.lo.x()) dx = n.lo.x() - p.x();
        else if (p.x() > n.hi.x()) dx = p.x() - n.hi.x();
        double dy = 0.0;
        if (p.y() < n.lo.y()) dy = n.lo.y() - p.y();
        else if (p.y() > n.hi.y()) dy = p.y() - n.hi.y();
        return dx * dx + dy * dy;
    }

    template <class ExactDistance>
    void searchNode(int nodeId, const Vector2d& p, ExactDistance&& exact,
                    int& best, double& bestS, double& bestD2) const {
        const Node& n = nodes[nodeId];
        if (boxDistance2(p, n) > bestD2) return;
        if (n.left < 0 && n.right < 0) {
            for (int i = n.begin; i < n.end; ++i) {
                int id = ids[i];
                double s = 0.0;
                double d2 = 0.0;
                exact(id, s, d2);
                if (d2 < bestD2) {
                    bestD2 = d2;
                    bestS = s;
                    best = id;
                }
            }
            return;
        }

        int first = n.left;
        int second = n.right;
        double dFirst = first >= 0 ? boxDistance2(p, nodes[first]) : 1e300;
        double dSecond = second >= 0 ? boxDistance2(p, nodes[second]) : 1e300;
        if (dSecond < dFirst) {
            std::swap(first, second);
            std::swap(dFirst, dSecond);
        }
        if (first >= 0 && dFirst <= bestD2) searchNode(first, p, exact, best, bestS, bestD2);
        if (second >= 0 && dSecond <= bestD2) searchNode(second, p, exact, best, bestS, bestD2);
    }

    template <class ExactDistance>
    int closest(const Vector2d& p, ExactDistance&& exact, double& s, double& dist2) const {
        if (empty()) {
            s = 0.0;
            dist2 = 1e300;
            return -1;
        }
        int best = -1;
        s = 0.0;
        dist2 = 1e300;
        searchNode(0, p, exact, best, s, dist2);
        return best;
    }

    template <class ExactDistance>
    int closest(const Vector2d& p, ExactDistance&& exact, int hint,
                double& s, double& dist2) const {
        if (empty()) {
            s = 0.0;
            dist2 = 1e300;
            return -1;
        }
        int best = -1;
        s = 0.0;
        dist2 = 1e300;
        if (hint >= 0 && hint < static_cast<int>(segLo.size())) {
            double hintS = 0.0;
            double hintD2 = 0.0;
            exact(hint, hintS, hintD2);
            best = hint;
            s = hintS;
            dist2 = hintD2;
        }
        searchNode(0, p, exact, best, s, dist2);
        return best;
    }
};

SolidALEMap::SolidALEMap()
    : solid_(nullptr),
      xa_(0.0),
      xb_(1.0),
      ya_(0.0),
      yb_(1.0),
      normalDecay_(0.07),
      wallMargin_(0.075),
      t0_(0.0),
      t1_(1e-12),
      referenceTree_(std::make_unique<SegmentTree>()),
      currentTree_(std::make_unique<SegmentTree>()) {}

SolidALEMap::~SolidALEMap() = default;
SolidALEMap::SolidALEMap(SolidALEMap&&) noexcept = default;
SolidALEMap& SolidALEMap::operator=(SolidALEMap&&) noexcept = default;

void SolidALEMap::setSolid(const ElasticSolid2D* solid) {
    solid_ = solid;
    if (solid_) referenceNodes_ = solid_->referenceNodes();
    else referenceNodes_.resize(0, 2);
    rebuildReferenceTree();
    rebuildMotionTree();
}

void SolidALEMap::setDomain(double xa, double xb, double ya, double yb) {
    xa_ = xa;
    xb_ = xb;
    ya_ = ya;
    yb_ = yb;
}

void SolidALEMap::setInfluence(double normalDecay, double wallMargin) {
    normalDecay_ = normalDecay;
    wallMargin_ = wallMargin;
}

void SolidALEMap::setReferenceNodes(const MatrixXd& referenceNodes) {
    referenceNodes_ = referenceNodes;
    rebuildReferenceTree();
}

void SolidALEMap::setCurrent(double time, const MatrixXd& nodes, const MatrixXd& velocities) {
    setMotion(time, time, nodes, velocities, nodes, velocities);
}

void SolidALEMap::setMotion(double ta, double tb, const MatrixXd& nodes0,
                            const MatrixXd& velocities0, const MatrixXd& nodes1,
                            const MatrixXd& velocities1) {
    t0_ = ta;
    t1_ = tb;
    nodes0_ = nodes0;
    nodes1_ = nodes1;
    velocities0_ = velocities0;
    velocities1_ = velocities1;
    rebuildMotionTree();
}

const MatrixXd& SolidALEMap::referenceNodes() const {
    return referenceNodes_;
}

double SolidALEMap::alpha(double time) const {
    double span = t1_ - t0_;
    if (std::abs(span) < 1e-14) return 1.0;
    return clamp01((time - t0_) / span);
}

Vector2d SolidALEMap::nodeAt(int i, double a) const {
    return ((1.0 - a) * nodes0_.row(i) + a * nodes1_.row(i)).transpose();
}

Vector2d SolidALEMap::nodeStepVelocity(int i) const {
    double span = t1_ - t0_;
    if (std::abs(span) >= 1e-14) return ((nodes1_.row(i) - nodes0_.row(i)) / span).transpose();
    return velocities1_.row(i).transpose();
}

double SolidALEMap::influence(double x, double y, double distance) const {
    double near = std::exp(-(distance / normalDecay_) * (distance / normalDecay_));
    double xWall = std::min(x - xa_, xb_ - x);
    double yWallBottom = y - ya_;
    double yWallTop = yb_ - y;
    double wall = smooth01(xWall / wallMargin_) *
                  smooth01(yWallBottom / wallMargin_) *
                  smooth01(yWallTop / wallMargin_);
    return near * wall;
}

void SolidALEMap::rebuildReferenceTree() {
    if (!referenceTree_) return;
    if (!solid_) {
        referenceTree_->clear();
        return;
    }
    const MatrixXd& X = referenceNodes_;
    const auto& segs = solid_->movingBoundarySegments();
    if (X.rows() == 0 || segs.empty()) {
        referenceTree_->clear();
        return;
    }

    std::vector<Vector2d> lo(segs.size());
    std::vector<Vector2d> hi(segs.size());
    std::vector<Vector2d> center(segs.size());
    for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
        Vector2d a = X.row(segs[i].a).transpose();
        Vector2d b = X.row(segs[i].b).transpose();
        lo[i] = a.cwiseMin(b);
        hi[i] = a.cwiseMax(b);
        center[i] = 0.5 * (a + b);
    }
    referenceTree_->build(lo, hi, center);
}

void SolidALEMap::rebuildMotionTree() {
    if (!currentTree_) return;
    if (!solid_ || nodes0_.rows() == 0 || nodes1_.rows() == 0) {
        currentTree_->clear();
        return;
    }
    const auto& segs = solid_->movingBoundarySegments();
    std::vector<Vector2d> lo(segs.size());
    std::vector<Vector2d> hi(segs.size());
    std::vector<Vector2d> center(segs.size());
    for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
        Vector2d a0 = nodes0_.row(segs[i].a).transpose();
        Vector2d b0 = nodes0_.row(segs[i].b).transpose();
        Vector2d a1 = nodes1_.row(segs[i].a).transpose();
        Vector2d b1 = nodes1_.row(segs[i].b).transpose();
        lo[i] = a0.cwiseMin(b0).cwiseMin(a1).cwiseMin(b1);
        hi[i] = a0.cwiseMax(b0).cwiseMax(a1).cwiseMax(b1);
        center[i] = 0.25 * (a0 + b0 + a1 + b1);
    }
    currentTree_->build(lo, hi, center);
}

int SolidALEMap::closestReferenceSegment(const Vector2d& point, double& s, double& dist2) const {
    if (!solid_) return -1;
    const MatrixXd& X = referenceNodes_;
    const auto& segs = solid_->movingBoundarySegments();
    return referenceTree_->closest(point, [&](int i, double& ss, double& d2) {
        Vector2d a = X.row(segs[i].a).transpose();
        Vector2d b = X.row(segs[i].b).transpose();
        ss = segmentProjection(point, a, b, d2);
    }, s, dist2);
}

int SolidALEMap::closestCurrentSegment(const Vector2d& point, double aTime,
                                       double& s, double& dist2) const {
    if (!solid_) return -1;
    const auto& segs = solid_->movingBoundarySegments();
    return currentTree_->closest(point, [&](int i, double& ss, double& d2) {
        Vector2d a = nodeAt(segs[i].a, aTime);
        Vector2d b = nodeAt(segs[i].b, aTime);
        ss = segmentProjection(point, a, b, d2);
    }, s, dist2);
}

int SolidALEMap::closestCurrentSegment(const Vector2d& point, double aTime,
                                       int segmentHint, double& s, double& dist2) const {
    if (!solid_) return -1;
    const auto& segs = solid_->movingBoundarySegments();
    return currentTree_->closest(point, [&](int i, double& ss, double& d2) {
        Vector2d a = nodeAt(segs[i].a, aTime);
        Vector2d b = nodeAt(segs[i].b, aTime);
        ss = segmentProjection(point, a, b, d2);
    }, segmentHint, s, dist2);
}

Vector2d SolidALEMap::refToPhys(const Vector2d& Xp, double time) const {
    if (!solid_ || nodes0_.rows() == 0) return Xp;
    double cutoff = 4.0 * normalDecay_;
    if (referenceTree_->rootDistance2(Xp) > cutoff * cutoff) return Xp;
    double s = 0.0;
    double d2 = 0.0;
    int idx = closestReferenceSegment(Xp, s, d2);
    if (idx < 0) return Xp;
    const auto& seg = solid_->movingBoundarySegments()[idx];
    const MatrixXd& X = referenceNodes_;
    double a = alpha(time);
    Vector2d ref = ((1.0 - s) * X.row(seg.a) + s * X.row(seg.b)).transpose();
    Vector2d cur = (1.0 - s) * nodeAt(seg.a, a) + s * nodeAt(seg.b, a);
    double dist = std::sqrt(d2);
    double w = (dist <= 1e-10) ? 1.0 : influence(Xp.x(), Xp.y(), dist);
    return Xp + w * (cur - ref);
}

Vector2d SolidALEMap::velocityAt(double x, double y, double time) const {
    int segmentHint = -1;
    return velocityAtCached(x, y, time, segmentHint);
}

Vector2d SolidALEMap::velocityAtCached(double x, double y, double time,
                                       int& segmentHint) const {
    if (!solid_ || nodes0_.rows() == 0) return Vector2d::Zero();
    Vector2d p(x, y);
    double cutoff = 4.0 * normalDecay_;
    if (currentTree_->rootDistance2(p) > cutoff * cutoff) return Vector2d::Zero();
    double s = 0.0;
    double d2 = 0.0;
    double a = alpha(time);
    int idx = closestCurrentSegment(p, a, segmentHint, s, d2);
    if (idx < 0) return Vector2d::Zero();
    segmentHint = idx;
    const auto& seg = solid_->movingBoundarySegments()[idx];
    Vector2d va = nodeStepVelocity(seg.a);
    Vector2d vb = nodeStepVelocity(seg.b);
    Vector2d vel = (1.0 - s) * va + s * vb;
    double dist = std::sqrt(d2);
    double w = (dist <= 1e-10) ? 1.0 : influence(x, y, dist);
    return w * vel;
}

double SolidALEMap::maxMeshSpeed(double time) const {
    (void)time;
    if (!solid_ || nodes0_.rows() == 0) return 0.0;
    double speed = 0.0;
    for (const auto& seg : solid_->movingBoundarySegments()) {
        speed = std::max(speed, nodeStepVelocity(seg.a).norm());
        speed = std::max(speed, nodeStepVelocity(seg.b).norm());
        if (velocities0_.rows() == nodes0_.rows()) {
            speed = std::max(speed, velocities0_.row(seg.a).norm());
            speed = std::max(speed, velocities0_.row(seg.b).norm());
        }
        if (velocities1_.rows() == nodes1_.rows()) {
            speed = std::max(speed, velocities1_.row(seg.a).norm());
            speed = std::max(speed, velocities1_.row(seg.b).norm());
        }
    }
    return speed;
}

double SolidALEMap::distanceToBoundary(double x, double y, double time) const {
    if (!solid_ || nodes0_.rows() == 0) return 1e300;
    double s = 0.0;
    double d2 = 0.0;
    closestCurrentSegment(Vector2d(x, y), alpha(time), s, d2);
    return std::sqrt(d2);
}

void overlaySolidMesh(std::vector<unsigned char>& image, int width, int height,
                      const ElasticSolid2D& solid,
                      double xa, double xb, double ya, double yb,
                      bool drawInterior) {
    auto toPix = [&](const Vector2d& p, int& px, int& py) {
        px = static_cast<int>(std::lround((p.x() - xa) / (xb - xa) * width));
        py = static_cast<int>(std::lround((yb - p.y()) / (yb - ya) * height));
    };

    const MatrixXd& x = solid.currentNodes();
    if (drawInterior) {
        const MatrixXi& elem = solid.elements();
        for (int e = 0; e < elem.rows(); ++e) {
            for (int k = 0; k < 3; ++k) {
                Vector2d a = x.row(elem(e, k)).transpose();
                Vector2d b = x.row(elem(e, (k + 1) % 3)).transpose();
                int ax = 0, ay = 0, bx = 0, by = 0;
                toPix(a, ax, ay);
                toPix(b, bx, by);
                drawColorLine(image, width, height, ax, ay, bx, by, 245, 245, 245, 0.72);
            }
        }
    }

    for (const auto& seg : solid.movingBoundarySegments()) {
        Vector2d a = x.row(seg.a).transpose();
        Vector2d b = x.row(seg.b).transpose();
        int ax = 0, ay = 0, bx = 0, by = 0;
        toPix(a, ax, ay);
        toPix(b, bx, by);
        drawColorLine(image, width, height, ax, ay, bx, by, 255, 210, 60, 0.90);
    }
}

void overlaySolidMesh(const std::string& path, const ElasticSolid2D& solid,
                      double xa, double xb, double ya, double yb, bool drawInterior) {
    int width = 0;
    int height = 0;
    std::vector<unsigned char> image;
    if (!readPPM(path, width, height, image)) return;
    overlaySolidMesh(image, width, height, solid, xa, xb, ya, yb, drawInterior);
    writePPM(path, width, height, image);
}

} // namespace euler_ale
