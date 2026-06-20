#include "ElasticSolid.h"

#include "Core.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace euler_ale {
namespace {

double clamp01(double s) {
    return std::clamp(s, 0.0, 1.0);
}

double smooth01(double s) {
    s = clamp01(s);
    return s * s * (3.0 - 2.0 * s);
}

double segmentProjection(const Vector2d& p, const Vector2d& a, const Vector2d& b,
                         double& dist2) {
    Vector2d ab = b - a;
    double den = ab.squaredNorm();
    double s = 0.0;
    if (den > 1e-30) s = clamp01((p - a).dot(ab) / den);
    Vector2d q = a + s * ab;
    dist2 = (p - q).squaredNorm();
    return s;
}

double orient2dSolid(const Vector2d& a, const Vector2d& b, const Vector2d& c) {
    return (b.x() - a.x()) * (c.y() - a.y()) -
           (b.y() - a.y()) * (c.x() - a.x());
}

SolidMeshQuality computeSolidQuality(const MatrixXd& nodes, const MatrixXi& elems) {
    SolidMeshQuality q;
    if (elems.rows() == 0) return q;
    q.minAngleDeg = 180.0;
    q.meanAngleDeg = 0.0;
    q.minArea = 1e300;
    q.maxArea = 0.0;
    q.minEdge = 1e300;
    long count = 0;
    for (int e = 0; e < elems.rows(); ++e) {
        Vector2d p[3] = {
            nodes.row(elems(e, 0)).transpose(),
            nodes.row(elems(e, 1)).transpose(),
            nodes.row(elems(e, 2)).transpose()
        };
        double area = 0.5 * std::abs(orient2dSolid(p[0], p[1], p[2]));
        q.minArea = std::min(q.minArea, area);
        q.maxArea = std::max(q.maxArea, area);
        for (int k = 0; k < 3; ++k) {
            Vector2d e0 = p[(k + 1) % 3] - p[k];
            Vector2d e1 = p[(k + 2) % 3] - p[k];
            double l0 = e0.norm();
            double l1 = e1.norm();
            q.minEdge = std::min({q.minEdge, l0, l1});
            double c = e0.dot(e1) / std::max(l0 * l1, 1e-300);
            double a = std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / M_PI;
            q.minAngleDeg = std::min(q.minAngleDeg, a);
            q.meanAngleDeg += a;
            ++count;
        }
    }
    if (count > 0) q.meanAngleDeg /= count;
    if (q.minArea == 1e300) q.minArea = 0.0;
    if (q.minEdge == 1e300) q.minEdge = 0.0;
    return q;
}

void drawColorLine(std::vector<unsigned char>& image, int width, int height,
                   int x0, int y0, int x1, int y1,
                   unsigned char r, unsigned char g, unsigned char b, double blend) {
    int dx = std::abs(x1 - x0);
    int dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    while (true) {
        if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) {
            size_t k = (static_cast<size_t>(y0) * width + x0) * 3;
            image[k] = static_cast<unsigned char>(blend * r + (1.0 - blend) * image[k]);
            image[k + 1] = static_cast<unsigned char>(blend * g + (1.0 - blend) * image[k + 1]);
            image[k + 2] = static_cast<unsigned char>(blend * b + (1.0 - blend) * image[k + 2]);
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
}

} // namespace

SolidMaterial::SolidMaterial()
    : density(600.0),
      thickness(1.0),
      young(1.5e4),
      poisson(0.34),
      damping(2.4),
      velocityLimit(1.8),
      displacementLimit(0.24) {}

ElasticSolid2D::ElasticSolid2D()
    : nx_(0), ny_(0), x0_(0.0), x1_(0.0), y0_(0.0), y1_(0.0) {}

int ElasticSolid2D::nodeIndex(int i, int j) const {
    return j * (nx_ + 1) + i;
}

void ElasticSolid2D::buildRectangularBeam(double x0, double x1, double y0, double y1,
                                          int nx, int ny, const SolidMaterial& material) {
    nx_ = std::max(1, nx);
    ny_ = std::max(1, ny);
    x0_ = x0;
    x1_ = x1;
    y0_ = y0;
    y1_ = y1;
    material_ = material;

    int nNode = (nx_ + 1) * (ny_ + 1);
    X_.resize(nNode, 2);
    x_.resize(nNode, 2);
    v_ = MatrixXd::Zero(nNode, 2);
    f_ = MatrixXd::Zero(nNode, 2);
    fixed_ = VectorXi::Zero(nNode);
    mass_ = VectorXd::Zero(nNode);
    tipNodes_.clear();

    for (int j = 0; j <= ny_; ++j) {
        double y = y0_ + (y1_ - y0_) * static_cast<double>(j) / ny_;
        for (int i = 0; i <= nx_; ++i) {
            double x = x0_ + (x1_ - x0_) * static_cast<double>(i) / nx_;
            int id = nodeIndex(i, j);
            X_.row(id) << x, y;
            if (j == 0) fixed_(id) = 1;
            if (j == ny_) tipNodes_.push_back(id);
        }
    }
    x_ = X_;

    elem_.resize(2 * nx_ * ny_, 3);
    int e = 0;
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int a = nodeIndex(i, j);
            int b = nodeIndex(i + 1, j);
            int c = nodeIndex(i, j + 1);
            int d = nodeIndex(i + 1, j + 1);
            if ((i + j) % 2 == 0) {
                elem_.row(e++) << a, b, d;
                elem_.row(e++) << a, d, c;
            } else {
                elem_.row(e++) << a, b, c;
                elem_.row(e++) << b, d, c;
            }
        }
    }

    appendBoundarySegments();
    assembleStiffness();
}

void ElasticSolid2D::buildRoundedRootBeam(double x0, double x1, double y0, double y1,
                                          double rootRadius, int nx, int ny,
                                          const SolidMaterial& material) {
    double width = std::max(1e-14, x1 - x0);
    double height = std::max(1e-14, y1 - y0);
    double r = std::clamp(rootRadius, 0.0, 0.49 * width);
    if (r <= 1e-14) {
        buildRectangularBeam(x0, x1, y0, y1, nx, ny, material);
        return;
    }

    nx_ = std::max(1, nx);
    ny_ = std::max(1, ny);
    x0_ = x0;
    x1_ = x1;
    y0_ = y0;
    y1_ = y1;
    material_ = material;

    const double hardQualityAngle = 18.0;

    int nNode = (nx_ + 1) * (ny_ + 1);
    X_.resize(nNode, 2);
    x_ = X_;
    v_ = MatrixXd::Zero(nNode, 2);
    f_ = MatrixXd::Zero(nNode, 2);
    fixed_ = VectorXi::Zero(nNode);
    mass_ = VectorXd::Zero(nNode);
    tipNodes_.clear();

    for (int j = 0; j <= ny_; ++j) {
        double y = y0_ + height * static_cast<double>(j) / ny_;
        double xl = x0_;
        double xr = x1_;
        if (y < y0_ + r) {
            double dy = y - (y0_ + r);
            double dx = std::sqrt(std::max(0.0, r * r - dy * dy));
            xl = x0_ - r + dx;
            xr = x1_ + r - dx;
        }
        for (int i = 0; i <= nx_; ++i) {
            double a = static_cast<double>(i) / nx_;
            int id = nodeIndex(i, j);
            X_.row(id) << (1.0 - a) * xl + a * xr, y;
            if (j == 0) fixed_(id) = 1;
            if (j == ny_) tipNodes_.push_back(id);
        }
    }
    x_ = X_;

    auto triMinAngle = [&](int ia, int ib, int ic) {
        Vector2d p[3] = {
            X_.row(ia).transpose(),
            X_.row(ib).transpose(),
            X_.row(ic).transpose()
        };
        double best = 180.0;
        for (int k = 0; k < 3; ++k) {
            Vector2d e0 = p[(k + 1) % 3] - p[k];
            Vector2d e1 = p[(k + 2) % 3] - p[k];
            double c = e0.dot(e1) / std::max(e0.norm() * e1.norm(), 1e-300);
            double a = std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / M_PI;
            best = std::min(best, a);
        }
        return best;
    };

    elem_.resize(2 * nx_ * ny_, 3);
    int e = 0;
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int a = nodeIndex(i, j);
            int b = nodeIndex(i + 1, j);
            int c = nodeIndex(i, j + 1);
            int d = nodeIndex(i + 1, j + 1);
            double qAD = std::min(triMinAngle(a, b, d), triMinAngle(a, d, c));
            double qBC = std::min(triMinAngle(a, b, c), triMinAngle(b, d, c));
            if (qAD >= qBC) {
                elem_.row(e++) << a, b, d;
                elem_.row(e++) << a, d, c;
            } else {
                elem_.row(e++) << a, b, c;
                elem_.row(e++) << b, d, c;
            }
        }
    }

    SolidMeshQuality q = meshQuality();
    if (q.minAngleDeg < hardQualityAngle) {
        std::ostringstream oss;
        oss << "buildRoundedRootBeam: failed structured solid mesh quality guard"
            << " min_angle=" << q.minAngleDeg
            << " mean_angle=" << q.meanAngleDeg
            << " min_edge=" << q.minEdge
            << " elems=" << elem_.rows();
        throw std::runtime_error(oss.str());
    }

    appendBoundarySegments();
    assembleStiffness();
}

void ElasticSolid2D::appendBoundarySegments() {
    boundarySegments_.clear();
    movingBoundarySegments_.clear();
    auto add = [&](int a, int b, int side, bool moving) {
        SolidBoundarySegment s{a, b, side};
        boundarySegments_.push_back(s);
        if (moving) movingBoundarySegments_.push_back(s);
    };
    for (int j = 0; j < ny_; ++j) {
        add(nodeIndex(0, j), nodeIndex(0, j + 1), SOLID_LEFT, true);
        add(nodeIndex(nx_, j), nodeIndex(nx_, j + 1), SOLID_RIGHT, true);
    }
    for (int i = 0; i < nx_; ++i) {
        add(nodeIndex(i, ny_), nodeIndex(i + 1, ny_), SOLID_TOP, true);
        add(nodeIndex(i, 0), nodeIndex(i + 1, 0), SOLID_BOTTOM, false);
    }
}

void ElasticSolid2D::assembleStiffness() {
    std::vector<Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(elem_.rows()) * 36);

    double E = material_.young;
    double nu = std::clamp(material_.poisson, -0.95, 0.49);
    Matrix3d D = Matrix3d::Zero();
    D << 1.0, nu, 0.0,
         nu, 1.0, 0.0,
         0.0, 0.0, 0.5 * (1.0 - nu);
    D *= E / std::max(1e-14, 1.0 - nu * nu);

    for (int e = 0; e < elem_.rows(); ++e) {
        int ids[3] = {elem_(e, 0), elem_(e, 1), elem_(e, 2)};
        Vector2d p[3] = {X_.row(ids[0]).transpose(),
                         X_.row(ids[1]).transpose(),
                         X_.row(ids[2]).transpose()};
        double twiceA = (p[1] - p[0]).x() * (p[2] - p[0]).y() -
                        (p[1] - p[0]).y() * (p[2] - p[0]).x();
        double area = 0.5 * std::abs(twiceA);
        if (area <= 1e-16) continue;

        double b[3] = {p[1].y() - p[2].y(), p[2].y() - p[0].y(), p[0].y() - p[1].y()};
        double c[3] = {p[2].x() - p[1].x(), p[0].x() - p[2].x(), p[1].x() - p[0].x()};
        Matrix<double, 3, 6> B = Matrix<double, 3, 6>::Zero();
        for (int i = 0; i < 3; ++i) {
            B(0, 2 * i) = b[i];
            B(1, 2 * i + 1) = c[i];
            B(2, 2 * i) = c[i];
            B(2, 2 * i + 1) = b[i];
        }
        B /= (2.0 * area);
        Matrix<double, 6, 6> Ke = material_.thickness * area * (B.transpose() * D * B);
        for (int a = 0; a < 3; ++a) {
            double m = material_.density * material_.thickness * area / 3.0;
            mass_(ids[a]) += m;
            for (int ia = 0; ia < 2; ++ia) {
                int row = 2 * ids[a] + ia;
                for (int bnode = 0; bnode < 3; ++bnode) {
                    for (int ib = 0; ib < 2; ++ib) {
                        int col = 2 * ids[bnode] + ib;
                        triplets.emplace_back(row, col, Ke(2 * a + ia, 2 * bnode + ib));
                    }
                }
            }
        }
    }

    for (int i = 0; i < mass_.size(); ++i) {
        mass_(i) = std::max(mass_(i), 1e-14);
    }
    stiffness_.resize(2 * X_.rows(), 2 * X_.rows());
    stiffness_.setFromTriplets(triplets.begin(), triplets.end());
}

void ElasticSolid2D::clearExternalForces() {
    f_.setZero();
}

int ElasticSolid2D::closestMovingSegment(const Vector2d& point, double& s, double& dist2) const {
    int best = -1;
    dist2 = 1e300;
    s = 0.0;
    for (int i = 0; i < static_cast<int>(movingBoundarySegments_.size()); ++i) {
        const SolidBoundarySegment& seg = movingBoundarySegments_[i];
        Vector2d a = x_.row(seg.a).transpose();
        Vector2d b = x_.row(seg.b).transpose();
        double d2 = 0.0;
        double ss = segmentProjection(point, a, b, d2);
        if (d2 < dist2) {
            dist2 = d2;
            s = ss;
            best = i;
        }
    }
    return best;
}

void ElasticSolid2D::addBoundaryTractionAt(const Vector2d& point, const Vector2d& traction,
                                           double measure) {
    double s = 0.0;
    double d2 = 0.0;
    int idx = closestMovingSegment(point, s, d2);
    if (idx < 0) return;
    const SolidBoundarySegment& seg = movingBoundarySegments_[idx];
    Vector2d nodal = measure * traction;
    f_.row(seg.a) += ((1.0 - s) * nodal).transpose();
    f_.row(seg.b) += (s * nodal).transpose();
}

void ElasticSolid2D::advanceExplicit(double dt) {
    if (X_.rows() == 0) return;
    VectorXd u(2 * X_.rows());
    for (int i = 0; i < X_.rows(); ++i) {
        u(2 * i) = x_(i, 0) - X_(i, 0);
        u(2 * i + 1) = x_(i, 1) - X_(i, 1);
    }
    VectorXd ku = stiffness_ * u;
    for (int i = 0; i < X_.rows(); ++i) {
        if (fixed_(i)) {
            x_.row(i) = X_.row(i);
            v_.row(i).setZero();
            continue;
        }
        Vector2d fint(ku(2 * i), ku(2 * i + 1));
        Vector2d fdamp = material_.damping * mass_(i) * v_.row(i).transpose();
        Vector2d a = (f_.row(i).transpose() - fint - fdamp) / mass_(i);
        v_.row(i) += (dt * a).transpose();

        double speed = v_.row(i).norm();
        if (material_.velocityLimit > 0.0 && speed > material_.velocityLimit) {
            v_.row(i) *= material_.velocityLimit / speed;
        }

        x_.row(i) += (dt * v_.row(i).transpose()).transpose();
        Vector2d disp = x_.row(i).transpose() - X_.row(i).transpose();
        double dn = disp.norm();
        if (material_.displacementLimit > 0.0 && dn > material_.displacementLimit) {
            disp *= material_.displacementLimit / dn;
            x_.row(i) = (X_.row(i).transpose() + disp).transpose();
            double outward = v_.row(i).dot(disp.transpose());
            if (outward > 0.0) {
                Vector2d dunit = disp / std::max(disp.norm(), 1e-14);
                v_.row(i) -= (outward * dunit).transpose();
            }
        }
    }
}

void ElasticSolid2D::setState(const MatrixXd& nodes, const MatrixXd& velocities) {
    if (nodes.rows() == x_.rows() && nodes.cols() == x_.cols()) x_ = nodes;
    if (velocities.rows() == v_.rows() && velocities.cols() == v_.cols()) v_ = velocities;
    for (int i = 0; i < fixed_.size(); ++i) {
        if (!fixed_(i)) continue;
        x_.row(i) = X_.row(i);
        v_.row(i).setZero();
    }
}

int ElasticSolid2D::numNodes() const {
    return static_cast<int>(X_.rows());
}

int ElasticSolid2D::numElements() const {
    return static_cast<int>(elem_.rows());
}

double ElasticSolid2D::totalMass() const {
    return mass_.sum();
}

double ElasticSolid2D::stableTimeStep(double cfl) const {
    if (elem_.rows() == 0) return 1e300;
    double hmin = 1e300;
    for (int e = 0; e < elem_.rows(); ++e) {
        for (int k = 0; k < 3; ++k) {
            Vector2d a = X_.row(elem_(e, k)).transpose();
            Vector2d b = X_.row(elem_(e, (k + 1) % 3)).transpose();
            hmin = std::min(hmin, (b - a).norm());
        }
    }
    double wave = std::sqrt(material_.young / std::max(1e-14, material_.density * (1.0 - material_.poisson * material_.poisson)));
    return cfl * hmin / std::max(wave, 1e-14);
}

double ElasticSolid2D::maxNodeSpeed() const {
    double speed = 0.0;
    for (int i = 0; i < v_.rows(); ++i) speed = std::max(speed, v_.row(i).norm());
    return speed;
}

double ElasticSolid2D::tipDisplacementX() const {
    double sum = 0.0;
    int cnt = 0;
    if (!tipNodes_.empty()) {
        for (int id : tipNodes_) {
            sum += x_(id, 0) - X_(id, 0);
            ++cnt;
        }
        return cnt ? sum / cnt : 0.0;
    }
    for (int i = 0; i <= nx_; ++i) {
        int id = nodeIndex(i, ny_);
        sum += x_(id, 0) - X_(id, 0);
        ++cnt;
    }
    return cnt ? sum / cnt : 0.0;
}

double ElasticSolid2D::tipVelocityX() const {
    double sum = 0.0;
    int cnt = 0;
    if (!tipNodes_.empty()) {
        for (int id : tipNodes_) {
            sum += v_(id, 0);
            ++cnt;
        }
        return cnt ? sum / cnt : 0.0;
    }
    for (int i = 0; i <= nx_; ++i) {
        int id = nodeIndex(i, ny_);
        sum += v_(id, 0);
        ++cnt;
    }
    return cnt ? sum / cnt : 0.0;
}

SolidMeshQuality ElasticSolid2D::meshQuality() const {
    return computeSolidQuality(X_, elem_);
}

const MatrixXd& ElasticSolid2D::referenceNodes() const {
    return X_;
}

const MatrixXd& ElasticSolid2D::currentNodes() const {
    return x_;
}

const MatrixXd& ElasticSolid2D::velocities() const {
    return v_;
}

const MatrixXi& ElasticSolid2D::elements() const {
    return elem_;
}

const std::vector<SolidBoundarySegment>& ElasticSolid2D::boundarySegments() const {
    return boundarySegments_;
}

const std::vector<SolidBoundarySegment>& ElasticSolid2D::movingBoundarySegments() const {
    return movingBoundarySegments_;
}

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
