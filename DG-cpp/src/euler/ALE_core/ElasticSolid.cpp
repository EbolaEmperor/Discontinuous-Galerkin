#include "ElasticSolid.h"

#include "Core.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

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

uint64_t solidEdgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

int solidEdgeA(uint64_t key) {
    return static_cast<int>(key >> 32);
}

int solidEdgeB(uint64_t key) {
    return static_cast<int>(key & 0xffffffffu);
}

bool barycentricInTriangle(const Vector2d& p, const Vector2d& a, const Vector2d& b,
                           const Vector2d& c, Vector3d& lam, double tol) {
    double den = orient2dSolid(a, b, c);
    if (std::abs(den) < 1e-30) return false;
    lam(0) = orient2dSolid(p, b, c) / den;
    lam(1) = orient2dSolid(a, p, c) / den;
    lam(2) = orient2dSolid(a, b, p) / den;
    return lam.minCoeff() >= -tol && lam.maxCoeff() <= 1.0 + tol;
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
        double signedTwiceArea = orient2dSolid(p[0], p[1], p[2]);
        if (signedTwiceArea <= 0.0) ++q.invertedElements;
        double area = 0.5 * std::abs(signedTwiceArea);
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
      damping(2.4) {}

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

void ElasticSolid2D::appendBoundarySegmentsFromMesh() {
    boundarySegments_.clear();
    movingBoundarySegments_.clear();
    tipNodes_.clear();

    SolidMeshQuality q = computeSolidQuality(X_, elem_);
    double tol = std::max(1e-10, 0.35 * std::max(q.minEdge, 1e-12));
    double centerX = 0.5 * (x0_ + x1_);

    fixed_ = VectorXi::Zero(X_.rows());
    for (int i = 0; i < X_.rows(); ++i) {
        if (X_(i, 1) <= y0_ + tol) fixed_(i) = 1;
        if (X_(i, 1) >= y1_ - tol) tipNodes_.push_back(i);
    }

    std::unordered_map<uint64_t, int> edgeCount;
    edgeCount.reserve(static_cast<size_t>(elem_.rows()) * 3);
    for (int e = 0; e < elem_.rows(); ++e) {
        for (int k = 0; k < 3; ++k) {
            int a = elem_(e, k);
            int b = elem_(e, (k + 1) % 3);
            ++edgeCount[solidEdgeKey(a, b)];
        }
    }

    auto add = [&](int a, int b, int side, bool moving) {
        SolidBoundarySegment s{a, b, side};
        boundarySegments_.push_back(s);
        if (moving) movingBoundarySegments_.push_back(s);
    };

    for (const auto& kv : edgeCount) {
        if (kv.second != 1) continue;
        int a = solidEdgeA(kv.first);
        int b = solidEdgeB(kv.first);
        Vector2d pa = X_.row(a).transpose();
        Vector2d pb = X_.row(b).transpose();
        Vector2d mid = 0.5 * (pa + pb);
        bool bottom = pa.y() <= y0_ + tol && pb.y() <= y0_ + tol;
        int side = SOLID_LEFT;
        bool moving = true;
        if (bottom) {
            side = SOLID_BOTTOM;
            moving = false;
        } else if (mid.y() >= y1_ - tol) {
            side = SOLID_TOP;
        } else if (mid.x() >= centerX) {
            side = SOLID_RIGHT;
        }
        add(a, b, side, moving);
    }
}

void ElasticSolid2D::resetReferenceMesh(const Mesh& referenceMesh,
                                        const SolidMaterial& material) {
    if (referenceMesh.node.cols() != 2 || referenceMesh.elem.cols() != 3 ||
        referenceMesh.node.rows() == 0 || referenceMesh.elem.rows() == 0) {
        throw std::runtime_error("ElasticSolid2D::resetReferenceMesh: invalid mesh");
    }

    nx_ = 0;
    ny_ = 0;
    material_ = material;
    X_ = referenceMesh.node;
    x_ = X_;
    v_ = MatrixXd::Zero(X_.rows(), 2);
    f_ = MatrixXd::Zero(X_.rows(), 2);
    elem_ = referenceMesh.elem;
    mass_ = VectorXd::Zero(X_.rows());
    fixed_ = VectorXi::Zero(X_.rows());

    x0_ = X_.col(0).minCoeff();
    x1_ = X_.col(0).maxCoeff();
    y0_ = X_.col(1).minCoeff();
    y1_ = X_.col(1).maxCoeff();

    appendBoundarySegmentsFromMesh();
    assembleStiffness();
}

void ElasticSolid2D::remeshToReferenceMesh(const Mesh& referenceMesh, bool enforceFixedNodes) {
    MatrixXd oldX = X_;
    MatrixXd oldx = x_;
    MatrixXd oldv = v_;
    MatrixXi oldElem = elem_;
    SolidMaterial material = material_;

    auto interpolateOld = [&](const Vector2d& p, MatrixXd& newX,
                              MatrixXd& newV, int row) {
        Vector3d lam = Vector3d::Zero();
        int containing = -1;
        const double tol = 1e-9;
        for (int e = 0; e < oldElem.rows(); ++e) {
            Vector2d a = oldX.row(oldElem(e, 0)).transpose();
            Vector2d b = oldX.row(oldElem(e, 1)).transpose();
            Vector2d c = oldX.row(oldElem(e, 2)).transpose();
            if (barycentricInTriangle(p, a, b, c, lam, tol)) {
                containing = e;
                break;
            }
        }

        if (containing >= 0) {
            Vector2d disp = Vector2d::Zero();
            Vector2d vel = Vector2d::Zero();
            for (int k = 0; k < 3; ++k) {
                int id = oldElem(containing, k);
                disp += lam(k) * (oldx.row(id).transpose() - oldX.row(id).transpose());
                vel += lam(k) * oldv.row(id).transpose();
            }
            newX.row(row) = (p + disp).transpose();
            newV.row(row) = vel.transpose();
            return;
        }

        int nearest = 0;
        double best = 1e300;
        for (int i = 0; i < oldX.rows(); ++i) {
            double d2 = (p - oldX.row(i).transpose()).squaredNorm();
            if (d2 < best) {
                best = d2;
                nearest = i;
            }
        }
        Vector2d disp = oldx.row(nearest).transpose() - oldX.row(nearest).transpose();
        newX.row(row) = (p + disp).transpose();
        newV.row(row) = oldv.row(nearest);
    };

    MatrixXd transferredX(referenceMesh.node.rows(), 2);
    MatrixXd transferredV(referenceMesh.node.rows(), 2);
    for (int i = 0; i < referenceMesh.node.rows(); ++i) {
        Vector2d p = referenceMesh.node.row(i).transpose();
        interpolateOld(p, transferredX, transferredV, i);
    }

    resetReferenceMesh(referenceMesh, material);
    x_ = transferredX;
    v_ = transferredV;
    if (enforceFixedNodes) {
        for (int i = 0; i < fixed_.size(); ++i) {
            if (!fixed_(i)) continue;
            x_.row(i) = X_.row(i);
            v_.row(i).setZero();
        }
    }
}

void ElasticSolid2D::remeshToCurrentMesh(const Mesh& currentMesh, bool enforceFixedNodes) {
    if (currentMesh.node.cols() != 2 || currentMesh.elem.cols() != 3 ||
        currentMesh.node.rows() == 0 || currentMesh.elem.rows() == 0) {
        throw std::runtime_error("ElasticSolid2D::remeshToCurrentMesh: invalid mesh");
    }

    MatrixXd oldX = X_;
    MatrixXd oldx = x_;
    MatrixXd oldv = v_;
    MatrixXi oldElem = elem_;
    SolidMaterial material = material_;

    auto interpolateOld = [&](const Vector2d& p, MatrixXd& newReference,
                              MatrixXd& newVelocity, int row) {
        Vector3d lam = Vector3d::Zero();
        int containing = -1;
        const double tol = 1e-9;
        for (int e = 0; e < oldElem.rows(); ++e) {
            Vector2d a = oldx.row(oldElem(e, 0)).transpose();
            Vector2d b = oldx.row(oldElem(e, 1)).transpose();
            Vector2d c = oldx.row(oldElem(e, 2)).transpose();
            if (barycentricInTriangle(p, a, b, c, lam, tol)) {
                containing = e;
                break;
            }
        }

        if (containing >= 0) {
            Vector2d ref = Vector2d::Zero();
            Vector2d vel = Vector2d::Zero();
            for (int k = 0; k < 3; ++k) {
                int id = oldElem(containing, k);
                ref += lam(k) * oldX.row(id).transpose();
                vel += lam(k) * oldv.row(id).transpose();
            }
            newReference.row(row) = ref.transpose();
            newVelocity.row(row) = vel.transpose();
            return;
        }

        int nearest = 0;
        double best = 1e300;
        for (int i = 0; i < oldx.rows(); ++i) {
            double d2 = (p - oldx.row(i).transpose()).squaredNorm();
            if (d2 < best) {
                best = d2;
                nearest = i;
            }
        }
        newReference.row(row) = oldX.row(nearest);
        newVelocity.row(row) = oldv.row(nearest);
    };

    MatrixXd transferredReference(currentMesh.node.rows(), 2);
    MatrixXd transferredVelocity(currentMesh.node.rows(), 2);
    for (int i = 0; i < currentMesh.node.rows(); ++i) {
        Vector2d p = currentMesh.node.row(i).transpose();
        interpolateOld(p, transferredReference, transferredVelocity, i);
    }

    Mesh referenceMesh;
    referenceMesh.node = transferredReference;
    referenceMesh.elem = currentMesh.elem;
    resetReferenceMesh(referenceMesh, material);
    x_ = currentMesh.node;
    v_ = transferredVelocity;
    if (enforceFixedNodes) {
        for (int i = 0; i < fixed_.size(); ++i) {
            if (!fixed_(i)) continue;
            x_.row(i) = X_.row(i);
            v_.row(i).setZero();
        }
    }
}

void ElasticSolid2D::setFixedNodesInDisk(const Vector2d& center, double radius,
                                         bool clearExisting) {
    if (clearExisting) fixed_ = VectorXi::Zero(X_.rows());
    double r2 = std::max(0.0, radius) * std::max(0.0, radius);
    for (int i = 0; i < X_.rows(); ++i) {
        Vector2d p = X_.row(i).transpose();
        if ((p - center).squaredNorm() <= r2) fixed_(i) = 1;
        if (!fixed_(i)) continue;
        x_.row(i) = X_.row(i);
        v_.row(i).setZero();
    }
}

void ElasticSolid2D::setAllBoundarySegmentsMoving() {
    movingBoundarySegments_ = boundarySegments_;
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

void ElasticSolid2D::smoothMovingBoundaryForces(int passes, double blend) {
    if (passes <= 0 || movingBoundarySegments_.empty() || f_.rows() == 0) return;
    blend = std::clamp(blend, 0.0, 1.0);
    if (blend <= 0.0) return;

    std::vector<std::vector<int>> adj(static_cast<size_t>(f_.rows()));
    std::vector<int> isBoundary(static_cast<size_t>(f_.rows()), 0);
    for (const SolidBoundarySegment& seg : movingBoundarySegments_) {
        if (seg.a < 0 || seg.b < 0 || seg.a >= f_.rows() || seg.b >= f_.rows() ||
            seg.a == seg.b) {
            continue;
        }
        adj[static_cast<size_t>(seg.a)].push_back(seg.b);
        adj[static_cast<size_t>(seg.b)].push_back(seg.a);
        isBoundary[static_cast<size_t>(seg.a)] = 1;
        isBoundary[static_cast<size_t>(seg.b)] = 1;
    }

    std::vector<int> nodes;
    nodes.reserve(movingBoundarySegments_.size());
    for (int i = 0; i < static_cast<int>(isBoundary.size()); ++i) {
        if (isBoundary[static_cast<size_t>(i)] && !adj[static_cast<size_t>(i)].empty()) {
            nodes.push_back(i);
        }
    }
    if (nodes.empty()) return;

    for (int pass = 0; pass < passes; ++pass) {
        MatrixXd old = f_;
        Vector2d totalBefore = Vector2d::Zero();
        Vector2d totalAfter = Vector2d::Zero();
        for (int id : nodes) totalBefore += old.row(id).transpose();

        for (int id : nodes) {
            Vector2d nbAvg = Vector2d::Zero();
            const std::vector<int>& nb = adj[static_cast<size_t>(id)];
            for (int j : nb) nbAvg += old.row(j).transpose();
            nbAvg /= static_cast<double>(nb.size());
            Vector2d filtered = (1.0 - blend) * old.row(id).transpose() + blend * nbAvg;
            f_.row(id) = filtered.transpose();
        }

        for (int id : nodes) totalAfter += f_.row(id).transpose();
        Vector2d correction = (totalBefore - totalAfter) / static_cast<double>(nodes.size());
        for (int id : nodes) f_.row(id) += correction.transpose();
    }
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
        x_.row(i) += (dt * v_.row(i).transpose()).transpose();
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

double ElasticSolid2D::strainEnergy() const {
    if (X_.rows() == 0) return 0.0;
    VectorXd u(2 * X_.rows());
    for (int i = 0; i < X_.rows(); ++i) {
        u(2 * i) = x_(i, 0) - X_(i, 0);
        u(2 * i + 1) = x_(i, 1) - X_(i, 1);
    }
    return 0.5 * u.dot(stiffness_ * u);
}

double ElasticSolid2D::kineticEnergy() const {
    double energy = 0.0;
    for (int i = 0; i < v_.rows(); ++i) {
        energy += 0.5 * mass_(i) * v_.row(i).squaredNorm();
    }
    return energy;
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

SolidMeshQuality ElasticSolid2D::currentMeshQuality() const {
    return computeSolidQuality(x_, elem_);
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

const MatrixXd& ElasticSolid2D::externalForces() const {
    return f_;
}

const VectorXd& ElasticSolid2D::lumpedMasses() const {
    return mass_;
}

const VectorXi& ElasticSolid2D::fixedMask() const {
    return fixed_;
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


double loadSolidFromFluidPressure(const Space& sp, const MatrixXd& U,
                                  ElasticSolid2D& solid, double pExt,
                                  double* meanPressure, double* drag, double* lift) {
    MatrixXd q1d;
    VectorXd w1d;
    sp.fem->quad1d(q1d, w1d);
    int nqe = static_cast<int>(w1d.size());
    int locDof = sp.fem->locDof;
    double fx = 0.0;
    double fy = 0.0;
    double area = 0.0;
    double pInt = 0.0;
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) != TAG_MOVING_WALL) continue;
        int tt = (sp.e2s(e, 0) >= 0) ? sp.e2s(e, 0) : sp.e2s(e, 1);
        if (tt < 0) continue;
        int n1 = sp.edge(e, 0);
        int n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, tt, n1, n2);
        Matrix<double, Dynamic, 4> Ue(locDof, 4);
        for (int i = 0; i < locDof; ++i) Ue.row(i) = U.row(sp.e2d(tt, i));
        for (int q = 0; q < nqe; ++q) {
            double l1 = q1d(q, 0);
            double l2 = q1d(q, 1);
            Vector3d lam = Vector3d::Zero();
            if (eo.et == 0) {
                lam(0) = (eo.dir == 0) ? l1 : l2;
                lam(1) = (eo.dir == 0) ? l2 : l1;
            } else if (eo.et == 1) {
                lam(1) = (eo.dir == 0) ? l1 : l2;
                lam(2) = (eo.dir == 0) ? l2 : l1;
            } else {
                lam(2) = (eo.dir == 0) ? l1 : l2;
                lam(0) = (eo.dir == 0) ? l2 : l1;
            }
            RowVectorXd phi = sp.fem->computeBasisValue_all(lam.transpose()).row(0);
            Vector4d Uq = (phi * Ue).transpose();
            double p = euler::pressure(Uq);
            Vector2d Pp = l1 * sp.mesh.node.row(n1).transpose() +
                          l2 * sp.mesh.node.row(n2).transpose();
            double ds = w1d(q) * eo.he;
            double dp = p - pExt;
            Vector2d traction = dp * eo.nout;
            solid.addBoundaryTractionAt(Pp, traction, ds);
            fx += traction.x() * ds;
            fy += traction.y() * ds;
            pInt += p * ds;
            area += ds;
        }
    }
    if (meanPressure) *meanPressure = (area > 0.0) ? pInt / area : pExt;
    if (drag) *drag = fx;
    if (lift) *lift = fy;
    return fx;
}

} // namespace euler_ale
