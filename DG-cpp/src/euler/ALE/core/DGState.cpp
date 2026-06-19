#include "DGState.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace euler_ale {
namespace {

Vector3d barycentricCoords(const Vector2d& p, const Vector2d& a,
                           const Vector2d& b, const Vector2d& c) {
    double det = (b.y() - c.y()) * (a.x() - c.x()) +
                 (c.x() - b.x()) * (a.y() - c.y());
    if (std::abs(det) < 1e-30) return Vector3d(-1.0, -1.0, -1.0);
    double l0 = ((b.y() - c.y()) * (p.x() - c.x()) +
                 (c.x() - b.x()) * (p.y() - c.y())) / det;
    double l1 = ((c.y() - a.y()) * (p.x() - c.x()) +
                 (a.x() - c.x()) * (p.y() - c.y())) / det;
    return Vector3d(l0, l1, 1.0 - l0 - l1);
}

} // namespace

void applyPrimitiveBounds(MatrixXd& U, double rhoFloor, double pFloor, double speedMax) {
    for (int i = 0; i < U.rows(); ++i) {
        double rho = U(i, 0);
        double mx = U(i, 1);
        double my = U(i, 2);
        double E = U(i, 3);
        if (!std::isfinite(rho) || !std::isfinite(mx) ||
            !std::isfinite(my) || !std::isfinite(E)) {
            U.row(i) = euler::primToCons(rhoFloor, 0.0, 0.0, pFloor).transpose();
            continue;
        }

        rho = std::max(rho, rhoFloor);
        double u = mx / rho;
        double v = my / rho;
        double speed = std::hypot(u, v);
        if (!std::isfinite(speed)) {
            u = 0.0;
            v = 0.0;
            speed = 0.0;
        }
        if (speed > speedMax) {
            double scale = speedMax / speed;
            u *= scale;
            v *= scale;
        }

        double kinetic = 0.5 * rho * (u * u + v * v);
        double p = (euler::GAMMA - 1.0) * (E - kinetic);
        if (!std::isfinite(p) || p < pFloor) p = pFloor;
        U.row(i) = euler::primToCons(rho, u, v, p).transpose();
    }
}

MatrixXd interpolateDGToSpace(Space& oldSp, const MatrixXd& oldU, Space& newSp) {
    int oldNt = oldSp.mesh.elem.rows();
    Vector2d lo = oldSp.mesh.node.colwise().minCoeff().transpose();
    Vector2d hi = oldSp.mesh.node.colwise().maxCoeff().transpose();
    Vector2d span = (hi - lo).cwiseMax(Vector2d(1e-12, 1e-12));
    int nx = std::clamp(static_cast<int>(std::sqrt(std::max(1, oldNt))), 32, 256);
    int ny = nx;
    std::vector<std::vector<int>> bins(static_cast<size_t>(nx * ny));
    auto clampIx = [&](double x) {
        return std::clamp(static_cast<int>((x - lo.x()) / span.x() * nx), 0, nx - 1);
    };
    auto clampIy = [&](double y) {
        return std::clamp(static_cast<int>((y - lo.y()) / span.y() * ny), 0, ny - 1);
    };
    auto binId = [&](int ix, int iy) { return iy * nx + ix; };

    std::vector<Vector2d> centroid(oldNt);
    for (int e = 0; e < oldNt; ++e) {
        Vector2d a = oldSp.mesh.node.row(oldSp.mesh.elem(e, 0)).transpose();
        Vector2d b = oldSp.mesh.node.row(oldSp.mesh.elem(e, 1)).transpose();
        Vector2d c = oldSp.mesh.node.row(oldSp.mesh.elem(e, 2)).transpose();
        Vector2d eLo = a.cwiseMin(b).cwiseMin(c);
        Vector2d eHi = a.cwiseMax(b).cwiseMax(c);
        centroid[e] = (a + b + c) / 3.0;
        int ix0 = clampIx(eLo.x()), ix1 = clampIx(eHi.x());
        int iy0 = clampIy(eLo.y()), iy1 = clampIy(eHi.y());
        for (int iy = iy0; iy <= iy1; ++iy)
            for (int ix = ix0; ix <= ix1; ++ix)
                bins[binId(ix, iy)].push_back(e);
    }

    auto locate = [&](const Vector2d& p, Vector3d& lamOut) {
        int ix = clampIx(p.x());
        int iy = clampIy(p.y());
        auto tryElem = [&](int e, double tol, Vector3d& lam) {
            Vector2d a = oldSp.mesh.node.row(oldSp.mesh.elem(e, 0)).transpose();
            Vector2d b = oldSp.mesh.node.row(oldSp.mesh.elem(e, 1)).transpose();
            Vector2d c = oldSp.mesh.node.row(oldSp.mesh.elem(e, 2)).transpose();
            lam = barycentricCoords(p, a, b, c);
            return lam.minCoeff() >= -tol && lam.maxCoeff() <= 1.0 + tol;
        };
        for (int ring = 0; ring <= 3; ++ring) {
            for (int yy = std::max(0, iy - ring); yy <= std::min(ny - 1, iy + ring); ++yy) {
                for (int xx = std::max(0, ix - ring); xx <= std::min(nx - 1, ix + ring); ++xx) {
                    if (ring > 0 && std::abs(xx - ix) < ring && std::abs(yy - iy) < ring) continue;
                    for (int e : bins[binId(xx, yy)]) {
                        Vector3d lam;
                        if (tryElem(e, 1e-9, lam)) {
                            lamOut = lam.cwiseMax(0.0);
                            double s = lamOut.sum();
                            if (s < 1e-14) lamOut = Vector3d(1.0, 0.0, 0.0);
                            else lamOut /= s;
                            return e;
                        }
                    }
                }
            }
        }
        double best = 1e300;
        int bestElem = 0;
        for (int e = 0; e < oldNt; ++e) {
            double d2 = (p - centroid[e]).squaredNorm();
            if (d2 < best) {
                best = d2;
                bestElem = e;
            }
        }
        Vector2d a = oldSp.mesh.node.row(oldSp.mesh.elem(bestElem, 0)).transpose();
        Vector2d b = oldSp.mesh.node.row(oldSp.mesh.elem(bestElem, 1)).transpose();
        Vector2d c = oldSp.mesh.node.row(oldSp.mesh.elem(bestElem, 2)).transpose();
        lamOut = barycentricCoords(p, a, b, c).cwiseMax(0.0);
        double s = lamOut.sum();
        if (s < 1e-14) lamOut = Vector3d(1.0, 0.0, 0.0);
        else lamOut /= s;
        return bestElem;
    };

    MatrixXd newU = MatrixXd::Zero(newSp.nDof, oldU.cols());
    MatrixXd dofLam = newSp.fem->lagrangeNodes();
    for (int e = 0; e < newSp.mesh.elem.rows(); ++e) {
        Vector2d a = newSp.mesh.node.row(newSp.mesh.elem(e, 0)).transpose();
        Vector2d b = newSp.mesh.node.row(newSp.mesh.elem(e, 1)).transpose();
        Vector2d c = newSp.mesh.node.row(newSp.mesh.elem(e, 2)).transpose();
        for (int i = 0; i < newSp.fem->locDof; ++i) {
            Vector3d lamNew = dofLam.row(i).transpose();
            Vector2d p = lamNew(0) * a + lamNew(1) * b + lamNew(2) * c;
            Vector3d lamOld;
            int oldElem = locate(p, lamOld);
            RowVectorXd phi = oldSp.fem->computeBasisValue_all(lamOld.transpose()).row(0);
            MatrixXd Ue(oldSp.fem->locDof, oldU.cols());
            for (int j = 0; j < oldSp.fem->locDof; ++j)
                Ue.row(j) = oldU.row(oldSp.e2d(oldElem, j));
            newU.row(newSp.e2d(e, i)) = phi * Ue;
        }
    }
    return newU;
}

} // namespace euler_ale
