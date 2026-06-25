#include "DGState.h"

#include <algorithm>
#include <cmath>
#include <queue>
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

struct CentroidKDTree {
    struct Node {
        int axis;
        int idx;
        double split;
        int left = -1;
        int right = -1;
        Vector2d boxLo;
        Vector2d boxHi;
    };

    std::vector<Node> nodes;
    const std::vector<Vector2d>* points = nullptr;

    int build(std::vector<int>& ids, int lo, int hi) {
        if (lo >= hi) return -1;
        int nodeIdx = static_cast<int>(nodes.size());
        nodes.emplace_back();
        Vector2d boxLo((*points)[ids[lo]]);
        Vector2d boxHi = boxLo;
        for (int i = lo + 1; i < hi; ++i) {
            const Vector2d& p = (*points)[ids[i]];
            boxLo = boxLo.cwiseMin(p);
            boxHi = boxHi.cwiseMax(p);
        }
        Vector2d ext = boxHi - boxLo;
        int axis = (ext.x() >= ext.y()) ? 0 : 1;
        int mid = (lo + hi) / 2;
        std::nth_element(ids.begin() + lo, ids.begin() + mid, ids.begin() + hi,
            [&](int a, int b) { return (*points)[a][axis] < (*points)[b][axis]; });
        Node& n = nodes[nodeIdx];
        n.axis = axis;
        n.idx = ids[mid];
        n.split = (*points)[ids[mid]][axis];
        n.boxLo = boxLo;
        n.boxHi = boxHi;
        n.left = build(ids, lo, mid);
        n.right = build(ids, mid + 1, hi);
        return nodeIdx;
    }

    void build(const std::vector<Vector2d>& pts) {
        points = &pts;
        nodes.clear();
        std::vector<int> ids(pts.size());
        for (int i = 0; i < static_cast<int>(pts.size()); ++i) ids[i] = i;
        if (!ids.empty()) {
            nodes.reserve(pts.size());
            build(ids, 0, static_cast<int>(ids.size()));
        }
    }

    double boxDist2(const Node& n, const Vector2d& q) const {
        double dx = std::max({n.boxLo.x() - q.x(), 0.0, q.x() - n.boxHi.x()});
        double dy = std::max({n.boxLo.y() - q.y(), 0.0, q.y() - n.boxHi.y()});
        return dx * dx + dy * dy;
    }

    void knn(const Vector2d& q, int k, std::vector<int>& out) const {
        out.clear();
        if (nodes.empty() || k <= 0) return;
        using Item = std::pair<double, int>;
        std::priority_queue<Item> heap;

        std::vector<int> stack;
        stack.reserve(64);
        stack.push_back(0);
        while (!stack.empty()) {
            int ni = stack.back();
            stack.pop_back();
            if (ni < 0) continue;
            const Node& n = nodes[ni];
            if (static_cast<int>(heap.size()) >= k && boxDist2(n, q) > heap.top().first) continue;
            const Vector2d& p = (*points)[n.idx];
            double d2 = (p - q).squaredNorm();
            if (static_cast<int>(heap.size()) < k) {
                heap.emplace(d2, n.idx);
            } else if (d2 < heap.top().first) {
                heap.pop();
                heap.emplace(d2, n.idx);
            }
            int near = n.left, far = n.right;
            if (q[n.axis] > n.split) std::swap(near, far);
            if (far >= 0) stack.push_back(far);
            if (near >= 0) stack.push_back(near);
        }
        out.resize(heap.size());
        for (int i = static_cast<int>(out.size()) - 1; i >= 0; --i) {
            out[i] = heap.top().second;
            heap.pop();
        }
    }
};

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
    std::vector<Vector2d> centroid(oldNt);
    for (int e = 0; e < oldNt; ++e) {
        Vector2d a = oldSp.mesh.node.row(oldSp.mesh.elem(e, 0)).transpose();
        Vector2d b = oldSp.mesh.node.row(oldSp.mesh.elem(e, 1)).transpose();
        Vector2d c = oldSp.mesh.node.row(oldSp.mesh.elem(e, 2)).transpose();
        centroid[e] = (a + b + c) / 3.0;
    }
    CentroidKDTree tree;
    tree.build(centroid);

    auto tryElem = [&](int e, const Vector2d& p, double tol, Vector3d& lam) {
        Vector2d a = oldSp.mesh.node.row(oldSp.mesh.elem(e, 0)).transpose();
        Vector2d b = oldSp.mesh.node.row(oldSp.mesh.elem(e, 1)).transpose();
        Vector2d c = oldSp.mesh.node.row(oldSp.mesh.elem(e, 2)).transpose();
        lam = barycentricCoords(p, a, b, c);
        return lam.minCoeff() >= -tol && lam.maxCoeff() <= 1.0 + tol;
    };

    auto locate = [&](const Vector2d& p, Vector3d& lamOut) {
        std::vector<int> nn;
        int kCap = std::min(oldNt, 64);
        tree.knn(p, kCap, nn);
        for (int e : nn) {
            Vector3d lam;
            if (tryElem(e, p, 1e-9, lam)) {
                lamOut = lam.cwiseMax(0.0);
                double s = lamOut.sum();
                if (s < 1e-14) lamOut = Vector3d(1.0, 0.0, 0.0);
                else lamOut /= s;
                return e;
            }
        }
        int bestElem = nn.empty() ? 0 : nn.front();
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
