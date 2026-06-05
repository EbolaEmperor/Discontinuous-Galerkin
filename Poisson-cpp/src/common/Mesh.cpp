#include "Mesh.h"
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>
#include <numeric>

void Mesh::getMesh(double h) {
    int N = std::round(1.0 / h);
    int nx = N + 1;
    int ny = N + 1;
    int nNode = nx * ny;
    int nElem = 2 * N * N;

    node.resize(nNode, 2);
    elem.resize(nElem, 3);

    // Nodes
    // MATLAB: [x,y]=meshgrid(0:h:1,0:h:1); node=[x(:),y(:)];
    for (int j = 0; j < nx; ++j) { // ix
        double x = j * h;
        for (int i = 0; i < ny; ++i) { // iy
            double y = i * h;
            int idx = i + j * ny;
            node(idx, 0) = x;
            node(idx, 1) = y;
        }
    }

    // Elements
    int k = 0;
    for (int j = 0; j < N; ++j) { // ix
        for (int i = 0; i < N; ++i) { // iy
            int n1 = i + j * ny;       // (x0, y0)
            int n4 = i + 1 + j * ny;   // (x0, y1)
            int n2 = i + (j + 1) * ny; // (x1, y0)
            int n3 = i + 1 + (j + 1) * ny; // (x1, y1)
            
            elem(k, 0) = n1; elem(k, 1) = n2; elem(k, 2) = n3;
            k++;
            elem(k, 0) = n1; elem(k, 1) = n3; elem(k, 2) = n4;
            k++;
        }
    }
}

static MatrixXd sortVerticesCCW(const MatrixXd& v) {
    Vector2d center = v.colwise().mean();
    std::vector<std::pair<double, int>> angIdx;
    int n = v.rows();
    angIdx.reserve(n);
    for (int i = 0; i < n; ++i) {
        double angle = std::atan2(v(i, 1) - center(1), v(i, 0) - center(0));
        angIdx.emplace_back(angle, i);
    }
    std::sort(angIdx.begin(), angIdx.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    MatrixXd sorted(n, 2);
    for (int i = 0; i < n; ++i) {
        sorted.row(i) = v.row(angIdx[i].second);
    }
    return sorted;
}

void Mesh::getPolygonMesh(const MatrixXd& vertices, double h) {
    MatrixXd verts = sortVerticesCCW(vertices);
    int nV = verts.rows();
    Vector2d center = verts.colwise().mean();

    node.resize(nV + 1, 2);
    node.topRows(nV) = verts;
    node.row(nV) = center.transpose();

    elem.resize(nV, 3);
    for (int i = 0; i < nV; ++i) {
        elem(i, 0) = i;
        elem(i, 1) = (i + 1) % nV;
        elem(i, 2) = nV; // center
    }

    int k = static_cast<int>(std::ceil(std::log2(1.0 / h)));
    for (int i = 0; i < k; ++i) {
        uniformRefine();
    }
}

void Mesh::getHexagonMesh(const Vector2d& center, double R, double h) {
    // Base: 6 equilateral triangles fanning out from the centre (hexagon side = R).
    const int nV = 6;
    node.resize(nV + 1, 2);
    for (int i = 0; i < nV; ++i) {
        double th = 2.0 * M_PI * i / nV;
        node(i, 0) = center(0) + R * std::cos(th);
        node(i, 1) = center(1) + R * std::sin(th);
    }
    node(nV, 0) = center(0);
    node(nV, 1) = center(1);
    elem.resize(nV, 3);
    for (int i = 0; i < nV; ++i) {
        elem(i, 0) = i;
        elem(i, 1) = (i + 1) % nV;
        elem(i, 2) = nV; // centre
    }
    // Refine until the edge (= R / 2^k) is <= h.
    int k = std::max(0, static_cast<int>(std::ceil(std::log2(R / std::max(h, 1e-12)))));
    for (int i = 0; i < k; ++i) uniformRefine();
}

void Mesh::getDiskMesh(const Vector2d& center, double R, double h) {
    getHexagonMesh(center, R, h); // hexagon inscribed in the circle of radius R

    // Smoothly map the hexagon onto the disk: keep the polar angle, rescale the
    // radius by R / r_hex(theta) where r_hex is the hexagon-boundary distance at
    // that angle. This is a positive-Jacobian map (no inversions) that distributes
    // the distortion (scale in [1, 1/cos30 = 1.1547]) over the whole radius, and
    // sends the hexagon boundary exactly onto the circle.
    const double apothem = R * std::cos(M_PI / 6.0); // hexagon inradius
    for (int i = 0; i < node.rows(); ++i) {
        double dx = node(i, 0) - center(0);
        double dy = node(i, 1) - center(1);
        double r = std::sqrt(dx * dx + dy * dy);
        if (r < 1e-14) continue; // centre stays put
        double th = std::atan2(dy, dx);
        // angle to the nearest hexagon face centre (faces at 30,90,150,... degrees)
        double a = th - M_PI / 6.0;
        a -= (M_PI / 3.0) * std::round(a / (M_PI / 3.0)); // wrap to [-30,30] degrees
        double r_hex = apothem / std::cos(a);
        double s = R / r_hex;
        node(i, 0) = center(0) + dx * s;
        node(i, 1) = center(1) + dy * s;
    }
}

void Mesh::uniformRefine() {
    int NT = elem.rows();
    int oldN = node.rows();

    std::vector<Vector2d> newNodes;
    newNodes.reserve(oldN + NT * 3 / 2);
    for (int i = 0; i < oldN; ++i) {
        newNodes.push_back(node.row(i));
    }

    std::map<std::pair<int, int>, int> midMap;
    auto getMid = [&](int a, int b) {
        int u = std::min(a, b);
        int v = std::max(a, b);
        std::pair<int, int> key{u, v};
        auto it = midMap.find(key);
        if (it != midMap.end()) return it->second;
        Vector2d mid = 0.5 * (node.row(u) + node.row(v));
        int idx = static_cast<int>(newNodes.size());
        newNodes.push_back(mid);
        midMap[key] = idx;
        return idx;
    };

    MatrixXi newElem(4 * NT, 3);
    for (int t = 0; t < NT; ++t) {
        int n1 = elem(t, 0);
        int n2 = elem(t, 1);
        int n3 = elem(t, 2);

        int m12 = getMid(n1, n2);
        int m23 = getMid(n2, n3);
        int m31 = getMid(n3, n1);

        int base = 4 * t;
        newElem(base, 0) = n1;  newElem(base, 1) = m12; newElem(base, 2) = m31;
        newElem(base + 1, 0) = m12; newElem(base + 1, 1) = n2; newElem(base + 1, 2) = m23;
        newElem(base + 2, 0) = m31; newElem(base + 2, 1) = m23; newElem(base + 2, 2) = n3;
        newElem(base + 3, 0) = m12; newElem(base + 3, 1) = m23; newElem(base + 3, 2) = m31;
    }

    node.resize(static_cast<int>(newNodes.size()), 2);
    for (int i = 0; i < static_cast<int>(newNodes.size()); ++i) {
        node.row(i) = newNodes[i];
    }
    elem = newElem;
}

VectorXi Mesh::findBdryNodes(const MatrixXd& pts) {
    int N = pts.rows();
    VectorXi isBd(N);
    double tol = 1e-10;
    for (int i = 0; i < N; ++i) {
        double x = pts(i, 0);
        double y = pts(i, 1);
        if (std::abs(x) < tol || std::abs(x - 1.0) < tol ||
            std::abs(y) < tol || std::abs(y - 1.0) < tol) {
            isBd(i) = 1;
        } else {
            isBd(i) = 0;
        }
    }
    return isBd;
}

void Mesh::getEdge2Side(MatrixXi& edge, MatrixXi& edge2side) {
    // Sort-based unique-edge construction (cache-friendly; no per-edge std::map
    // node allocations).  Output is byte-identical to the previous map version:
    // edges ascending by (min,max) node id; side 0 = the element traversing the
    // edge in increasing node order (u<v, "left"), side 1 = the other ("right").
    int NT = elem.rows();
    struct Rec { int a, b, t, side; };
    std::vector<Rec> recs;
    recs.reserve(static_cast<size_t>(NT) * 3);
    for (int t = 0; t < NT; ++t) {
        int nodes[3] = {elem(t, 0), elem(t, 1), elem(t, 2)};
        for (int i = 0; i < 3; ++i) {
            int u = nodes[i], v = nodes[(i + 1) % 3];           // CCW directed edge u->v
            recs.push_back({std::min(u, v), std::max(u, v), t, (u < v) ? 0 : 1});
        }
    }
    std::sort(recs.begin(), recs.end(), [](const Rec& x, const Rec& y) {
        return x.a != y.a ? x.a < y.a : x.b < y.b;
    });
    int n = static_cast<int>(recs.size());
    int NE = 0;                                                 // count distinct edges
    for (int i = 0; i < n; ) {
        int j = i; while (j < n && recs[j].a == recs[i].a && recs[j].b == recs[i].b) ++j;
        ++NE; i = j;
    }
    edge.resize(NE, 2);
    edge2side.resize(NE, 2);
    int k = 0;
    for (int i = 0; i < n; ) {
        int j = i; while (j < n && recs[j].a == recs[i].a && recs[j].b == recs[i].b) ++j;
        edge(k, 0) = recs[i].a; edge(k, 1) = recs[i].b;
        edge2side(k, 0) = -1; edge2side(k, 1) = -1;
        for (int r = i; r < j; ++r) edge2side(k, recs[r].side) = recs[r].t;
        ++k; i = j;
    }
}
