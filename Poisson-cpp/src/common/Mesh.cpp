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
    // Map from sorted edge (u, v) to {elem_left, elem_right}
    std::map<std::pair<int, int>, std::pair<int, int>> edgeMap;
    
    int NT = elem.rows();
    for (int t = 0; t < NT; ++t) {
        // Edges: [1,2], [2,3], [3,1] (0-based: [0,1], [1,2], [2,0])
        int nodes[3] = {elem(t, 0), elem(t, 1), elem(t, 2)};
        
        for (int i = 0; i < 3; ++i) {
            int u = nodes[i];
            int v = nodes[(i + 1) % 3];
            // Edge u -> v is CCW for element t.
            
            int n1 = std::min(u, v);
            int n2 = std::max(u, v);
            std::pair<int, int> key = {n1, n2};
            
            if (edgeMap.find(key) == edgeMap.end()) {
                edgeMap[key] = {-1, -1};
            }
            
            if (u < v) {
                edgeMap[key].first = t;
            } else {
                edgeMap[key].second = t;
            }
        }
    }
    
    int NE = edgeMap.size();
    edge.resize(NE, 2);
    edge2side.resize(NE, 2);
    
    int k = 0;
    for (auto const& [key, val] : edgeMap) {
        edge(k, 0) = key.first;
        edge(k, 1) = key.second;
        edge2side(k, 0) = val.first;
        edge2side(k, 1) = val.second;
        k++;
    }
}
