#include "Mesh.h"
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>

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

