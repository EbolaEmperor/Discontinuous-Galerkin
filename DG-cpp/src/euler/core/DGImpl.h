#ifndef EULER_DG_IMPL_H
#define EULER_DG_IMPL_H

#include "DG.h"
#include "Quadrature.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <thread>
#include <vector>

namespace euler {

namespace {
template <class Worker>
void parallel_ranges(int n, const Worker& worker) {
    unsigned hw = std::thread::hardware_concurrency();
    unsigned nt = std::max(1u, std::min<unsigned>(hw ? hw : 1u, (unsigned)std::max(1, n)));
    if (nt <= 1 || n < 64) { worker(0, n); return; }
    int chunk = (n + (int)nt - 1) / (int)nt;
    std::vector<std::thread> ths;
    for (unsigned t = 0; t < nt; ++t) {
        int lo = (int)t * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        ths.emplace_back([lo, hi, &worker] { worker(lo, hi); });
    }
    for (auto& th : ths) th.join();
}

std::pair<int, int> edgeTypeDir(int a, int b) {
    if ((a == 0 && b == 1) || (a == 1 && b == 0)) return {0, (a > b) ? 1 : 0};
    if ((a == 1 && b == 2) || (a == 2 && b == 1)) return {1, (a > b) ? 1 : 0};
    return {2, (a == 0 && b == 2) ? 1 : 0};
}

struct EdgeOnElem {
    int et, dir;
    Vector2d nout;
    double he;
};

inline EdgeOnElem edgeOnElem(const Mesh& mesh, int t, int n1, int n2) {
    int a = -1, b = -1, w = -1;
    for (int k = 0; k < 3; ++k) {
        int v = mesh.elem(t, k);
        if (v == n1) a = k; else if (v == n2) b = k; else w = k;
    }
    Vector2d p1 = mesh.node.row(n1), p2 = mesh.node.row(n2), pw = mesh.node.row(mesh.elem(t, w));
    Vector2d evec = p2 - p1;
    double he = evec.norm();
    Vector2d nrm(evec.y() / he, -evec.x() / he);
    if (nrm.dot(pw - 0.5 * (p1 + p2)) > 0) nrm = -nrm;
    auto td = edgeTypeDir(a, b);
    return {td.first, td.second, nrm, he};
}
} // namespace

struct EulerDG::Impl {
    MatrixXd quadL;
    VectorXd wv;
    std::vector<RowVectorXd> phiV;
    std::vector<MatrixXd>    dphiV;
    MatrixXd quad1d;
    VectorXd w1d;
    int nqe = 0;
    std::vector<std::vector<std::vector<RowVectorXd>>> ephi;
    std::vector<std::vector<std::vector<MatrixXd>>>    edphi;
    MatrixXd Pdrop;
    MatrixXd Mref;
    struct EG { bool interior; int ta, tb, et_a, dir_a, et_b, dir_b, n1, n2; double nx, ny, he; };
    std::vector<EG> eg;
    struct EE { int nb, et, dir, et_nb, dir_nb, ei, n1, n2; double nx, ny, he; };
    std::vector<std::array<EE, 3>> ee;
};

} // namespace euler

#endif
