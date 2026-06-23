#include "Core.h"

#include "ElasticSolid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace euler_ale {

using namespace Eigen;

const int TAG_INTERIOR = 0;
const int TAG_MOVING_WALL = 1;
const int TAG_OUTFLOW = 2;
const int TAG_SLIP_WALL = 3;
const int TAG_EXACT = 9;

namespace {
template <class Worker>
void parallelRanges(int n, const Worker& worker) {
    unsigned hw = std::thread::hardware_concurrency();
    unsigned nt = std::max(1u, std::min<unsigned>(hw ? hw : 1u, static_cast<unsigned>(std::max(1, n))));
    if (nt <= 1 || n < 128) {
        worker(0, n);
        return;
    }
    int chunk = (n + static_cast<int>(nt) - 1) / static_cast<int>(nt);
    std::vector<std::thread> threads;
    threads.reserve(nt);
    for (unsigned i = 0; i < nt; ++i) {
        int lo = static_cast<int>(i) * chunk;
        int hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        threads.emplace_back([lo, hi, &worker]() { worker(lo, hi); });
    }
    for (auto& thread : threads) thread.join();
}
} // namespace

double triArea(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0));
    Vector2d p1 = m.node.row(m.elem(t, 1));
    Vector2d p2 = m.node.row(m.elem(t, 2));
    return 0.5 * std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
                          (p2.x() - p0.x()) * (p1.y() - p0.y()));
}

double longestEdge(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0));
    Vector2d p1 = m.node.row(m.elem(t, 1));
    Vector2d p2 = m.node.row(m.elem(t, 2));
    return std::max({(p1 - p0).norm(), (p2 - p1).norm(), (p0 - p2).norm()});
}

double hCFL(const Mesh& m, int t) {
    return 2.0 * triArea(m, t) / std::max(longestEdge(m, t), 1e-300);
}

std::pair<int, int> edgeTypeDir(int a, int b) {
    if ((a == 0 && b == 1) || (a == 1 && b == 0)) return {0, (a > b) ? 1 : 0};
    if ((a == 1 && b == 2) || (a == 2 && b == 1)) return {1, (a > b) ? 1 : 0};
    return {2, (a == 0 && b == 2) ? 1 : 0};
}

EdgeOnElem edgeOnElem(const Mesh& mesh, int t, int n1, int n2) {
    int a = -1, b = -1, w = -1;
    for (int k = 0; k < 3; ++k) {
        int v = mesh.elem(t, k);
        if (v == n1) a = k;
        else if (v == n2) b = k;
        else w = k;
    }
    Vector2d p1 = mesh.node.row(n1);
    Vector2d p2 = mesh.node.row(n2);
    Vector2d pw = mesh.node.row(mesh.elem(t, w));
    Vector2d evec = p2 - p1;
    double he = evec.norm();
    Vector2d nrm(evec.y() / he, -evec.x() / he);
    if (nrm.dot(pw - 0.5 * (p1 + p2)) > 0.0) nrm = -nrm;
    auto td = edgeTypeDir(a, b);
    return {td.first, td.second, nrm, he};
}


// ---------------------------------------------------------------------------
// Variable-component NVB forest implementation.
// ---------------------------------------------------------------------------
struct ALEAdaptiveForest::Impl {
public:
    Impl(const Mesh& baseMesh, int order, int nComp)
        : ord_(order), locDof_((order + 1) * (order + 2) / 2), nComp_(nComp) {
        int nB = static_cast<int>(baseMesh.node.rows());
        int NT = static_cast<int>(baseMesh.elem.rows());
        nBaseVerts_ = nB;
        vertsRef_.reserve(nB * 4);
        for (int i = 0; i < nB; ++i) vertsRef_.push_back(baseMesh.node.row(i).transpose());

        tris_.reserve(NT * 8);
        coef_.reserve(NT * 8);
        for (int t = 0; t < NT; ++t) {
            int g[3] = {baseMesh.elem(t, 0), baseMesh.elem(t, 1), baseMesh.elem(t, 2)};
            double L[3];
            for (int e = 0; e < 3; ++e)
                L[e] = (vertsRef_[g[e]] - vertsRef_[g[(e + 1) % 3]]).norm();
            int le = 0;
            if (L[1] > L[le]) le = 1;
            if (L[2] > L[le]) le = 2;
            Tri tr;
            tr.v = {g[(le + 2) % 3], g[le], g[(le + 1) % 3]};
            tris_.push_back(tr);
            coef_.emplace_back(MatrixXd::Zero(locDof_, nComp_));
        }
        buildTransferOps();
        leafOrder_.resize(NT);
        for (int t = 0; t < NT; ++t) leafOrder_[t] = t;
    }

    int numLeaves() const { return static_cast<int>(leafOrder_.size()); }
    int gen(int flatElem) const { return tris_[leafOrder_[flatElem]].gen; }

    void syncFromState(const MatrixXd& U, const MatrixXi& e2d, int locDof) {
        for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
            int t = leafOrder_[i];
            for (int k = 0; k < locDof; ++k) coef_[t].row(k) = U.row(e2d(i, k));
        }
    }

    MatrixXd gatherState(const MatrixXi& e2d, int locDof, int nDof) const {
        MatrixXd U = MatrixXd::Zero(nDof, nComp_);
        for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
            int t = leafOrder_[i];
            for (int k = 0; k < locDof; ++k) U.row(e2d(i, k)) = coef_[t].row(k);
        }
        return U;
    }

    void buildMesh(Mesh& mesh, const std::function<Vector2d(const Vector2d&)>& map) {
        leafOrder_.clear();
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t)
            if (tris_[t].leaf()) leafOrder_.push_back(t);

        std::unordered_map<int, int> remap;
        remap.reserve(leafOrder_.size() * 2);
        int nN = 0;
        for (int t : leafOrder_)
            for (int k = 0; k < 3; ++k) {
                int v = tris_[t].v[k];
                if (remap.find(v) == remap.end()) remap.emplace(v, nN++);
            }

        mesh.node.resize(nN, 2);
        for (const auto& kv : remap) mesh.node.row(kv.second) = map(vertsRef_[kv.first]).transpose();
        mesh.elem.resize(static_cast<int>(leafOrder_.size()), 3);
        for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
            const auto& v = tris_[leafOrder_[i]].v;
            mesh.elem(i, 0) = remap[v[0]];
            mesh.elem(i, 1) = remap[v[1]];
            mesh.elem(i, 2) = remap[v[2]];
        }
    }

    std::pair<int, int> adapt(const std::vector<int>& flag, int maxGen) {
        std::unordered_set<int> refineSet, coarsenSet;
        for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
            if (flag[i] > 0) refineSet.insert(leafOrder_[i]);
            else if (flag[i] < 0) coarsenSet.insert(leafOrder_[i]);
        }

        int nCoarsen = 0;
        {
            std::map<int, std::vector<int>> vert2leaf;
            for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
                if (!tris_[t].leaf()) continue;
                for (int k = 0; k < 3; ++k) vert2leaf[tris_[t].v[k]].push_back(t);
            }
            std::set<int> parentsToMerge;
            for (auto& kv : vert2leaf) {
                int m = kv.first;
                if (m < nBaseVerts_) continue;
                bool ok = true;
                for (int L : kv.second) {
                    if (tris_[L].v[0] != m || !coarsenSet.count(L) || tris_[L].parent < 0) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;
                for (int L : kv.second) parentsToMerge.insert(tris_[L].parent);
            }
            for (int p : parentsToMerge) {
                int c0 = tris_[p].child[0], c1 = tris_[p].child[1];
                if (c0 < 0 || c1 < 0 || !tris_[c0].leaf() || !tris_[c1].leaf()) continue;
                int wc0 = tris_[c0].whichChild, wc1 = tris_[c1].whichChild;
                coef_[p] = 0.5 * (Tcoarsen_[wc0] * coef_[c0] + Tcoarsen_[wc1] * coef_[c1]);
                tris_[c0].alive = false;
                tris_[c1].alive = false;
                tris_[p].child = {-1, -1};
                ++nCoarsen;
            }
        }

        int nRefine = 0;
        {
            std::map<std::pair<int, int>, std::array<int, 2>> e2l;
            rebuildLeafEdges(e2l);
            std::vector<int> todo;
            for (int t : refineSet)
                if (tris_[t].leaf() && tris_[t].gen < maxGen) todo.push_back(t);
            std::sort(todo.begin(), todo.end());
            for (int t : todo) {
                if (!tris_[t].leaf()) continue;
                refineElement(t, maxGen, e2l);
                ++nRefine;
            }
        }
        return {nRefine, nCoarsen};
    }

private:
    struct Tri {
        std::array<int, 3> v;
        std::array<int, 2> child{-1, -1};
        int parent = -1;
        int whichChild = -1;
        int gen = 0;
        bool alive = true;
        bool leaf() const { return alive && child[0] < 0; }
    };

    int ord_ = 0, locDof_ = 0, nComp_ = 0, nBaseVerts_ = 0;
    std::vector<Vector2d> vertsRef_;
    std::vector<Tri> tris_;
    std::vector<MatrixXd> coef_;
    std::map<std::pair<int, int>, int> midmap_;
    std::vector<int> leafOrder_;
    Matrix3d B_[2];
    MatrixXd Trestrict_[2], Tcoarsen_[2];

    static std::pair<int, int> key(int a, int b) {
        return {std::min(a, b), std::max(a, b)};
    }

    void buildTransferOps() {
        B_[0] << 0, 1, 0,
                 0.5, 0, 1,
                 0.5, 0, 0;
        B_[1] << 0, 0, 1,
                 0.5, 0, 0,
                 0.5, 1, 0;

        Mesh dummy;
        dummy.node.resize(3, 2);
        dummy.node << 0, 0, 1, 0, 0, 1;
        dummy.elem.resize(1, 3);
        dummy.elem << 0, 1, 2;
        FEM fem(ord_, dummy, false);
        MatrixXd quadL;
        VectorXd w;
        fem.quad2d(quadL, w);
        MatrixXd phiSelf = fem.computeBasisValue_all(quadL);
        MatrixXd MrefInv = (phiSelf.transpose() * w.asDiagonal() * phiSelf).inverse();
        for (int c = 0; c < 2; ++c) {
            MatrixXd mapped(quadL.rows(), 3);
            for (int q = 0; q < quadL.rows(); ++q)
                mapped.row(q) = (B_[c] * quadL.row(q).transpose()).transpose();
            MatrixXd phiPar = fem.computeBasisValue_all(mapped);
            Trestrict_[c] = MrefInv * (phiSelf.transpose() * w.asDiagonal() * phiPar);
            Tcoarsen_[c] = MrefInv * (phiPar.transpose() * w.asDiagonal() * phiSelf);
        }
    }

    int midpoint(int a, int b) {
        auto k = key(a, b);
        auto it = midmap_.find(k);
        if (it != midmap_.end()) return it->second;
        int m = static_cast<int>(vertsRef_.size());
        vertsRef_.push_back(0.5 * (vertsRef_[a] + vertsRef_[b]));
        midmap_[k] = m;
        return m;
    }

    static void edgeAdd(std::map<std::pair<int, int>, std::array<int, 2>>& e2l,
                        std::pair<int, int> k, int leaf) {
        auto it = e2l.find(k);
        if (it == e2l.end()) it = e2l.emplace(k, std::array<int, 2>{-1, -1}).first;
        auto& s = it->second;
        if (s[0] == leaf || s[1] == leaf) return;
        if (s[0] < 0) s[0] = leaf;
        else if (s[1] < 0) s[1] = leaf;
        else throw std::runtime_error("ALEAdaptiveForest: non-manifold edge");
    }

    static void edgeRemove(std::map<std::pair<int, int>, std::array<int, 2>>& e2l,
                           std::pair<int, int> k, int leaf) {
        auto it = e2l.find(k);
        if (it == e2l.end()) return;
        if (it->second[0] == leaf) it->second[0] = -1;
        else if (it->second[1] == leaf) it->second[1] = -1;
        if (it->second[0] < 0 && it->second[1] < 0) e2l.erase(it);
    }

    void rebuildLeafEdges(std::map<std::pair<int, int>, std::array<int, 2>>& e2l) const {
        e2l.clear();
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
            if (!tris_[t].leaf()) continue;
            const auto& v = tris_[t].v;
            edgeAdd(e2l, key(v[0], v[1]), t);
            edgeAdd(e2l, key(v[1], v[2]), t);
            edgeAdd(e2l, key(v[2], v[0]), t);
        }
    }

    int neighborAcross(int t, int a, int b,
                       const std::map<std::pair<int, int>, std::array<int, 2>>& e2l) const {
        auto it = e2l.find(key(a, b));
        if (it == e2l.end()) return -1;
        if (it->second[0] == t) return it->second[1];
        if (it->second[1] == t) return it->second[0];
        return -1;
    }

    void refineElement(int t, int maxGen,
                       std::map<std::pair<int, int>, std::array<int, 2>>& e2l) {
        static thread_local int depth = 0;
        struct Guard {
            int& d;
            explicit Guard(int& x) : d(x) { ++d; }
            ~Guard() { --d; }
        } guard(depth);
        if (depth > 4096) throw std::runtime_error("ALEAdaptiveForest: NVB recursion too deep");

        if (!tris_[t].leaf()) return;
        int a = -1, b = -1, F = -1;
        while (true) {
            if (!tris_[t].leaf()) return;
            a = tris_[t].v[1];
            b = tris_[t].v[2];
            F = neighborAcross(t, a, b, e2l);
            if (F < 0) break;
            auto re = std::array<int, 2>{tris_[F].v[1], tris_[F].v[2]};
            if (key(re[0], re[1]) == key(a, b)) break;
            refineElement(F, INT_MAX, e2l);
        }
        if (!tris_[t].leaf()) return;
        int m = midpoint(a, b);
        auto split = [&](int s) {
            int p = tris_[s].v[0], va = tris_[s].v[1], vb = tris_[s].v[2];
            int c0 = static_cast<int>(tris_.size());
            Tri t0;
            t0.v = {m, p, va};
            t0.parent = s;
            t0.whichChild = 0;
            t0.gen = tris_[s].gen + 1;
            tris_.push_back(t0);
            coef_.push_back(Trestrict_[0] * coef_[s]);
            int c1 = static_cast<int>(tris_.size());
            Tri t1;
            t1.v = {m, vb, p};
            t1.parent = s;
            t1.whichChild = 1;
            t1.gen = tris_[s].gen + 1;
            tris_.push_back(t1);
            coef_.push_back(Trestrict_[1] * coef_[s]);
            tris_[s].child = {c0, c1};
            edgeRemove(e2l, key(p, va), s);
            edgeRemove(e2l, key(va, vb), s);
            edgeRemove(e2l, key(vb, p), s);
            edgeAdd(e2l, key(m, p), c0);
            edgeAdd(e2l, key(p, va), c0);
            edgeAdd(e2l, key(va, m), c0);
            edgeAdd(e2l, key(m, vb), c1);
            edgeAdd(e2l, key(vb, p), c1);
            edgeAdd(e2l, key(p, m), c1);
        };
        split(t);
        if (F >= 0 && tris_[F].leaf()) split(F);
        (void)maxGen;
    }
};


ALEAdaptiveForest::ALEAdaptiveForest(const Mesh& baseMesh, int order, int nComp)
    : impl_(std::make_unique<Impl>(baseMesh, order, nComp)) {}

ALEAdaptiveForest::~ALEAdaptiveForest() = default;
ALEAdaptiveForest::ALEAdaptiveForest(ALEAdaptiveForest&&) noexcept = default;
ALEAdaptiveForest& ALEAdaptiveForest::operator=(ALEAdaptiveForest&&) noexcept = default;

int ALEAdaptiveForest::numLeaves() const { return impl_->numLeaves(); }
int ALEAdaptiveForest::gen(int flatElem) const { return impl_->gen(flatElem); }
void ALEAdaptiveForest::syncFromState(const MatrixXd& U, const MatrixXi& e2d, int locDof) {
    impl_->syncFromState(U, e2d, locDof);
}
MatrixXd ALEAdaptiveForest::gatherState(const MatrixXi& e2d, int locDof, int nDof) const {
    return impl_->gatherState(e2d, locDof, nDof);
}
void ALEAdaptiveForest::buildMesh(Mesh& mesh, const RefMapFn& map, double time) {
    impl_->buildMesh(mesh, [&](const Vector2d& X) { return map(X, time); });
}
std::pair<int, int> ALEAdaptiveForest::adapt(const std::vector<int>& flag, int maxGen) {
    return impl_->adapt(flag, maxGen);
}

Vector4d aleRusanov(const Vector4d& UL, const Vector4d& UR, double nx, double ny, double wn) {
    Vector4d FL = euler::normalFlux(UL, nx, ny) - wn * UL;
    Vector4d FR = euler::normalFlux(UR, nx, ny) - wn * UR;
    auto relSpeed = [&](const Vector4d& U) {
        double rho = std::max(U(0), 1e-300);
        double un = (U(1) * nx + U(2) * ny) / rho;
        return std::abs(un - wn) + euler::soundSpeed(U);
    };
    double a = std::max(relSpeed(UL), relSpeed(UR));
    return 0.5 * (FL + FR) - 0.5 * a * (UR - UL);
}

Vector4d movingWallGhost(const Vector4d& Um, double nx, double ny, double wn) {
    double rho = std::max(Um(0), 1e-14);
    double u = Um(1) / rho;
    double v = Um(2) / rho;
    double p = euler::pressure(Um);
    double un = u * nx + v * ny;
    double ug = u + 2.0 * (wn - un) * nx;
    double vg = v + 2.0 * (wn - un) * ny;
    return euler::primToCons(rho, ug, vg, std::max(p, 1e-12));
}

Vector4d characteristicPressureOutletGhost(const Vector4d& Um, double nx, double ny,
                                           double wn, const Vector4d& farPrimitive) {
    const double g = euler::GAMMA;
    const double gm1 = g - 1.0;

    double rho = std::max(Um(0), 1e-12);
    double u = Um(1) / rho;
    double v = Um(2) / rho;
    double p = std::max(euler::pressure(Um), 1e-12);
    double un = u * nx + v * ny;
    double relUn = un - wn;
    double c = std::sqrt(g * p / rho);
    if (relUn >= c) return Um;

    double rhoInf = std::max(farPrimitive(0), 1e-12);
    double uInf = farPrimitive(1);
    double vInf = farPrimitive(2);
    double pInf = std::max(farPrimitive(3), 1e-12);
    if (relUn <= 0.0) {
        return euler::primToCons(rhoInf, uInf, vInf, pInf);
    }

    double unInf = uInf * nx + vInf * ny;
    double cInf = std::sqrt(g * pInf / rhoInf);
    double rp = relUn + 2.0 * c / gm1;
    double rm = (unInf - wn) - 2.0 * cInf / gm1;
    double relUnB = 0.5 * (rp + rm);
    double cB = 0.25 * gm1 * (rp - rm);
    if (!std::isfinite(cB) || cB <= 1e-8) {
        return euler::primToCons(rhoInf, uInf, vInf, pInf);
    }

    double entropy = p / std::pow(rho, g);
    double rhoB = std::pow(std::max(cB * cB, 1e-16) / (g * entropy), 1.0 / gm1);
    double pB = rhoB * cB * cB / g;
    double tx = -ny;
    double ty = nx;
    double ut = u * tx + v * ty;
    double unB = relUnB + wn;
    double uB = unB * nx + ut * tx;
    double vB = unB * ny + ut * ty;
    return euler::primToCons(std::max(rhoB, 1e-12), uB, vB, std::max(pB, 1e-12));
}

struct ALEReferenceCache {
    int order = 0;
    int locDof = 0;
    MatrixXd quadL;
    MatrixXd quad1d;
    MatrixXd Mref;
    MatrixXd MrefInv;
    VectorXd wv;
    VectorXd w1d;
    int nqe = 0;
    std::vector<RowVectorXd> phiV;
    std::vector<MatrixXd> dphiV;
    std::vector<std::vector<std::vector<RowVectorXd>>> ephi;
};

const ALEReferenceCache& referenceCache(int order) {
    static std::map<int, std::unique_ptr<ALEReferenceCache>> caches;
    auto it = caches.find(order);
    if (it != caches.end()) return *it->second;

    auto cache = std::make_unique<ALEReferenceCache>();
    cache->order = order;
    Mesh dummy;
    dummy.node.resize(3, 2);
    dummy.node << 0.0, 0.0,
                  1.0, 0.0,
                  0.0, 1.0;
    dummy.elem.resize(1, 3);
    dummy.elem << 0, 1, 2;
    FEM fem(order, dummy, false);
    cache->locDof = fem.locDof;

    fem.quad2d(cache->quadL, cache->wv);
    int nqv = static_cast<int>(cache->wv.size());
    cache->phiV.resize(nqv);
    cache->dphiV.resize(nqv);
    cache->Mref = MatrixXd::Zero(cache->locDof, cache->locDof);
    for (int q = 0; q < nqv; ++q) {
        cache->phiV[q] = fem.computeBasisValue_all(cache->quadL.row(q)).row(0);
        cache->dphiV[q] = fem.computeBasisDlam_all(cache->quadL.row(q).transpose());
        cache->Mref.noalias() += cache->wv(q) *
                                  (cache->phiV[q].transpose() * cache->phiV[q]);
    }
    cache->MrefInv = cache->Mref.inverse();

    fem.quad1d(cache->quad1d, cache->w1d);
    cache->nqe = static_cast<int>(cache->w1d.size());
    cache->ephi.assign(3, std::vector<std::vector<RowVectorXd>>(
                              2, std::vector<RowVectorXd>(cache->nqe)));
    for (int et = 0; et < 3; ++et) {
        int i1 = 0, i2 = 1;
        if (et == 1) {
            i1 = 1;
            i2 = 2;
        } else if (et == 2) {
            i1 = 2;
            i2 = 0;
        }
        for (int dir = 0; dir < 2; ++dir) {
            for (int q = 0; q < cache->nqe; ++q) {
                double l1 = cache->quad1d(q, 0), l2 = cache->quad1d(q, 1);
                Vector3d lam = Vector3d::Zero();
                if (dir == 0) {
                    lam(i1) = l1;
                    lam(i2) = l2;
                } else {
                    lam(i1) = l2;
                    lam(i2) = l1;
                }
                cache->ephi[et][dir][q] =
                    fem.computeBasisValue_all(lam.transpose()).row(0);
            }
        }
    }

    auto [inserted, _] = caches.emplace(order, std::move(cache));
    return *inserted->second;
}

class ALEEulerDG {
public:
    ALEEulerDG(FEM& fem, Mesh& mesh, const MatrixXi& e2d,
               const MatrixXi& edge, const MatrixXi& e2s, const VectorXi& tag)
        : fem_(fem), mesh_(mesh), e2d_(e2d), edge_(edge), e2s_(e2s), tag_(tag),
          NT_(static_cast<int>(mesh.elem.rows())), NE_(static_cast<int>(edge.rows())),
          locDof_(fem.locDof), nDof_(static_cast<int>(e2d.maxCoeff()) + 1),
          ref_(referenceCache(fem.ord)) {
        precompute();
    }

    void applyMass(const MatrixXd& U, MatrixXd& MU) const {
        MU.resize(nDof_, 4);
        parallelRanges(NT_, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Ue(locDof_, 4);
            for (int t = lo; t < hi; ++t) {
                for (int i = 0; i < locDof_; ++i) Ue.row(i) = U.row(e2d_(t, i));
                Ue = (fem_.area(t) * ref_.Mref) * Ue;
                for (int i = 0; i < locDof_; ++i) MU.row(e2d_(t, i)) = Ue.row(i);
            }
        });
    }

    void applyMassInverse(MatrixXd& R) const {
        parallelRanges(NT_, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Re(locDof_, 4);
            for (int t = lo; t < hi; ++t) {
                for (int i = 0; i < locDof_; ++i) Re.row(i) = R.row(e2d_(t, i));
                Re = (ref_.MrefInv / fem_.area(t)) * Re;
                for (int i = 0; i < locDof_; ++i) R.row(e2d_(t, i)) = Re.row(i);
            }
        });
    }

    void residual(const MatrixXd& U, double time, const MeshVelocityFn& meshVel,
                  const ALEBCFn& bc, MatrixXd& R) const {
        auto velocityAt = [&](double x, double y, int) {
            return meshVel(x, y, time);
        };
        residualImpl(U, time, bc, R, velocityAt);
    }

    void residual(const MatrixXd& U, double time, const SolidALEMap& meshMap,
                  const ALEBCFn& bc, MatrixXd& R,
                  std::vector<int>& volumeHints, std::vector<int>& edgeHints) const {
        const int nqv = static_cast<int>(ref_.wv.size());
        const size_t nVolHints = static_cast<size_t>(NT_) * static_cast<size_t>(nqv);
        const size_t nEdgeHints = static_cast<size_t>(NT_) * 3u * static_cast<size_t>(ref_.nqe);
        if (volumeHints.size() != nVolHints) volumeHints.assign(nVolHints, -1);
        if (edgeHints.size() != nEdgeHints) edgeHints.assign(nEdgeHints, -1);
        auto velocityAt = [&](double x, double y, int hintSlot) {
            if (hintSlot >= 0 && hintSlot < static_cast<int>(volumeHints.size()))
                return meshMap.velocityAtCached(x, y, time, volumeHints[hintSlot]);
            int edgeSlot = hintSlot - static_cast<int>(volumeHints.size());
            return meshMap.velocityAtCached(x, y, time, edgeHints[edgeSlot]);
        };
        residualImpl(U, time, bc, R, velocityAt);
    }

private:
    template <class VelocityAt>
    void residualImpl(const MatrixXd& U, double time, const ALEBCFn& bc,
                      MatrixXd& R, VelocityAt&& velocityAt) const {
        R.setZero(nDof_, 4);
        const int nqv = static_cast<int>(ref_.wv.size());
        parallelRanges(NT_, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Ue(locDof_, 4), Unb(locDof_, 4), Re(locDof_, 4);
            MatrixXd G(2, locDof_);
            for (int tt = lo; tt < hi; ++tt) {
                for (int i = 0; i < locDof_; ++i) Ue.row(i) = U.row(e2d_(tt, i));
                Re.setZero();
                Vector2d p0 = mesh_.node.row(mesh_.elem(tt, 0));
                Vector2d p1 = mesh_.node.row(mesh_.elem(tt, 1));
                Vector2d p2 = mesh_.node.row(mesh_.elem(tt, 2));
                double area = fem_.area(tt);
                for (int q = 0; q < nqv; ++q) {
                    Vector4d Uq = (ref_.phiV[q] * Ue).transpose();
                    Vector4d Fx, Fy;
                    euler::fluxes(Uq, Fx, Fy);
                    Vector3d lam = ref_.quadL.row(q).transpose();
                    Vector2d xq = lam(0) * p0 + lam(1) * p1 + lam(2) * p2;
                    int hintSlot = tt * nqv + q;
                    Vector2d wq = velocityAt(xq.x(), xq.y(), hintSlot);
                    Fx.noalias() -= wq.x() * Uq;
                    Fy.noalias() -= wq.y() * Uq;
                    G.noalias() = fem_.Dlam[tt] * ref_.dphiV[q];
                    double wa = ref_.wv(q) * area;
                    Re.noalias() += (wa * G.row(0).transpose()) * Fx.transpose();
                    Re.noalias() += (wa * G.row(1).transpose()) * Fy.transpose();
                }

                for (int k = 0; k < 3; ++k) {
                    const EE& r = ee_[tt][k];
                    bool interior = (r.nb != -1);
                    if (interior)
                        for (int i = 0; i < locDof_; ++i) Unb.row(i) = U.row(e2d_(r.nb, i));
                    int tag = tag_.size() ? tag_(r.ei) : 0;
                    for (int q = 0; q < ref_.nqe; ++q) {
                        const RowVectorXd& pa = ref_.ephi[r.et][r.dir][q];
                        Vector4d Um = (pa * Ue).transpose();
                        double l1 = ref_.quad1d(q, 0), l2 = ref_.quad1d(q, 1);
                        Vector2d Pp = l1 * mesh_.node.row(r.n1).transpose() +
                                      l2 * mesh_.node.row(r.n2).transpose();
                        int hintSlot = NT_ * nqv + (tt * 3 + k) * ref_.nqe + q;
                        Vector2d w = velocityAt(Pp.x(), Pp.y(), hintSlot);
                        double wn = w.x() * r.nx + w.y() * r.ny;
                        Vector4d Up;
                        if (interior) {
                            Up = (ref_.ephi[r.et_nb][r.dir_nb][q] * Unb).transpose();
                        } else {
                            Up = bc ? bc(Pp.x(), Pp.y(), time, Um, r.nx, r.ny, tag, wn) : Um;
                        }
                        Vector4d Hn = aleRusanov(Um, Up, r.nx, r.ny, wn);
                        double whe = ref_.w1d(q) * r.he;
                        Re.noalias() -= (whe * pa.transpose()) * Hn.transpose();
                    }
                }
                for (int i = 0; i < locDof_; ++i) R.row(e2d_(tt, i)) = Re.row(i);
            }
        });
    }

    struct EE {
        int nb = -1, et = 0, dir = 0, et_nb = 0, dir_nb = 0, ei = -1, n1 = -1, n2 = -1;
        double nx = 0.0, ny = 0.0, he = 0.0;
    };

    FEM& fem_;
    Mesh& mesh_;
    const MatrixXi& e2d_;
    const MatrixXi& edge_;
    const MatrixXi& e2s_;
    VectorXi tag_;
    int NT_ = 0, NE_ = 0, locDof_ = 0, nDof_ = 0;
    const ALEReferenceCache& ref_;
    std::vector<std::array<EE, 3>> ee_;

    void precompute() {
        auto pkey = [](int a, int b) {
            return (static_cast<int64_t>(std::min(a, b)) << 32) |
                   static_cast<uint32_t>(std::max(a, b));
        };
        std::unordered_map<int64_t, int> emap;
        emap.reserve(static_cast<size_t>(NE_) * 2);
        for (int e = 0; e < NE_; ++e) emap[pkey(edge_(e, 0), edge_(e, 1))] = e;
        ee_.assign(NT_, std::array<EE, 3>{});
        const int loc[3][2] = {{0, 1}, {1, 2}, {2, 0}};
        for (int t = 0; t < NT_; ++t) {
            for (int k = 0; k < 3; ++k) {
                int va = mesh_.elem(t, loc[k][0]), vb = mesh_.elem(t, loc[k][1]);
                int e = emap[pkey(va, vb)];
                int n1 = edge_(e, 0), n2 = edge_(e, 1);
                EdgeOnElem es = edgeOnElem(mesh_, t, n1, n2);
                int nb = (e2s_(e, 0) == t) ? e2s_(e, 1) : e2s_(e, 0);
                EE r;
                r.nb = nb;
                r.et = es.et;
                r.dir = es.dir;
                r.ei = e;
                r.n1 = n1;
                r.n2 = n2;
                r.nx = es.nout.x();
                r.ny = es.nout.y();
                r.he = es.he;
                if (nb != -1) {
                    EdgeOnElem en = edgeOnElem(mesh_, nb, n1, n2);
                    r.et_nb = en.et;
                    r.dir_nb = en.dir;
                }
                ee_[t][k] = r;
            }
        }
    }
};


void rebuildSpace(ALEAdaptiveForest& forest, int ord, const RefMapFn& refMap,
                  double t, const Tagger& tagger, Space& sp) {
    forest.buildMesh(sp.mesh, refMap, t);
    sp.fem = std::make_unique<FEM>(ord, sp.mesh, false);
    sp.fem->getDOF(sp.mesh, sp.e2d, sp.nDof);
    sp.mesh.getEdge2Side(sp.edge, sp.e2s);
    sp.minH = 1e300;
    for (int k = 0; k < sp.mesh.elem.rows(); ++k) sp.minH = std::min(sp.minH, hCFL(sp.mesh, k));
    sp.tag = VectorXi::Zero(sp.edge.rows());
    sp.volumeVelocityHints.clear();
    sp.edgeVelocityHints.clear();
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.e2s(e, 0) != -1 && sp.e2s(e, 1) != -1) continue;
        Vector2d m = 0.5 * (sp.mesh.node.row(sp.edge(e, 0)) +
                            sp.mesh.node.row(sp.edge(e, 1))).transpose();
        sp.tag(e) = tagger(m.x(), m.y(), t);
    }
}

void updateStaticSpace(StaticSpace& ss, const RefMapFn& refMap, double t) {
    if (!ss.space.fem) throw std::runtime_error("updateStaticSpace: uninitialized StaticSpace");
    ss.space.mesh.node.resize(ss.referenceMesh.node.rows(), 2);
    for (int i = 0; i < ss.referenceMesh.node.rows(); ++i) {
        Vector2d X = ss.referenceMesh.node.row(i).transpose();
        ss.space.mesh.node.row(i) = refMap(X, t).transpose();
    }
    gradbasis_my(ss.space.mesh, ss.space.fem->Dlam, ss.space.fem->area);
    ss.space.minH = 1e300;
    for (int k = 0; k < ss.space.mesh.elem.rows(); ++k)
        ss.space.minH = std::min(ss.space.minH, hCFL(ss.space.mesh, k));
}

void initializeStaticSpace(ALEAdaptiveForest& forest, int ord, const RefMapFn& refMap,
                           double t, const Tagger& tagger, StaticSpace& ss) {
    RefMapFn identity = [](const Vector2d& X, double) { return X; };
    forest.buildMesh(ss.referenceMesh, identity, 0.0);
    ss.space.mesh = ss.referenceMesh;
    ss.space.fem = std::make_unique<FEM>(ord, ss.space.mesh, false);
    ss.space.fem->getDOF(ss.space.mesh, ss.space.e2d, ss.space.nDof);
    ss.space.mesh.getEdge2Side(ss.space.edge, ss.space.e2s);
    ss.space.volumeVelocityHints.clear();
    ss.space.edgeVelocityHints.clear();
    updateStaticSpace(ss, refMap, t);

    ss.space.tag = VectorXi::Zero(ss.space.edge.rows());
    for (int e = 0; e < ss.space.edge.rows(); ++e) {
        if (ss.space.e2s(e, 0) != -1 && ss.space.e2s(e, 1) != -1) continue;
        Vector2d m = 0.5 * (ss.space.mesh.node.row(ss.space.edge(e, 0)) +
                            ss.space.mesh.node.row(ss.space.edge(e, 1))).transpose();
        ss.space.tag(e) = tagger(m.x(), m.y(), t);
    }
}

MatrixXd advanceOne(ALEAdaptiveForest& forest, int ord, const RefMapFn& refMap,
                    double t, double dt, const Tagger& tagger,
                    const MeshVelocityFn& meshVel, const ALEBCFn& bc,
                    const MatrixXd& U) {
    Space sn, s1;
    rebuildSpace(forest, ord, refMap, t, tagger, sn);
    rebuildSpace(forest, ord, refMap, t + dt, tagger, s1);
    ALEEulerDG dgN(*sn.fem, sn.mesh, sn.e2d, sn.edge, sn.e2s, sn.tag);
    ALEEulerDG dg1(*s1.fem, s1.mesh, s1.e2d, s1.edge, s1.e2s, s1.tag);

    MatrixXd Qn, R1, R2;
    dgN.applyMass(U, Qn);
    dgN.residual(U, t, meshVel, bc, R1);
    MatrixXd Q1 = Qn + dt * R1;
    MatrixXd U1 = Q1;
    dg1.applyMassInverse(U1);
    dg1.residual(U1, t + dt, meshVel, bc, R2);
    MatrixXd Qnew = 0.5 * Qn + 0.5 * (Q1 + dt * R2);
    MatrixXd Unew = Qnew;
    dg1.applyMassInverse(Unew);
    return Unew;
}

MatrixXd advanceOneStatic(StaticSpace& ssN, StaticSpace& ss1,
                          const RefMapFn& refMap, double t, double dt,
                          const MeshVelocityFn& meshVel, const ALEBCFn& bc,
                          const MatrixXd& U) {
    updateStaticSpace(ssN, refMap, t);
    updateStaticSpace(ss1, refMap, t + dt);
    Space& sn = ssN.space;
    Space& s1 = ss1.space;
    ALEEulerDG dgN(*sn.fem, sn.mesh, sn.e2d, sn.edge, sn.e2s, sn.tag);
    ALEEulerDG dg1(*s1.fem, s1.mesh, s1.e2d, s1.edge, s1.e2s, s1.tag);

    MatrixXd Qn, R1, R2;
    dgN.applyMass(U, Qn);
    dgN.residual(U, t, meshVel, bc, R1);
    MatrixXd Q1 = Qn + dt * R1;
    MatrixXd U1 = Q1;
    dg1.applyMassInverse(U1);
    dg1.residual(U1, t + dt, meshVel, bc, R2);
    MatrixXd Qnew = 0.5 * Qn + 0.5 * (Q1 + dt * R2);
    MatrixXd Unew = Qnew;
    dg1.applyMassInverse(Unew);
    return Unew;
}

MatrixXd advanceOneStatic(StaticSpace& ssN, StaticSpace& ss1,
                          const RefMapFn& refMap, double t, double dt,
                          const SolidALEMap& meshMap, const ALEBCFn& bc,
                          const MatrixXd& U) {
    updateStaticSpace(ssN, refMap, t);
    updateStaticSpace(ss1, refMap, t + dt);
    Space& sn = ssN.space;
    Space& s1 = ss1.space;
    ALEEulerDG dgN(*sn.fem, sn.mesh, sn.e2d, sn.edge, sn.e2s, sn.tag);
    ALEEulerDG dg1(*s1.fem, s1.mesh, s1.e2d, s1.edge, s1.e2s, s1.tag);

    MatrixXd Qn, R1, R2;
    dgN.applyMass(U, Qn);
    dgN.residual(U, t, meshMap, bc, R1, sn.volumeVelocityHints, sn.edgeVelocityHints);
    MatrixXd Q1 = Qn + dt * R1;
    MatrixXd U1 = Q1;
    dg1.applyMassInverse(U1);
    dg1.residual(U1, t + dt, meshMap, bc, R2, s1.volumeVelocityHints, s1.edgeVelocityHints);
    MatrixXd Qnew = 0.5 * Qn + 0.5 * (Q1 + dt * R2);
    MatrixXd Unew = Qnew;
    dg1.applyMassInverse(Unew);
    return Unew;
}

double estimateDt(const Space& sp, const MatrixXd& U, const MaxMeshSpeedFn& maxMeshSpeed,
                  double t, int ord, double cfl) {
    double h = (sp.minH > 0.0) ? sp.minH : 1e300;
    if (!(sp.minH > 0.0))
        for (int k = 0; k < sp.mesh.elem.rows(); ++k) h = std::min(h, hCFL(sp.mesh, k));
    double lam = 1e-12;
    for (int i = 0; i < U.rows(); ++i) {
        Vector4d Ui = U.row(i).transpose();
        double rho = std::max(Ui(0), 1e-12);
        double spd = std::hypot(Ui(1), Ui(2)) / rho + euler::soundSpeed(Ui);
        lam = std::max(lam, spd);
    }
    lam += maxMeshSpeed(t);
    return cfl * h / ((2.0 * ord + 1.0) * lam);
}

std::vector<int> computeAMRFlags(const Space& sp, const MatrixXd& U, const ALEAdaptiveForest& forest,
                                 int maxGen, double thRef, double thCrs, int bufferLayers,
                                 bool allowCoarsen) {
    int NT = sp.mesh.elem.rows();
    int locDof = sp.fem->locDof;
    VectorXd rbar(NT), eta(NT);
    MatrixXd dphiC = sp.fem->computeBasisDlam_all(Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
    VectorXd rblk(locDof);
    for (int t = 0; t < NT; ++t) {
        double rb = 0.0;
        for (int i = 0; i < locDof; ++i) {
            rblk(i) = U(sp.e2d(t, i), 0);
            rb += rblk(i);
        }
        rb /= locDof;
        rbar(t) = rb;
        Vector2d g = sp.fem->Dlam[t] * (dphiC * rblk);
        eta(t) = longestEdge(sp.mesh, t) * g.norm() / std::max(rb, 1e-9);
    }
    for (int e = 0; e < sp.e2s.rows(); ++e) {
        int a = sp.e2s(e, 0), b = sp.e2s(e, 1);
        if (a < 0 || b < 0) continue;
        double j = std::abs(rbar(a) - rbar(b)) / std::max(std::min(rbar(a), rbar(b)), 1e-9);
        eta(a) = std::max(eta(a), j);
        eta(b) = std::max(eta(b), j);
    }
    std::vector<int> flag(NT, 0);
    for (int t = 0; t < NT; ++t) {
        if (eta(t) > thRef && forest.gen(t) < maxGen) flag[t] = 1;
        else if (allowCoarsen && eta(t) < thCrs && forest.gen(t) > 0) flag[t] = -1;
    }
    std::vector<std::vector<int>> nbr(NT);
    for (int e = 0; e < sp.e2s.rows(); ++e) {
        int a = sp.e2s(e, 0), b = sp.e2s(e, 1);
        if (a >= 0 && b >= 0) {
            nbr[a].push_back(b);
            nbr[b].push_back(a);
        }
    }
    std::vector<int> front;
    for (int t = 0; t < NT; ++t)
        if (flag[t] == 1) front.push_back(t);
    for (int layer = 0; layer < bufferLayers; ++layer) {
        std::vector<int> next;
        for (int t : front) {
            for (int n : nbr[t]) {
                if (flag[n] != 1 && forest.gen(n) < maxGen) {
                    flag[n] = 1;
                    next.push_back(n);
                }
            }
        }
        front.swap(next);
        if (front.empty()) break;
    }
    for (int t = 0; t < NT; ++t) {
        if (flag[t] != -1) continue;
        for (int n : nbr[t]) {
            if (flag[n] == 1) {
                flag[t] = 0;
                break;
            }
        }
    }
    return flag;
}

bool readPPM(const std::string& path, int& W, int& H, std::vector<unsigned char>& img) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string magic;
    int maxv = 0;
    in >> magic >> W >> H >> maxv;
    in.get();
    if (magic != "P6" || W <= 0 || H <= 0) return false;
    img.resize(static_cast<size_t>(W) * H * 3);
    in.read(reinterpret_cast<char*>(img.data()), static_cast<std::streamsize>(img.size()));
    return true;
}

void writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img) {
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
}

std::vector<unsigned char> downsampleImage(const std::vector<unsigned char>& src,
                                           int srcW, int srcH, int factor,
                                           int& dstW, int& dstH) {
    if (factor <= 1) {
        dstW = srcW;
        dstH = srcH;
        return src;
    }
    if (srcW <= 0 || srcH <= 0 || srcW % factor != 0 || srcH % factor != 0) {
        dstW = 0;
        dstH = 0;
        return {};
    }

    dstW = srcW / factor;
    dstH = srcH / factor;
    std::vector<unsigned char> dst(static_cast<size_t>(dstW) * dstH * 3, 255);
    int samples = factor * factor;
    for (int y = 0; y < dstH; ++y) {
        for (int x = 0; x < dstW; ++x) {
            unsigned int acc[3] = {0, 0, 0};
            for (int sy = 0; sy < factor; ++sy) {
                int srcY = y * factor + sy;
                for (int sx = 0; sx < factor; ++sx) {
                    int srcX = x * factor + sx;
                    size_t srcIdx = (static_cast<size_t>(srcY) * srcW + srcX) * 3;
                    acc[0] += src[srcIdx];
                    acc[1] += src[srcIdx + 1];
                    acc[2] += src[srcIdx + 2];
                }
            }
            size_t dstIdx = (static_cast<size_t>(y) * dstW + x) * 3;
            dst[dstIdx] = static_cast<unsigned char>((acc[0] + samples / 2) / samples);
            dst[dstIdx + 1] = static_cast<unsigned char>((acc[1] + samples / 2) / samples);
            dst[dstIdx + 2] = static_cast<unsigned char>((acc[2] + samples / 2) / samples);
        }
    }
    return dst;
}

bool downsamplePPM(const std::string& srcPath, const std::string& dstPath,
                   int factor) {
    int srcW = 0, srcH = 0;
    std::vector<unsigned char> src;
    if (!readPPM(srcPath, srcW, srcH, src)) return false;
    int dstW = 0, dstH = 0;
    std::vector<unsigned char> dst = downsampleImage(src, srcW, srcH, factor, dstW, dstH);
    if (dst.empty() || dstW <= 0 || dstH <= 0) return false;
    writePPM(dstPath, dstW, dstH, dst);
    return true;
}

void drawLine(std::vector<unsigned char>& img, int W, int H, int x0, int y0, int x1, int y1) {
    int dx = std::abs(x1 - x0), dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1, err = dx + dy;
    while (true) {
        if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) {
            size_t k = (static_cast<size_t>(y0) * W + x0) * 3;
            img[k] = static_cast<unsigned char>(0.55 * 30 + 0.45 * img[k]);
            img[k + 1] = static_cast<unsigned char>(0.55 * 220 + 0.45 * img[k + 1]);
            img[k + 2] = static_cast<unsigned char>(0.55 * 255 + 0.45 * img[k + 2]);
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

void overlayMesh(std::vector<unsigned char>& img, int W, int H,
                 const Mesh& mesh, double xa, double xb, double ya, double yb) {
    auto toPix = [&](const Vector2d& p, int& px, int& py) {
        px = static_cast<int>(std::lround((p.x() - xa) / (xb - xa) * W));
        py = static_cast<int>(std::lround((yb - p.y()) / (yb - ya) * H));
    };
    for (int t = 0; t < mesh.elem.rows(); ++t) {
        for (int k = 0; k < 3; ++k) {
            Vector2d A = mesh.node.row(mesh.elem(t, k));
            Vector2d B = mesh.node.row(mesh.elem(t, (k + 1) % 3));
            int ax, ay, bx, by;
            toPix(A, ax, ay);
            toPix(B, bx, by);
            drawLine(img, W, H, ax, ay, bx, by);
        }
    }
}

void overlayMesh(const std::string& path, const Mesh& mesh,
                 double xa, double xb, double ya, double yb) {
    int W = 0, H = 0;
    std::vector<unsigned char> img;
    if (!readPPM(path, W, H, img)) return;
    overlayMesh(img, W, H, mesh, xa, xb, ya, yb);
    writePPM(path, W, H, img);
}

} // namespace euler_ale
