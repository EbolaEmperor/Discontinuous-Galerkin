#include "AdaptiveMesh.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace euler {

// ===========================================================================
// Construction: detect the criss-cross base topology, label it for NVB, and
// precompute the constant solution-transfer operators.
// ===========================================================================
AdaptiveForest::AdaptiveForest(const Mesh& baseMesh, int order)
    : ord_(order), locDof_((order + 1) * (order + 2) / 2) {
    int nB = static_cast<int>(baseMesh.node.rows());
    int NT = static_cast<int>(baseMesh.elem.rows());
    nBaseVerts_ = nB;
    verts_.reserve(nB * 2);
    for (int i = 0; i < nB; ++i) verts_.push_back(baseMesh.node.row(i).transpose());

    // NVB labelling of each base triangle: refinement edge = LONGEST edge (for a
    // criss-cross right-triangle mesh this is the shared cell diagonal, so both
    // triangles of a cell pick the same shared edge -> the mesh is "compatibly
    // labelled" and the closure recursion terminates).  Store as (peak, e0, e1)
    // with peak = vertex opposite the refinement edge.
    tris_.reserve(NT * 4);
    coef_.reserve(NT * 4);
    for (int t = 0; t < NT; ++t) {
        int g[3] = {baseMesh.elem(t, 0), baseMesh.elem(t, 1), baseMesh.elem(t, 2)};
        double L[3];
        for (int e = 0; e < 3; ++e)
            L[e] = (verts_[g[e]] - verts_[g[(e + 1) % 3]]).norm();   // edge (e, e+1)
        int le = 0; if (L[1] > L[le]) le = 1; if (L[2] > L[le]) le = 2;  // longest edge index
        int a = g[le], b = g[(le + 1) % 3], peak = g[(le + 2) % 3];      // refedge (a,b), opposite peak
        Tri tr; tr.v = {peak, a, b}; tr.gen = 0; tr.parent = -1;
        tris_.push_back(tr);
        coef_.emplace_back(MatrixXd::Zero(locDof_, 4));
    }
    buildTransferOps();

    // initial flat mesh = the base leaves
    leafOrder_.resize(NT);
    for (int t = 0; t < NT; ++t) leafOrder_[t] = t;
}

// ---------------------------------------------------------------------------
// Constant transfer operators.  The child->parent barycentric maps B_c and the
// child:parent area ratio (=1/2 for a bisection) are fixed, so restriction and
// L2-projection are element-independent matrices.  child0=(m,v0,v1) uses B_0,
// child1=(m,v2,v0) uses B_1, with m the midpoint of the refinement edge (v1,v2).
// ---------------------------------------------------------------------------
void AdaptiveForest::buildTransferOps() {
    // parent barycentric coords of each child vertex (columns = child v0,v1,v2):
    //   v0=peak->(1,0,0)  v1->(0,1,0)  v2->(0,0,1)  m=mid(v1,v2)->(0,.5,.5)
    B_[0] << 0, 1, 0,                  // child0 = (m, v0, v1)
             0.5, 0, 1,
             0.5, 0, 0;
    B_[1] << 0, 0, 1,                  // child1 = (m, v2, v0)
             0.5, 0, 0,
             0.5, 1, 0;

    Mesh dummy; dummy.node.resize(3, 2); dummy.node << 0, 0, 1, 0, 0, 1;
    dummy.elem.resize(1, 3); dummy.elem << 0, 1, 2;
    FEM fem(ord_, dummy);
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);            // nq x 3, weights
    int nq = static_cast<int>(w.size());
    MatrixXd phiSelf = fem.computeBasisValue_all(quadL);        // nq x locDof
    MatrixXd Mref = phiSelf.transpose() * w.asDiagonal() * phiSelf;
    MatrixXd MrefInv = Mref.inverse();
    for (int c = 0; c < 2; ++c) {
        MatrixXd mapped(nq, 3);
        for (int q = 0; q < nq; ++q) mapped.row(q) = (B_[c] * quadL.row(q).transpose()).transpose();
        MatrixXd phiPar = fem.computeBasisValue_all(mapped);   // parent basis at the same physical pts
        // childCoeff = MrefInv (phiSelf^T W phiPar) parentCoeff  (exact restriction)
        Trestrict_[c] = MrefInv * (phiSelf.transpose() * w.asDiagonal() * phiPar);
        // parentLoad += MrefInv (phiPar^T W phiSelf) childCoeff   (L2 projection load)
        Tcoarsen_[c] = MrefInv * (phiPar.transpose() * w.asDiagonal() * phiSelf);
    }
}

int AdaptiveForest::midpoint(int a, int b) {
    auto k = key(a, b);
    auto it = midmap_.find(k);
    if (it != midmap_.end()) return it->second;
    int m = static_cast<int>(verts_.size());
    verts_.push_back(0.5 * (verts_[a] + verts_[b]));
    midmap_[k] = m;
    return m;
}

// ---------------------------------------------------------------------------
// Low-level: split leaf t along its refinement edge (v1,v2) at midpoint m,
// creating its two NVB children, transferring the solution EXACTLY, and keeping
// the leaf-edge adjacency map consistent.
// ---------------------------------------------------------------------------
namespace {
inline void edgeAdd(std::map<std::pair<int, int>, std::array<int, 2>>& e2l,
                    std::pair<int, int> k, int leaf) {
    // NB: map operator[] value-initializes a NEW std::array to {0,0}, which would
    // masquerade as "occupied by leaf 0"; insert with the {-1,-1} empty sentinel.
    auto it = e2l.find(k);
    if (it == e2l.end()) it = e2l.emplace(k, std::array<int, 2>{-1, -1}).first;
    auto& s = it->second;
    if (s[0] == leaf || s[1] == leaf) return;              // idempotent (defensive)
    if (s[0] < 0) s[0] = leaf;
    else if (s[1] < 0) s[1] = leaf;
    else {
        std::ostringstream os;
        os << "AdaptiveForest: >2 leaves on edge (" << k.first << "," << k.second
           << "): existing {" << s[0] << "," << s[1] << "} + new " << leaf;
        throw std::runtime_error(os.str());
    }
}
inline void edgeRemove(std::map<std::pair<int, int>, std::array<int, 2>>& e2l,
                       std::pair<int, int> k, int leaf) {
    auto it = e2l.find(k);
    if (it == e2l.end()) return;
    if (it->second[0] == leaf) it->second[0] = -1;
    else if (it->second[1] == leaf) it->second[1] = -1;
    if (it->second[0] < 0 && it->second[1] < 0) e2l.erase(it);
}
} // namespace

void AdaptiveForest::refineElement(int t, int maxGen,
                                   std::map<std::pair<int, int>, std::array<int, 2>>& e2l) {
    // closure depth guard: NVB on a compatibly-labelled mesh terminates quickly;
    // a runaway means a labelling/topology bug -> fail loudly rather than hang.
    // RAII so the counter is restored even if a deeper call throws.
    static thread_local int depth = 0;
    struct DepthGuard { int& d; DepthGuard(int& x) : d(x) { ++d; } ~DepthGuard() { --d; } } guard(depth);
    if (depth > 4096) throw std::runtime_error("AdaptiveForest: NVB closure recursion too deep");

    if (!tris_[t].leaf()) return;
    int a = -1, b = -1, F = -1;
    while (true) {
        // A neighbour's conformity closure can wrap back and bisect t itself (when
        // t's REFINEMENT happens via one of its legs); if so t is already refined
        // -- the outer call's intent is satisfied, so bail before double-splitting.
        if (!tris_[t].leaf()) return;
        a = tris_[t].v[1]; b = tris_[t].v[2];
        F = neighborAcross(t, a, b, e2l);
        if (F < 0) break;                                  // boundary refinement edge -> compatible
        auto re = refEdge(F);
        if (key(re[0], re[1]) == key(a, b)) break;         // neighbour shares it as ITS refedge
        refineElement(F, INT_MAX, e2l);                    // make F compatible first (closure)
    }
    if (!tris_[t].leaf()) return;                          // t was bisected during the closure above
    // (a,b) is now compatibly divisible: split t and (if interior) F at one midpoint.
    int m = midpoint(a, b);
    auto split = [&](int s) {
        int p = tris_[s].v[0], va = tris_[s].v[1], vb = tris_[s].v[2];
        int c0 = static_cast<int>(tris_.size());
        Tri t0; t0.v = {m, p, va}; t0.parent = s; t0.whichChild = 0; t0.gen = tris_[s].gen + 1;
        tris_.push_back(t0); coef_.push_back(Trestrict_[0] * coef_[s]);
        int c1 = static_cast<int>(tris_.size());
        Tri t1; t1.v = {m, vb, p}; t1.parent = s; t1.whichChild = 1; t1.gen = tris_[s].gen + 1;
        tris_.push_back(t1); coef_.push_back(Trestrict_[1] * coef_[s]);
        tris_[s].child = {c0, c1};
        // adjacency: parent off its 3 edges, children on their 3 edges
        edgeRemove(e2l, key(p, va), s); edgeRemove(e2l, key(va, vb), s); edgeRemove(e2l, key(vb, p), s);
        edgeAdd(e2l, key(m, p), c0); edgeAdd(e2l, key(p, va), c0); edgeAdd(e2l, key(va, m), c0);
        edgeAdd(e2l, key(m, vb), c1); edgeAdd(e2l, key(vb, p), c1); edgeAdd(e2l, key(p, m), c1);
    };
    split(t);
    if (F >= 0 && tris_[F].leaf()) split(F);               // F stays a leaf through t's split
}

int AdaptiveForest::neighborAcross(int t, int a, int b,
                                   const std::map<std::pair<int, int>, std::array<int, 2>>& e2l) const {
    auto it = e2l.find(key(a, b));
    if (it == e2l.end()) return -1;
    if (it->second[0] == t) return it->second[1];
    if (it->second[1] == t) return it->second[0];
    return -1;
}

void AdaptiveForest::rebuildLeafEdges(std::map<std::pair<int, int>, std::array<int, 2>>& e2l) const {
    e2l.clear();
    for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
        if (!tris_[t].leaf()) continue;
        const auto& v = tris_[t].v;
        edgeAdd(e2l, key(v[0], v[1]), t);
        edgeAdd(e2l, key(v[1], v[2]), t);
        edgeAdd(e2l, key(v[2], v[0]), t);
    }
}

// ===========================================================================
// State sync / gather (couple the forest leaf coefficients to the DG state).
// ===========================================================================
void AdaptiveForest::syncFromState(const MatrixXd& U, const MatrixXi& e2d, int locDof) {
    for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
        int t = leafOrder_[i];
        for (int k = 0; k < locDof; ++k) coef_[t].row(k) = U.row(e2d(i, k));
    }
}

MatrixXd AdaptiveForest::gatherState(const MatrixXi& e2d, int locDof, int nDof) const {
    MatrixXd U(nDof, 4);
    for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
        int t = leafOrder_[i];
        for (int k = 0; k < locDof; ++k) U.row(e2d(i, k)) = coef_[t].row(k);
    }
    return U;
}

// ===========================================================================
// One adaptation pass: coarsen (one level) then refine (one level).
// flag[flatElem]: +1 refine, -1 coarsen, 0 hold.
// ===========================================================================
std::pair<int, int> AdaptiveForest::adapt(const std::vector<int>& flag, int maxGen) {
    // Map the flat-element flags onto stable forest tri ids BEFORE any mutation.
    std::unordered_set<int> refineSet, coarsenSet;
    for (int i = 0; i < static_cast<int>(leafOrder_.size()); ++i) {
        if (flag[i] > 0) refineSet.insert(leafOrder_[i]);
        else if (flag[i] < 0) coarsenSet.insert(leafOrder_[i]);
    }

    // ---------------- COARSEN (one level) ----------------
    // A midpoint m is removable iff every current leaf incident to m has m as its
    // peak (nothing finer hangs off m) AND all those leaves are flagged coarsen.
    // Then un-bisect every parent pair around m together (conformity preserved).
    int nCoarsen = 0;
    {
        std::map<int, std::vector<int>> vert2leaf;
        for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
            if (!tris_[t].leaf()) continue;
            for (int k = 0; k < 3; ++k) vert2leaf[tris_[t].v[k]].push_back(t);
        }
        std::set<int> parentsToMerge;
        for (auto& [m, leaves] : vert2leaf) {
            if (m < nBaseVerts_) continue;                 // base vertices are never removable
            bool ok = true;
            for (int L : leaves)
                if (tris_[L].v[0] != m || !coarsenSet.count(L) || tris_[L].parent < 0) { ok = false; break; }
            if (!ok) continue;
            for (int L : leaves) parentsToMerge.insert(tris_[L].parent);
        }
        for (int p : parentsToMerge) {
            int c0 = tris_[p].child[0], c1 = tris_[p].child[1];
            if (c0 < 0 || c1 < 0 || !tris_[c0].leaf() || !tris_[c1].leaf()) continue;  // defensive
            // whichChild tells us which transfer matrix each child uses
            int wc0 = tris_[c0].whichChild, wc1 = tris_[c1].whichChild;
            MatrixXd pc = 0.5 * (Tcoarsen_[wc0] * coef_[c0] + Tcoarsen_[wc1] * coef_[c1]);
            coef_[p] = pc;
            tris_[c0].alive = false; tris_[c1].alive = false;  // orphan the old children
            tris_[p].child = {-1, -1};                      // parent becomes a leaf again
            ++nCoarsen;
        }
    }

    // ---------------- REFINE (one level) ----------------
    int nRefine = 0;
    {
        std::map<std::pair<int, int>, std::array<int, 2>> e2l;
        rebuildLeafEdges(e2l);
        // snapshot the still-live flagged leaves (coarsening cannot have touched a
        // refine-flagged cell: coarsening needs ALL incident leaves flagged coarsen)
        std::vector<int> todo;
        for (int t : refineSet)
            if (tris_[t].leaf() && tris_[t].gen < maxGen) todo.push_back(t);
        std::sort(todo.begin(), todo.end());
        for (int t : todo) {
            if (!tris_[t].leaf()) continue;                 // already refined by a neighbour's closure
            refineElement(t, maxGen, e2l);
            ++nRefine;
        }
    }
    return {nRefine, nCoarsen};
}

// ===========================================================================
// Compact the live leaves into a flat conforming Mesh + record leaf order.
// ===========================================================================
void AdaptiveForest::buildMesh(Mesh& mesh) {
    leafOrder_.clear();
    for (int t = 0; t < static_cast<int>(tris_.size()); ++t)
        if (tris_[t].leaf()) leafOrder_.push_back(t);

    std::unordered_map<int, int> remap;                     // forest vid -> compact node id
    remap.reserve(leafOrder_.size() * 2);
    int nN = 0;
    for (int t : leafOrder_)
        for (int k = 0; k < 3; ++k) {
            auto it = remap.find(tris_[t].v[k]);
            if (it == remap.end()) remap.emplace(tris_[t].v[k], nN++);
        }
    mesh.node.resize(nN, 2);
    for (auto& [vid, nid] : remap) mesh.node.row(nid) = verts_[vid].transpose();
    int NT = static_cast<int>(leafOrder_.size());
    mesh.elem.resize(NT, 3);
    for (int i = 0; i < NT; ++i) {
        const auto& v = tris_[leafOrder_[i]].v;
        mesh.elem(i, 0) = remap[v[0]]; mesh.elem(i, 1) = remap[v[1]]; mesh.elem(i, 2) = remap[v[2]];
    }
}

// ===========================================================================
// Cheap structural conformity check: no edge may carry >2 leaves.  (The driver
// does the strong geometric test -- every 1-sided edge lies on the domain
// boundary -- which actually detects hanging nodes.)
// ===========================================================================
bool AdaptiveForest::checkConforming(std::string& msg) const {
    std::map<std::pair<int, int>, int> cnt;
    for (int t = 0; t < static_cast<int>(tris_.size()); ++t) {
        if (!tris_[t].leaf()) continue;
        const auto& v = tris_[t].v;
        cnt[key(v[0], v[1])]++; cnt[key(v[1], v[2])]++; cnt[key(v[2], v[0])]++;
    }
    int bad = 0;
    for (auto& [e, c] : cnt) if (c > 2) ++bad;
    if (bad) { std::ostringstream os; os << bad << " edges carry >2 leaves"; msg = os.str(); return false; }
    msg = "ok";
    return true;
}

} // namespace euler
