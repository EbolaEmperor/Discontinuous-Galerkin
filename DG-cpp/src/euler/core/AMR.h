#ifndef ADAPTIVE_MESH_H
#define ADAPTIVE_MESH_H

#include <Eigen/Dense>
#include <array>
#include <map>
#include <vector>

#include "FEM.h"
#include "Mesh.h"

// ===========================================================================
// Conforming h-adaptive triangulation by Newest-Vertex Bisection (NVB), with
// local high-order DG solution transfer, for the DG-Euler solver.
//
// WHY NVB (and not non-conforming "hanging node" AMR): the DG solver builds
// its face stencil from Mesh::getEdge2Side, which assumes a CONFORMING mesh
// (every interior edge shared by EXACTLY two triangles).  A hanging node would
// make an interior edge look one-sided and be silently treated as a boundary.
// NVB keeps the mesh conforming through arbitrary refine/coarsen, so the entire
// proven flux / artificial-viscosity / limiter / IMEX machinery is reused
// byte-for-byte; only the mesh changes between remeshes.
//
// FOREST.  Every triangle is a node of a binary refinement forest.  A triangle
// is stored as (v[0]=peak, v[1], v[2]) with its REFINEMENT EDGE = (v[1],v[2])
// (the edge opposite the newest vertex / peak).  Bisecting it inserts the
// midpoint m of (v[1],v[2]) and creates two children that share the segment
// peak->m:
//        child0 = (m, v0, v1)   (ref edge (v0,v1))
//        child1 = (m, v2, v0)   (ref edge (v2,v0))
// so each child's new refinement edge is one of the two "legs" adjacent to the
// old peak -- the standard NVB child rule.  m becomes the peak of both children.
//
// BASE LABELLING.  makeRectMesh splits each cell [a,b,c]+[a,c,d] across the
// diagonal (a,c).  We label BOTH triangles with refinement edge = that diagonal
// (peak = the off-diagonal vertex).  Then every mesh edge is the refinement edge
// of BOTH its triangles (the diagonals) or of NEITHER (the grid lines), i.e. the
// mesh is "compatibly labelled" (Stevenson) -> the conformity-closure recursion
// terminates and the first bisection of any cell is a local red 1->4 split.
//
// CONFORMITY CLOSURE (refine).  To bisect a triangle T we must bisect its
// refinement-edge neighbour F simultaneously.  If F's refinement edge is not the
// shared edge we first (recursively) refine F; NVB guarantees that after F is
// bisected the child of F adjacent to the shared edge has it as ITS refinement
// edge, so a bounded recursion makes the edge "compatibly divisible" and we
// bisect both sides at one shared midpoint -> no hanging node, ever.
//
// COARSEN.  A midpoint m is removable iff EVERY current leaf incident to m has m
// as its peak (so nothing finer hangs off m) AND all are flagged for coarsening;
// then all parent pairs around m are un-bisected together.  This is exactly the
// inverse of bisection and keeps the mesh conforming.
//
// SOLUTION TRANSFER.  dP_k coefficients ride on the leaves.  Refinement restricts
// the parent polynomial to each child EXACTLY (a degree-k polynomial restricted
// to a sub-triangle is still degree k) -- zero error, conservative.  Coarsening
// L2-projects the child polynomials onto the parent space -- conservative (the
// constants are in the space).  Both reduce to CONSTANT precomputed matrices
// (T_restrict[c], T_coarsen[c]) because the child->parent barycentric maps B_c
// and the child/parent area ratio (=1/2 for a bisection) are fixed.
// ===========================================================================

namespace euler {

using namespace Eigen;

class AdaptiveForest {
public:
    // Build the forest from a structured criss-cross base mesh (from
    // makeRectMesh) and the DG order (for the transfer operators).  The base
    // mesh's two-triangles-per-cell topology is detected and NVB-labelled.
    AdaptiveForest(const Mesh& baseMesh, int order);

    // Number of live leaves (== rows of the current flat mesh).
    int numLeaves() const { return static_cast<int>(leafOrder_.size()); }
    // Generation (number of bisections from the base) of flat element i.
    int gen(int flatElem) const { return tris_[leafOrder_[flatElem]].gen; }

    // Copy the current DG state into the leaf coefficient store.  U is nDof x 4
    // and e2d(flatElem,i) the dof of local node i; flat element order matches the
    // last buildMesh().  Must be called before adapt() so the indicator and the
    // transfer act on the live solution.
    void syncFromState(const MatrixXd& U, const MatrixXi& elem2dof, int locDof);

    // One adaptation pass.  flag[flatElem]: +1 refine, -1 coarsen, 0 hold.
    // Coarsening is applied first (one level), then refinement (one level).
    // Returns the number of (refine, coarsen) operations applied.
    std::pair<int, int> adapt(const std::vector<int>& flag, int maxGen);

    // Rebuild the flat conforming mesh from the live leaves (dedup referenced
    // vertices, compact).  Records the leaf<->flat-element correspondence.
    void buildMesh(Mesh& mesh);

    // Assemble the nDof x 4 conservative state for the current flat mesh from the
    // leaf coefficients (call after buildMesh + getDOF on the new mesh).
    MatrixXd gatherState(const MatrixXi& elem2dof, int locDof, int nDof) const;

    // Conformity self-test: returns true iff the flat mesh has no hanging nodes
    // (every interior edge shared by exactly two leaves).  msg gets a diagnostic.
    bool checkConforming(std::string& msg) const;

private:
    struct Tri {
        std::array<int, 3> v;          // forest vertex ids; v[0]=peak, refedge=(v[1],v[2])
        std::array<int, 2> child{-1, -1};
        int parent = -1;
        int whichChild = -1;           // 0/1: which child of parent (selects transfer matrix)
        int gen = 0;
        bool alive = true;             // false once a coarsen un-bisects this triangle's parent
        // a live leaf is the truth of the current mesh; a non-alive node is an
        // orphaned former child (its parent was merged back into a leaf).
        bool leaf() const { return alive && child[0] < 0; }
    };

    // ---- forest state ----
    int ord_, locDof_, nBaseVerts_;
    std::vector<Vector2d> verts_;                 // forest vertices (grows monotonically)
    std::vector<Tri> tris_;                       // triangle pool
    std::vector<Matrix<double, Dynamic, 4>> coef_;// per-tri dP_k coeffs (valid on leaves)
    std::map<std::pair<int, int>, int> midmap_;   // (sorted vid pair) -> midpoint vid
    std::vector<int> leafOrder_;                  // flat element i -> tri id

    // ---- transfer operators (constant, precomputed from the order) ----
    Matrix3d B_[2];                               // child c bary -> parent bary
    MatrixXd Trestrict_[2];                       // childCoeff = Trestrict[c] * parentCoeff
    MatrixXd Tcoarsen_[2];                         // parentCoeff += 0.5*Tcoarsen[c]*childCoeff[c]

    // ---- helpers ----
    void buildTransferOps();
    int  midpoint(int a, int b);                  // get/create midpoint vertex of (a,b)
    void bisect(int t);                            // split leaf t along its refinement edge
    void refineElement(int t, int maxGen,          // NVB refine with conformity closure
                       std::map<std::pair<int, int>, std::array<int, 2>>& e2leaf);
    int  neighborAcross(int t, int a, int b,       // current leaf sharing edge (a,b) with t
                        const std::map<std::pair<int, int>, std::array<int, 2>>& e2leaf) const;
    void rebuildLeafEdges(std::map<std::pair<int, int>, std::array<int, 2>>& e2leaf) const;
    std::array<int, 2> refEdge(int t) const { return {tris_[t].v[1], tris_[t].v[2]}; }
    static std::pair<int, int> key(int a, int b) { return {std::min(a, b), std::max(a, b)}; }
};

} // namespace euler

#endif
