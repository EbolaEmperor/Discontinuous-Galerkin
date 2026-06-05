#ifndef NS_MESHGEN_H
#define NS_MESHGEN_H

#include <Eigen/Dense>
#include <vector>
#include "Mesh.h"

using namespace Eigen;

// ---------------------------------------------------------------------------
// Unstructured high-quality triangular mesh generator for the classic
// "flow past a cylinder" domain: an axis-aligned rectangle with a circular
// hole (the cylinder) removed.  The mesh is produced from a single target
// element size h by a DistMesh-style force-equilibrium iteration (Persson &
// Strang, 2004) on top of a from-scratch Bowyer-Watson Delaunay triangulator.
//
// The element size is graded: ~ h on the cylinder surface, growing smoothly to
// ~ size_far away from it, so the near-wake (where the von Karman street lives)
// is well resolved while the far field stays cheap.
// ---------------------------------------------------------------------------

// Rectangle [xa,xb] x [ya,yb] with a disk of radius r centred at (cx,cy) removed.
struct CylinderGeom {
    double xa, xb, ya, yb;   // outer rectangle
    double cx, cy, r;        // cylinder centre + radius

    // Signed distance to the fluid domain (negative inside the fluid, i.e.
    // inside the rectangle AND outside the disk).  d = max(d_rect, -d_circ).
    double sdist(double x, double y) const {
        double dRect = -std::min(std::min(y - ya, yb - y), std::min(x - xa, xb - x));
        double dCirc = std::hypot(x - cx, y - cy) - r;
        return std::max(dRect, -dCirc);
    }
};

// Boundary-edge classification used by the NS solver to apply the right BC.
enum BdryTag {
    BD_INTERIOR = 0, // shared by two elements
    BD_INFLOW   = 1, // left  wall  x = xa : Dirichlet velocity (free stream)
    BD_OUTFLOW  = 2, // right wall  x = xb : do-nothing (Neumann u, Dirichlet p=0)
    BD_WALL     = 3, // top/bottom y=ya,yb : slip (v.n = 0, free tangential)
    BD_CYL      = 4  // cylinder surface    : no-slip (u = 0)
};

// Generate the mesh into `mesh` (CCW-oriented triangles).  `h` is the target
// edge length on the cylinder; `sizeFarRatio` is size_far / h in the far field;
// `gradeRate` controls how fast the size grows with distance from the cylinder.
// Returns the number of DistMesh iterations actually taken.
int generateCylinderMesh(Mesh& mesh, const CylinderGeom& g, double h,
                         double sizeFarRatio = 6.0, double gradeRate = 0.25,
                         int maxIter = 600, bool verbose = true);

// Classify every edge of `mesh` (rows of `edge`) as one of BdryTag.  Boundary
// edges (one neighbour == -1) are matched to the nearest piece of the geometry.
VectorXi classifyEdges(const Mesh& mesh, const MatrixXi& edge,
                       const MatrixXi& edge2side, const CylinderGeom& g);

// Stand-alone Bowyer-Watson Delaunay triangulation of a 2-D point cloud.
// Fills `tris` with CCW triangles (indices into `pts`).  Triangles touching the
// auxiliary super-triangle are discarded.
void delaunayTriangulate(const std::vector<Vector2d>& pts,
                         std::vector<Vector3i>& tris);

// Mesh quality report: minimum/mean triangle angle (degrees), area extremes.
void meshQuality(const Mesh& mesh, double& minAngleDeg, double& meanAngleDeg,
                 double& minArea, double& maxArea);

#endif
