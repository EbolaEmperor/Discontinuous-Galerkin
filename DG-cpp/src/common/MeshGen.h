#ifndef COMMON_MESHGEN_H
#define COMMON_MESHGEN_H

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include "Mesh.h"

using namespace Eigen;

// ---------------------------------------------------------------------------
// Generic unstructured triangular mesh generator.
//
// The core routine below meshes any 2-D domain described by a signed-distance
// function and a target-size function.  Geometry-specific helpers such as
// generateCylinderMesh and generateBowlMesh are thin wrappers around this
// generic DistMesh-style engine.
// ---------------------------------------------------------------------------

struct DistanceMeshSpec {
    // Axis-aligned sampling box that contains the domain.
    double xa = 0.0, xb = 1.0, ya = 0.0, yb = 1.0;

    // Baseline element size used for tolerances and legacy uniform seeding.
    double h0 = 0.1;

    // Optional initial hex-lattice spacing.  Set this to the finest requested
    // size when targetSize is much smaller than h0 near local features.
    double seedH = 0.0;

    // Negative inside the fluid domain, zero on the boundary, positive outside.
    std::function<double(double, double)> signedDistance;

    // Desired local edge length.
    std::function<double(double, double)> targetSize;

    // Points held fixed by the mesh smoother.  Put corners, exact boundary
    // samples, and any geometry anchors here.
    std::vector<Vector2d> fixedPoints;

    // Optional projection for escaped interior points.  If omitted, the generic
    // numerical signed-distance gradient projection is used.
    std::function<Vector2d(const Vector2d&)> projectInside;
};

// Generate a CCW triangular mesh into `mesh`.  Returns the number of
// force-equilibrium iterations actually taken.
int generateDistanceMesh(Mesh& mesh, const DistanceMeshSpec& spec,
                         int maxIter = 600, bool verbose = true);

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
    double sdist(double x, double y) const;
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

// ---------------------------------------------------------------------------
// Closed circular "bowl" domain: a single disk of radius R centred at (cx,cy),
// with no inflow/outflow -- the whole boundary is a stationary no-slip wall.
// Used by the spoon-stirred-soup FSI demo.  The interior is meshed nearly
// uniformly (the shed vortex dipole roams the whole bowl), with optional mild
// grading that keeps elements finest inside the stirring radius rFocus.
// ---------------------------------------------------------------------------
struct BowlGeom {
    double cx, cy, R;
    // Signed distance to the fluid domain (negative inside the disk).
    double sdist(double x, double y) const;
};

// Boundary tag for the bowl: interior, or the no-slip circular wall.
enum BowlBdryTag { BD_BOWL_INTERIOR = 0, BD_BOWL_WALL = 1 };

// Generate a high-quality triangular mesh of the disk by the same DistMesh
// force-equilibrium iteration used for the cylinder domain.  The mesh is
// ADAPTIVELY refined to an annular band [bandLo, bandHi] (absolute radii) -- the
// region swept by the stirring spoon and traversed by the shed dipole: inside the
// band the target edge length is `hFine`, and outside it grows with distance from
// the band (rate `gradeRate`) up to `hFine*farRatio` in the quiescent core / rim.
// This is the same distance-graded "adaptive" sizing the cylinder mesh uses,
// applied to the spoon's known trajectory instead of a cylinder wall.
int generateBowlMesh(Mesh& mesh, const BowlGeom& g, double hFine,
                     double farRatio = 1.9, double gradeRate = 0.2,
                     double bandLo = 0.0, double bandHi = 1e9,
                     int maxIter = 600, bool verbose = true);

// Classify every edge of a bowl mesh: BD_BOWL_WALL on the boundary, else
// BD_BOWL_INTERIOR.
VectorXi classifyBowlEdges(const Mesh& mesh, const MatrixXi& edge,
                           const MatrixXi& edge2side, const BowlGeom& g);

// Stand-alone Bowyer-Watson Delaunay triangulation of a 2-D point cloud.
// Fills `tris` with CCW triangles (indices into `pts`).  Triangles touching the
// auxiliary super-triangle are discarded.
void delaunayTriangulate(const std::vector<Vector2d>& pts,
                         std::vector<Vector3i>& tris);

// Mesh quality report: minimum/mean triangle angle (degrees), area extremes.
void meshQuality(const Mesh& mesh, double& minAngleDeg, double& meanAngleDeg,
                 double& minArea, double& maxArea);

#endif
