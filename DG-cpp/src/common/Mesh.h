#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

class Mesh {
public:
    MatrixXd node; // N_node x 2
    MatrixXi elem; // N_elem x 3

    Mesh() {}
    
    // Generates square mesh on [0,1]^2 with spacing h
    void getMesh(double h);
    // Generates polygon fan mesh with regular refinements toward target h
    void getPolygonMesh(const MatrixXd& vertices, double h);
    // Regular-hexagon domain: 6 equilateral fan triangles, uniformly refined until
    // the element edge <= h (so ~ one element per h, consistent with getMesh).
    void getHexagonMesh(const Vector2d& center, double circumradius, double h);
    // Disk domain, high-quality mesh auto-generated from the element size h: builds
    // the hexagon mesh (equilateral base) and smoothly maps it onto the circle --
    // each node keeps its angle and its radius is rescaled by R / r_hex(theta), so
    // the (<=15%) distortion is spread over the whole radius and the boundary lands
    // exactly on the circle.
    void getDiskMesh(const Vector2d& center, double radius, double h);
    // One step of uniform refinement (splits each triangle into 4)
    void uniformRefine();
    
    // Returns indices of boundary nodes (returns 1 if boundary, 0 otherwise)
    VectorXi findBdryNodes(const MatrixXd& nodes);

    // Returns edge information
    // edge: N_edge x 2 (node indices)
    // edge2side: N_edge x 2 (element indices)
    void getEdge2Side(MatrixXi& edge, MatrixXi& edge2side);
};

#endif
