#ifndef FEM_H
#define FEM_H

#include <Eigen/Dense>
#include <vector>
#include "Mesh.h"

using namespace Eigen;

class FEM {
public:
    int ord;
    int locDof;
    MatrixXd coef; // locDof x locDof
    std::vector<MatrixXd> Dlam; // NT items, each 2x3
    VectorXd area; // NT

    FEM(int order, Mesh& mesh);

    void getDOF(const Mesh& mesh, MatrixXi& elem2dof, int& nDof);
    
    // Basis functions on reference element
    // lam: n_points x 3
    // Returns n_points x locDof
    MatrixXd computeBasisValue_all(const MatrixXd& lam);

    // Derivatives of basis functions w.r.t barycentric coordinates
    // lam: 1 x 3 (single point)
    // Returns 3 x locDof
    MatrixXd computeBasisDlam_all(const Vector3d& lam);
    
    // Helper to get quadrature points for this order
    // Returns matrix nq x 3 and weights nq x 1
    void quad2d(MatrixXd& quadL, VectorXd& w);
    void quad1d(MatrixXd& quadL, VectorXd& w);

private:
    MatrixXd initBasis();
    MatrixXd getSpan(const MatrixXd& lam);
};

// Computes gradients of barycentric coordinates and element areas
void gradbasis_my(const Mesh& mesh, std::vector<MatrixXd>& Dlam, VectorXd& area);

#endif

