#ifndef VECFEM_H
#define VECFEM_H

#include <Eigen/Dense>
#include <vector>
#include "Mesh.h"
#include "FEM.h"
#include "utils.h"

class VecFEM {
public:
    int ord;
    int locDof;
    int scalarLocDof;
    std::vector<Eigen::MatrixXd> Dlam; // NT items, each 2x3
    Eigen::VectorXd area;              // NT

    VecFEM(int order, const Mesh& mesh);

    void getDOF(const Mesh& mesh, Eigen::MatrixXi& elem2dof, int& nDof) const;

    // 2 x locDof (first scalar block -> x-component, second -> y-component)
    Eigen::MatrixXd computeBasisValue_all(const Eigen::Vector3d& lam) const;

    // Divergence shape functions: [d/dx phi_1 ... d/dx phi_n, d/dy phi_1 ... d/dy phi_n]
    Eigen::RowVectorXd computeBasisDiv_all(int tid, const Eigen::Vector3d& lam) const;

    void quad2d(Eigen::MatrixXd& quadL, Eigen::VectorXd& w) const;
    void quad1d(Eigen::MatrixXd& quadL, Eigen::VectorXd& w) const;
};

#endif
