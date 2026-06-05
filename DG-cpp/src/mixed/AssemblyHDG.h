#ifndef ASSEMBLY_HDG_H
#define ASSEMBLY_HDG_H

#include <Eigen/Sparse>
#include "VecFEM.h"
#include "FEM.h"
#include "Mesh.h"
#include "ExactSolution.h"

Eigen::SparseMatrix<double> assembleMass(VecFEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof);

Eigen::SparseMatrix<double> assembleDivMass(VecFEM& femSigma, FEM& femU, const Mesh& mesh,
                                            const Eigen::MatrixXi& elem2dofSigma, const Eigen::MatrixXi& elem2dofU);

Eigen::SparseMatrix<double> assembleGradMass(FEM& femU, VecFEM& femSigma, const Mesh& mesh,
                                             const Eigen::MatrixXi& elem2dofU, const Eigen::MatrixXi& elem2dofSigma);

Eigen::SparseMatrix<double> assembleIP_HDG(VecFEM& femSigma, FEM& femU, const Mesh& mesh,
                                           const Eigen::MatrixXi& elem2dofSigma, const Eigen::MatrixXi& elem2dofU,
                                           const Eigen::MatrixXi& edge, const Eigen::MatrixXi& edge2side,
                                           double alpha);

Eigen::VectorXd assembleWeakBDC(VecFEM& femSigma, FEM& femU, const Mesh& mesh,
                                const Eigen::MatrixXi& elem2dofSigma, const Eigen::MatrixXi& elem2dofU,
                                const Eigen::MatrixXi& edge, const Eigen::MatrixXi& edge2side,
                                double alpha, const ExactSolution& sol);

double l2ErrorScalar(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
                     const Eigen::VectorXd& uh, const ExactSolution& sol);

double l2ErrorVector(VecFEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
                     const Eigen::VectorXd& sigmah, const ExactSolution& sol);

#endif
