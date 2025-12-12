#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <Eigen/Sparse>
#include "FEM.h"
#include "Mesh.h"
#include "ExactSolution.h"

using namespace Eigen;

SparseMatrix<double> assembleK_Poi2D(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof);

SparseMatrix<double> assembleIP_Poi2D(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, 
                                      const MatrixXi& edge, const MatrixXi& edge2side, 
                                      double sigma, double beta);

VectorXd assembleLoadVector(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, const ExactSolution& sol);

void interpStrongBDC(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, 
                     const ExactSolution& sol, VectorXd& c, std::vector<int>& freeDof);

void getH1Err(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, const VectorXd& c, 
              const ExactSolution& sol, double& errH1, double& errL2);

#endif

