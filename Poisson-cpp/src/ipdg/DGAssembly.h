#ifndef DG_ASSEMBLY_H
#define DG_ASSEMBLY_H

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
                     const ExactSolution& sol, VectorXd& c, std::vector<int>& freeDof,
                     const MatrixXd* polygonVertices = nullptr);

void getH1Err(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, const VectorXd& c, 
              const ExactSolution& sol, double& errH1, double& errL2);

// Biharmonic C0-IPCG helpers
SparseMatrix<double> assembleK_Bihar2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof);
SparseMatrix<double> assembleIP_Bihar2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                        const MatrixXi& edge, const MatrixXi& edge2side,
                                        double sigma, double beta);
VectorXd assembleLoadVectorBihar(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                 const ExactSolution& sol);
VectorXd modifyLoadVectorNitsche_Bihar2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                         const MatrixXi& edge, const MatrixXi& edge2side,
                                         double sigma, double beta, const ExactSolution& sol);

#endif
