#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <Eigen/Dense>
#include <functional>
#include "VecFEM.h"
#include "FEM.h"
#include "Mesh.h"

Eigen::VectorXd solveLocalPoisson(FEM& femStar, const VecFEM& femSigma, FEM& femU,
                                  const Mesh& mesh, const Eigen::MatrixXi& elem2dofStar,
                                  const Eigen::MatrixXi& elem2dofSigma, const Eigen::MatrixXi& elem2dofU,
                                  const Eigen::VectorXd& sigmah, const Eigen::VectorXd& uh,
                                  const std::function<double(const Eigen::Vector2d&)>& f);

#endif
