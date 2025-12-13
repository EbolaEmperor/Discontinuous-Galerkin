#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <Eigen/Dense>

using namespace Eigen;

class Quadrature {
public:
    static void quadpts(int order, MatrixXd& lambda, VectorXd& weight);
    static void quadpts1_my(int n, MatrixXd& x, VectorXd& w);
    static void quadpts2_my(int order, MatrixXd& lambda, VectorXd& weight);
};

#endif
