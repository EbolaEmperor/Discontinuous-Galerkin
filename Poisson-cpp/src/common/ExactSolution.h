#ifndef EXACTSOLUTION_H
#define EXACTSOLUTION_H

#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

class ExactSolution {
public:
    double center;
    ExactSolution(double c = 0.3) : center(c) {}

    VectorXd u_exact(const MatrixXd& p) const {
        VectorXd val(p.rows());
        for(int i=0; i<p.rows(); ++i) {
            val(i) = std::sin(M_PI*(p(i,0)-center)) * std::sin(M_PI*(p(i,1)-center));
        }
        return val;
    }

    MatrixXd grad_u_exact(const MatrixXd& p) const {
        MatrixXd val(p.rows(), 2);
        for(int i=0; i<p.rows(); ++i) {
            val(i,0) = M_PI * std::cos(M_PI*(p(i,0)-center)) * std::sin(M_PI*(p(i,1)-center));
            val(i,1) = M_PI * std::sin(M_PI*(p(i,0)-center)) * std::cos(M_PI*(p(i,1)-center));
        }
        return val;
    }
    
    VectorXd laplace_u_exact(const MatrixXd& p) const {
        VectorXd val(p.rows());
        double f = -2 * M_PI * M_PI;
        for(int i=0; i<p.rows(); ++i) {
             val(i) = f * std::sin(M_PI*(p(i,0)-center)) * std::sin(M_PI*(p(i,1)-center));
        }
        return val;
    }

    VectorXd lap_lap_u_exact(const MatrixXd& p) const {
        VectorXd val(p.rows());
        double f = 4 * M_PI * M_PI * M_PI * M_PI;
        for (int i = 0; i < p.rows(); ++i) {
            val(i) = f * std::sin(M_PI * (p(i, 0) - center)) * std::sin(M_PI * (p(i, 1) - center));
        }
        return val;
    }
};

#endif
