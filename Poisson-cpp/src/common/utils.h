#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>

using namespace Eigen;

// Generates multi-indices for homogeneous polynomials of degree k in 3 variables
// Returns a 3 x N matrix where N = (k+1)(k+2)/2
inline MatrixXi numSplit3(int k) {
    int N = (k + 1) * (k + 2) / 2;
    MatrixXi number(3, N);
    int cnt = 0;
    for (int i = 1; i <= k + 1; ++i) {
        int n_indices = i;
        for (int j = 0; j < n_indices; ++j) {
            number(0, cnt + j) = k - i + 1;
            number(1, cnt + j) = i - 1 - j;
            number(2, cnt + j) = j;
        }
        cnt += i;
    }
    return number;
}

inline MatrixXd polyBasisHomo3D(int k, const MatrixXd& p) {
    MatrixXi number = numSplit3(k); // 3 x N_basis
    int n_points = p.rows();
    int n_basis = number.cols();
    MatrixXd val(n_points, n_basis);

    for (int i = 0; i < n_points; ++i) {
        for (int j = 0; j < n_basis; ++j) {
            double v = 1.0;
            if (number(0, j) > 0) v *= std::pow(p(i, 0), number(0, j));
            if (number(1, j) > 0) v *= std::pow(p(i, 1), number(1, j));
            if (number(2, j) > 0) v *= std::pow(p(i, 2), number(2, j));
            val(i, j) = v;
        }
    }
    return val;
}

inline MatrixXd polyBasisHomoGrad3D(int k, const Vector3d& p) {
    MatrixXi number = numSplit3(k); // 3 x N_basis
    int n_basis = number.cols();
    MatrixXd val(3, n_basis);

    for (int j = 0; j < n_basis; ++j) {
        // d/dlam1
        double v1 = number(0, j);
        if (v1 != 0) {
            v1 *= std::pow(p(0), std::abs(number(0, j) - 1));
            v1 *= std::pow(p(1), number(1, j));
            v1 *= std::pow(p(2), number(2, j));
        } else {
            v1 = 0;
        }

        // d/dlam2
        double v2 = number(1, j);
        if (v2 != 0) {
            v2 *= std::pow(p(0), number(0, j));
            v2 *= std::pow(p(1), std::abs(number(1, j) - 1));
            v2 *= std::pow(p(2), number(2, j));
        } else {
            v2 = 0;
        }

        // d/dlam3
        double v3 = number(2, j);
        if (v3 != 0) {
            v3 *= std::pow(p(0), number(0, j));
            v3 *= std::pow(p(1), number(1, j));
            v3 *= std::pow(p(2), std::abs(number(2, j) - 1));
        } else {
            v3 = 0;
        }
        
        val(0, j) = v1;
        val(1, j) = v2;
        val(2, j) = v3;
    }
    return val;
}

#endif

