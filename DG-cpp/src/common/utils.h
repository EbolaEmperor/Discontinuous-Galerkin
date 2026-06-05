#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

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

inline MatrixXd polyBasisHomoHess3D(int k, const Vector3d& p) {
    MatrixXi number = numSplit3(k); // 3 x N_basis
    int n_basis = number.cols();
    MatrixXd H(6, n_basis);

    auto pow_int = [](double base, int exp) {
        if (exp <= 0) return 1.0;
        return std::pow(base, exp);
    };

    double p1 = p(0), p2 = p(1), p3 = p(2);
    for (int j = 0; j < n_basis; ++j) {
        int a = number(0, j);
        int b = number(1, j);
        int c = number(2, j);

        double H11 = a * std::max(a - 1, 0) * pow_int(p1, a - 2) * pow_int(p2, b) * pow_int(p3, c);
        double H22 = b * std::max(b - 1, 0) * pow_int(p1, a) * pow_int(p2, b - 2) * pow_int(p3, c);
        double H33 = c * std::max(c - 1, 0) * pow_int(p1, a) * pow_int(p2, b) * pow_int(p3, c - 2);

        double H12 = a * b * pow_int(p1, a - 1) * pow_int(p2, b - 1) * pow_int(p3, c);
        double H13 = a * c * pow_int(p1, a - 1) * pow_int(p2, b) * pow_int(p3, c - 1);
        double H23 = b * c * pow_int(p1, a) * pow_int(p2, b - 1) * pow_int(p3, c - 1);

        H(0, j) = H11;
        H(1, j) = H22;
        H(2, j) = H33;
        H(3, j) = H12;
        H(4, j) = H13;
        H(5, j) = H23;
    }
    return H;
}

#endif
