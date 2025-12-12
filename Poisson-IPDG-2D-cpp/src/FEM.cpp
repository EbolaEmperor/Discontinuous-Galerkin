#include "FEM.h"
#include "utils.h"
#include "Quadrature.h"
#include <iostream>

void gradbasis_my(const Mesh& mesh, std::vector<MatrixXd>& Dlam, VectorXd& area) {
    int NT = mesh.elem.rows();
    Dlam.resize(NT);
    area.resize(NT);
    
    for (int t = 0; t < NT; ++t) {
        // vertices indices
        int i1 = mesh.elem(t, 0);
        int i2 = mesh.elem(t, 1);
        int i3 = mesh.elem(t, 2);
        
        Vector2d p1 = mesh.node.row(i1);
        Vector2d p2 = mesh.node.row(i2);
        Vector2d p3 = mesh.node.row(i3);
        
        Vector2d ve1 = p3 - p2;
        Vector2d ve2 = p1 - p3;
        Vector2d ve3 = p2 - p1;
        
        double at = 0.5 * (-ve3(0)*ve2(1) + ve3(1)*ve2(0));
        area(t) = at;
        
        MatrixXd Dl(2, 3);
        // Dlambda(:,1) = [-ve1.y; ve1.x] / (2*area)
        // Dlambda(:,2) = [-ve2.y; ve2.x] / (2*area)
        // Dlambda(:,3) = [-ve3.y; ve3.x] / (2*area)
        
        Dl(0, 0) = -ve1(1); Dl(1, 0) = ve1(0);
        Dl(0, 1) = -ve2(1); Dl(1, 1) = ve2(0);
        Dl(0, 2) = -ve3(1); Dl(1, 2) = ve3(0);
        
        Dl /= (2 * at);
        Dlam[t] = Dl;
        
        if (at < 0) {
            area(t) = -at;
            // MATLAB reverses sign of area but Dlambda is correct as is?
            // "The sign of Dlambda is always right since signed area is used"
            // So we divide by signed area 'at'.
            // And store absolute area in 'area'.
        }
    }
}

FEM::FEM(int order, Mesh& mesh) : ord(order) {
    locDof = (ord + 1) * (ord + 2) / 2;
    coef = initBasis();
    gradbasis_my(mesh, Dlam, area);
}

void FEM::getDOF(const Mesh& mesh, MatrixXi& elem2dof, int& nDof) {
    int NT = mesh.elem.rows();
    nDof = NT * locDof;
    elem2dof.resize(NT, locDof);
    
    int idx = 0;
    // For DG, DOFs are local to element, just numbered sequentially
    // But usually we want structure [locDof x NT] in MATLAB or similar?
    // MATLAB code: elem2dof = 1:NT*locDof; reshape...
    // elem2dof(t, :) are indices for element t.
    // 1-based in MATLAB, 0-based here.
    
    for (int t = 0; t < NT; ++t) {
        for (int i = 0; i < locDof; ++i) {
            elem2dof(t, i) = idx++;
        }
    }
}

MatrixXd FEM::initBasis() {
    // Computes coefficients for basis functions
    // Uses pnt as interpolation nodes
    int N = locDof;
    MatrixXd A(N, N);
    
    MatrixXd pnt(N, 3);
    if (ord == 0) {
        pnt << 1.0/3, 1.0/3, 1.0/3;
    } else {
        int ind = 0;
        for (int i = 0; i <= ord; ++i) {
            for (int j = 0; j <= ord - i; ++j) {
                pnt(ind, 0) = 1.0 - (double)(i + j) / ord;
                pnt(ind, 1) = (double)i / ord;
                pnt(ind, 2) = (double)j / ord;
                ind++;
            }
        }
    }
    
    for (int i = 0; i < N; ++i) {
        A.row(i) = getSpan(pnt.row(i));
    }
    
    // coef = inv(A)
    return A.inverse();
}

MatrixXd FEM::getSpan(const MatrixXd& lam) {
    // lam is 1 x 3 or N x 3
    // polyBasisHomo3D handles N x 3
    return polyBasisHomo3D(ord, lam);
}

MatrixXd FEM::computeBasisValue_all(const MatrixXd& lam) {
    MatrixXd span = getSpan(lam);
    return span * coef;
}

MatrixXd FEM::computeBasisDlam_all(const Vector3d& lam) {
    // Returns 3 x locDof
    // polyBasisHomoGrad3D returns 3 x locDof (if coef was identity)
    // We multiply by coef.
    // polyBasisHomoGrad3D(ord, lam) returns 3 x N_basis matrix.
    // coef is N_basis x locDof (N_basis == locDof).
    // result = poly_grad * coef
    
    MatrixXd poly_grad = polyBasisHomoGrad3D(ord, lam);
    return poly_grad * coef;
}

void FEM::quad2d(MatrixXd& quadL, VectorXd& w) {
    Quadrature::quadpts2_my(2 * (ord + 1), quadL, w);
}

void FEM::quad1d(MatrixXd& quadL, VectorXd& w) {
    Quadrature::quadpts1_my(2 * (ord + 1), quadL, w);
}

