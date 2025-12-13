#include "FEM.h"
#include "utils.h"
#include "Quadrature.h"
#include <iostream>
#include <map>
#include <cmath>

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
    buildR();
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

MatrixXd FEM::lagrangeNodes() const {
    if (ord == 0) {
        MatrixXd pnt(1, 3);
        pnt << 1.0 / 3, 1.0 / 3, 1.0 / 3;
        return pnt;
    }

    MatrixXd pnt(locDof, 3);
    // Vertices in fixed order
    pnt << 1.0, 0.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;

    if (ord > 1) {
        int nEdgeDof = ord - 1;
        int idx = 3;
        for (int k = 1; k <= nEdgeDof; ++k) {
            double t = static_cast<double>(k) / ord;
            // Edge 1: [2,3] => lam1 = 0
            pnt.row(idx++) << 0.0, 1.0 - t, t;
        }
        for (int k = 1; k <= nEdgeDof; ++k) {
            double t = static_cast<double>(k) / ord;
            // Edge 2: [3,1] => lam2 = 0
            pnt.row(idx++) << t, 0.0, 1.0 - t;
        }
        for (int k = 1; k <= nEdgeDof; ++k) {
            double t = static_cast<double>(k) / ord;
            // Edge 3: [1,2] => lam3 = 0
            pnt.row(idx++) << 1.0 - t, t, 0.0;
        }

        if (ord > 2) {
            for (int i = 1; i <= ord; ++i) {
                for (int j = 1; j <= ord - i; ++j) {
                    int k = ord - i - j;
                    if (k >= 1) {
                        pnt.row(idx++) << 1.0 - static_cast<double>(i + j) / ord,
                                          static_cast<double>(i) / ord,
                                          static_cast<double>(j) / ord;
                    }
                }
            }
        }
    }
    return pnt;
}

MatrixXd FEM::initBasis() {
    // Computes coefficients for basis functions using structured Lagrange nodes
    int N = locDof;
    MatrixXd A(N, N);

    MatrixXd pnt = lagrangeNodes();

    for (int i = 0; i < N; ++i) {
        A.row(i) = getSpan(pnt.row(i));
    }

    return A.inverse();
}

void FEM::getConformingDOF(const Mesh& mesh, MatrixXi& elem2dof, int& nDof, MatrixXd& dofCoords) const {
    int NT = mesh.elem.rows();
    int nNode = mesh.node.rows();
    int nEdgeDof = std::max(0, ord - 1);

    elem2dof.resize(NT, locDof);

    // 1. Vertices: directly use existing mesh node indices
    elem2dof.block(0, 0, NT, 3) = mesh.elem;

    // Prepare storage for global DOF coordinates
    std::vector<Vector2d> coords;
    coords.reserve(static_cast<size_t>(NT * locDof));
    for (int i = 0; i < nNode; ++i) {
        coords.push_back(mesh.node.row(i));
    }

    int currentIdx = nNode;

    // 2. Edges
    if (nEdgeDof > 0) {
        std::map<std::pair<int, int>, int> edgeMap;
        std::vector<std::pair<int, int>> edgeList;

        // Local edge node indices (0-based) following order [2,3], [3,1], [1,2]
        const int edgeNodes[3][2] = {{1, 2}, {2, 0}, {0, 1}};

        // Build unique edges
        for (int t = 0; t < NT; ++t) {
            for (int e = 0; e < 3; ++e) {
                int a = mesh.elem(t, edgeNodes[e][0]);
                int b = mesh.elem(t, edgeNodes[e][1]);
                int u = std::min(a, b);
                int v = std::max(a, b);
                std::pair<int, int> key{u, v};
                if (edgeMap.find(key) == edgeMap.end()) {
                    int idx = static_cast<int>(edgeList.size());
                    edgeMap[key] = idx;
                    edgeList.push_back(key);
                }
            }
        }

        // Create edge DOF coordinates in global direction (small -> big)
        for (const auto& e : edgeList) {
            const Vector2d p1 = mesh.node.row(e.first);
            const Vector2d p2 = mesh.node.row(e.second);
            for (int k = 1; k <= nEdgeDof; ++k) {
                double t = static_cast<double>(k) / ord;
                coords.push_back((1.0 - t) * p1 + t * p2);
            }
        }

        // Assign edge DOFs to elements respecting orientation
        for (int t = 0; t < NT; ++t) {
            for (int e = 0; e < 3; ++e) {
                int a = mesh.elem(t, edgeNodes[e][0]);
                int b = mesh.elem(t, edgeNodes[e][1]);
                int u = std::min(a, b);
                int v = std::max(a, b);
                int gEdgeIdx = edgeMap[{u, v}];
                bool isDirect = (a < b); // matches global small->big if true

                int base = currentIdx + gEdgeIdx * nEdgeDof;
                int colStart = 3 + e * nEdgeDof;
                for (int k = 0; k < nEdgeDof; ++k) {
                    int offset = isDirect ? k : (nEdgeDof - 1 - k);
                    elem2dof(t, colStart + k) = base + offset;
                }
            }
        }

        currentIdx += static_cast<int>(edgeList.size()) * nEdgeDof;
    }

    // 3. Face (interior) nodes
    if (ord > 2) {
        std::vector<Vector3d> internalLam;
        for (int i = 1; i <= ord; ++i) {
            for (int j = 1; j <= ord - i; ++j) {
                int k = ord - i - j;
                if (k >= 1) {
                    internalLam.emplace_back(1.0 - static_cast<double>(i + j) / ord,
                                             static_cast<double>(i) / ord,
                                             static_cast<double>(j) / ord);
                }
            }
        }
        int nInternal = static_cast<int>(internalLam.size());
        int startCol = 3 + 3 * nEdgeDof;

        for (int t = 0; t < NT; ++t) {
            int i1 = mesh.elem(t, 0);
            int i2 = mesh.elem(t, 1);
            int i3 = mesh.elem(t, 2);
            const Vector2d& p1 = mesh.node.row(i1);
            const Vector2d& p2 = mesh.node.row(i2);
            const Vector2d& p3 = mesh.node.row(i3);

            for (int k = 0; k < nInternal; ++k) {
                const Vector3d& lam = internalLam[k];
                coords.push_back(lam(0) * p1 + lam(1) * p2 + lam(2) * p3);
                elem2dof(t, startCol + k) = currentIdx + k;
            }
            currentIdx += nInternal;
        }
    }

    nDof = currentIdx;
    dofCoords.resize(nDof, 2);
    for (int i = 0; i < nDof; ++i) {
        dofCoords(i, 0) = coords[i](0);
        dofCoords(i, 1) = coords[i](1);
    }
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

MatrixXd FEM::computeBasisGrad_all(int tid, const Vector3d& lam) {
    const MatrixXd dphi_dl = computeBasisDlam_all(lam);
    return Dlam[tid] * dphi_dl;
}

MatrixXd FEM::computeBasisHlam_all(const Vector3d& lam) {
    MatrixXd poly_hess = polyBasisHomoHess3D(ord, lam);
    return poly_hess * coef;
}

MatrixXd FEM::computeBasisHessian_all(int tid, const Vector3d& lam) {
    MatrixXd Hlam = computeBasisHlam_all(lam); // 6 x locDof
    return R[tid] * Hlam; // 3 x locDof, with sqrt(2) scaling in the 3rd row
}

RowVectorXd FEM::computeBasisDirectedDiff2_all(int tid, const Vector3d& lam, const Vector2d& dir) {
    MatrixXd Hlam = computeBasisHlam_all(lam); // 6 x locDof
    Vector3d w = Dlam[tid].transpose() * dir;   // barycentric directional derivatives

    double w1 = w(0), w2 = w(1), w3 = w(2);
    RowVectorXd val = (w1 * w1) * Hlam.row(0)
                    + (w2 * w2) * Hlam.row(1)
                    + (w3 * w3) * Hlam.row(2)
                    + 2.0 * (w1 * w2) * Hlam.row(3)
                    + 2.0 * (w1 * w3) * Hlam.row(4)
                    + 2.0 * (w2 * w3) * Hlam.row(5);
    return val;
}

void FEM::quad2d(MatrixXd& quadL, VectorXd& w) {
    Quadrature::quadpts2_my(2 * (ord + 1), quadL, w);
}

void FEM::quad1d(MatrixXd& quadL, VectorXd& w) {
    Quadrature::quadpts1_my(2 * (ord + 1), quadL, w);
}

void FEM::buildR() {
    int NT = static_cast<int>(Dlam.size());
    R.resize(NT);
    const double sqrt2 = std::sqrt(2.0);
    for (int t = 0; t < NT; ++t) {
        const MatrixXd& Dl = Dlam[t]; // 2 x 3
        double a1 = Dl(0, 0), a2 = Dl(0, 1), a3 = Dl(0, 2);
        double b1 = Dl(1, 0), b2 = Dl(1, 1), b3 = Dl(1, 2);

        MatrixXd Rt(3, 6);
        // Hxx coefficients
        Rt(0, 0) = a1 * a1; Rt(0, 1) = a2 * a2; Rt(0, 2) = a3 * a3;
        Rt(0, 3) = 2 * a1 * a2; Rt(0, 4) = 2 * a1 * a3; Rt(0, 5) = 2 * a2 * a3;
        // Hyy coefficients
        Rt(1, 0) = b1 * b1; Rt(1, 1) = b2 * b2; Rt(1, 2) = b3 * b3;
        Rt(1, 3) = 2 * b1 * b2; Rt(1, 4) = 2 * b1 * b3; Rt(1, 5) = 2 * b2 * b3;
        // Hxy coefficients (with sqrt(2) scaling)
        Rt(2, 0) = a1 * b1; Rt(2, 1) = a2 * b2; Rt(2, 2) = a3 * b3;
        Rt(2, 3) = a1 * b2 + a2 * b1;
        Rt(2, 4) = a1 * b3 + a3 * b1;
        Rt(2, 5) = a2 * b3 + a3 * b2;
        Rt.row(2) *= sqrt2;
        R[t] = Rt;
    }
}
