#include "DGAssembly.h"

SparseMatrix<double> assembleK_Poi2D(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    
    std::vector<Triplet<double>> triplets;
    triplets.reserve(NT * locDof * locDof);
    
    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = w.size();
    
    // Precompute reference gradients dphi/dlam at quad points
    std::vector<MatrixXd> dphi_dlam_q(nq);
    for(int q=0; q<nq; ++q) {
        dphi_dlam_q[q] = fem.computeBasisDlam_all(quadL.row(q));
    }
    
    MatrixXd At(locDof, locDof);
    MatrixXd grad_phi(2, locDof);  // Reuse matrix
    for (int t = 0; t < NT; ++t) {
        At.setZero();
        const MatrixXd& Dlam = fem.Dlam[t]; // 2 x 3
        double area_t = fem.area(t);
        
        for (int q = 0; q < nq; ++q) {
            // grad_phi = Dlam * dphi_dlam  (2x3 * 3xlocDof = 2xlocDof)
            grad_phi.noalias() = Dlam * dphi_dlam_q[q];
            
            double weight = w(q) * area_t;
            // At += weight * grad_phi^T * grad_phi
            // More efficient: compute outer product directly
            At.noalias() += weight * grad_phi.transpose() * grad_phi;
        }
        
        for (int i = 0; i < locDof; ++i) {
            for (int j = 0; j < locDof; ++j) {
                triplets.push_back(Triplet<double>(elem2dof(t, i), elem2dof(t, j), At(i, j)));
            }
        }
    }
    
    SparseMatrix<double> A(nDof, nDof);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

SparseMatrix<double> assembleIP_Poi2D(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, 
                                      const MatrixXi& edge, const MatrixXi& edge2side, 
                                      double sigma, double beta) {
    int NE = edge.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    
    std::vector<Triplet<double>> triplets;
    // Estimate triplets count. Assume most edges are internal.
    // Each edge has 4 blocks of locDof*locDof.
    triplets.reserve(NE * 4 * locDof * locDof);
    
    MatrixXd quad1d;
    VectorXd w1d;
    fem.quad1d(quad1d, w1d);
    int nq = w1d.size();
    
    // Precompute basis functions for all possible edge configurations
    // In a triangle, edges can be (0,1), (1,2), or (2,0)
    // For each edge type, we need both forward and reverse order
    // phi_edge_q[edge_type][direction][q] where direction: 0=forward, 1=reverse
    std::vector<std::vector<std::vector<MatrixXd>>> phi_edge_q(3);  // 3 edge types, 2 directions
    std::vector<std::vector<std::vector<MatrixXd>>> dphi_edge_q(3); // 3 edge types, 2 directions
    
    for (int edge_type = 0; edge_type < 3; ++edge_type) {
        phi_edge_q[edge_type].resize(2);  // 2 directions
        dphi_edge_q[edge_type].resize(2);
        
        int idx1, idx2;
        if (edge_type == 0) { idx1 = 0; idx2 = 1; }      // edge (0,1)
        else if (edge_type == 1) { idx1 = 1; idx2 = 2; }  // edge (1,2)
        else { idx1 = 2; idx2 = 0; }                      // edge (2,0)
        
        for (int dir = 0; dir < 2; ++dir) {
            phi_edge_q[edge_type][dir].resize(nq);
            dphi_edge_q[edge_type][dir].resize(nq);
            
            for (int q = 0; q < nq; ++q) {
                double lam_e1 = quad1d(q, 0);
                double lam_e2 = quad1d(q, 1);
                
                Vector3d lam = Vector3d::Zero();
                if (dir == 0) {
                    // Forward direction: (idx1, idx2)
                    lam(idx1) = lam_e1;
                    lam(idx2) = lam_e2;
                } else {
                    // Reverse direction: (idx2, idx1) - swap lam_e1 and lam_e2
                    lam(idx1) = lam_e2;
                    lam(idx2) = lam_e1;
                }
                
                phi_edge_q[edge_type][dir][q] = fem.computeBasisValue_all(lam.transpose()); // 1 x locDof
                dphi_edge_q[edge_type][dir][q] = fem.computeBasisDlam_all(lam); // 3 x locDof
            }
        }
    }
    
    MatrixXd M11(locDof, locDof);
    MatrixXd M12(locDof, locDof);
    MatrixXd M21(locDof, locDof);
    MatrixXd M22(locDof, locDof);

    for (int e = 0; e < NE; ++e) {
        int t1 = edge2side(e, 0);
        int t2 = edge2side(e, 1);
        
        // Skip boundary edges for IP
        if (t1 == -1 || t2 == -1) continue;
        
        // Vertices of edge
        int n1 = edge(e, 0);
        int n2 = edge(e, 1);
        
        const Vector2d& p1 = mesh.node.row(n1);
        const Vector2d& p2 = mesh.node.row(n2);
        Vector2d evec = p2 - p1;
        double he = evec.norm();
        Vector2d nvec; nvec << evec(1), -evec(0); // Normal (dy, -dx)
        nvec /= he;
        
        // Find vertex indices - optimized: check all at once
        int idx1_1 = -1, idx1_2 = -1;
        int idx2_1 = -1, idx2_2 = -1;
        
        const auto& elem1 = mesh.elem.row(t1);
        const auto& elem2 = mesh.elem.row(t2);
        for(int k=0; k<3; ++k) {
            if (elem1(k) == n1) idx1_1 = k;
            else if (elem1(k) == n2) idx1_2 = k;
            if (elem2(k) == n1) idx2_1 = k;
            else if (elem2(k) == n2) idx2_2 = k;
        }
        
        // Local matrices
        M11.setZero(); M12.setZero();
        M21.setZero(); M22.setZero();
        
        // Determine edge type and direction for t1 and t2
        // Edge type: 0=(0,1), 1=(1,2), 2=(2,0)
        // Direction: 0=forward, 1=reverse
        // Forward: (0,1), (1,2), (2,0)
        // Reverse: (1,0), (2,1), (0,2)
        auto getEdgeTypeAndDir = [](int idx1, int idx2) -> std::pair<int, int> {
            if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0)) {
                return {0, (idx1 > idx2) ? 1 : 0};
            } else if ((idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1)) {
                return {1, (idx1 > idx2) ? 1 : 0};
            } else { // (2,0) or (0,2)
                return {2, (idx1 == 0 && idx2 == 2) ? 1 : 0};
            }
        };
        
        auto [edge_type1, dir1] = getEdgeTypeAndDir(idx1_1, idx1_2);
        auto [edge_type2, dir2] = getEdgeTypeAndDir(idx2_1, idx2_2);
        
        // Precompute constants
        double penalty = sigma / he;
        double w_he = w1d(0) * he;  // Will be updated per quad point
        
        // Reusable matrices
        MatrixXd g1(2, locDof), g2(2, locDof);
        MatrixXd ngrad1(1, locDof), ngrad2(1, locDof);
        MatrixXd temp1(locDof, locDof), temp2(locDof, locDof);
        
        const MatrixXd& Dlam1 = fem.Dlam[t1];
        const MatrixXd& Dlam2 = fem.Dlam[t2];
        
        for (int q = 0; q < nq; ++q) {
            // Use precomputed basis functions with correct direction
            const MatrixXd& phi1 = phi_edge_q[edge_type1][dir1][q]; // 1 x locDof
            const MatrixXd& phi2 = phi_edge_q[edge_type2][dir2][q]; // 1 x locDof
            
            const MatrixXd& dphi1_dl = dphi_edge_q[edge_type1][dir1][q]; // 3 x locDof
            const MatrixXd& dphi2_dl = dphi_edge_q[edge_type2][dir2][q]; // 3 x locDof
            
            g1.noalias() = Dlam1 * dphi1_dl; // 2 x locDof
            g2.noalias() = Dlam2 * dphi2_dl; // 2 x locDof
            
            ngrad1.noalias() = nvec.transpose() * g1;
            ngrad2.noalias() = nvec.transpose() * g2;
            
            // Jumps and averages
            w_he = w1d(q) * he;
            
            // Precompute outer products that are reused
            MatrixXd phi1T_phi1 = phi1.transpose() * phi1;  // locDof x locDof
            MatrixXd phi1T_phi2 = phi1.transpose() * phi2;  // locDof x locDof
            MatrixXd phi2T_phi1 = phi2.transpose() * phi1;  // locDof x locDof
            MatrixXd phi2T_phi2 = phi2.transpose() * phi2;  // locDof x locDof
            
            // M11 += w_q * (-beta * 0.5 * ngrad1^T * phi1 - 0.5 * phi1^T * ngrad1 + penalty * phi1^T * phi1)
            temp1.noalias() = ngrad1.transpose() * phi1;
            temp2.noalias() = phi1.transpose() * ngrad1;
            M11.noalias() += w_he * (-beta * 0.5 * temp1 - 0.5 * temp2 + penalty * phi1T_phi1);
                          
            // M12 += w_q * (beta * 0.5 * ngrad1^T * phi2 - 0.5 * phi1^T * ngrad2 - penalty * phi1^T * phi2)
            temp1.noalias() = ngrad1.transpose() * phi2;
            temp2.noalias() = phi1.transpose() * ngrad2;
            M12.noalias() += w_he * (beta * 0.5 * temp1 - 0.5 * temp2 - penalty * phi1T_phi2);
                          
            // M21 += w_q * (-beta * 0.5 * ngrad2^T * phi1 + 0.5 * phi2^T * ngrad1 - penalty * phi2^T * phi1)
            temp1.noalias() = ngrad2.transpose() * phi1;
            temp2.noalias() = phi2.transpose() * ngrad1;
            M21.noalias() += w_he * (-beta * 0.5 * temp1 + 0.5 * temp2 - penalty * phi2T_phi1);
                          
            // M22 += w_q * (beta * 0.5 * ngrad2^T * phi2 + 0.5 * phi2^T * ngrad2 + penalty * phi2^T * phi2)
            temp1.noalias() = ngrad2.transpose() * phi2;
            temp2.noalias() = phi2.transpose() * ngrad2;
            M22.noalias() += w_he * (beta * 0.5 * temp1 + 0.5 * temp2 + penalty * phi2T_phi2);
        }
        
        // Assemble
        for(int i=0; i<locDof; ++i) {
            for(int j=0; j<locDof; ++j) {
                triplets.push_back(Triplet<double>(elem2dof(t1, i), elem2dof(t1, j), M11(i, j)));
                triplets.push_back(Triplet<double>(elem2dof(t1, i), elem2dof(t2, j), M12(i, j)));
                triplets.push_back(Triplet<double>(elem2dof(t2, i), elem2dof(t1, j), M21(i, j)));
                triplets.push_back(Triplet<double>(elem2dof(t2, i), elem2dof(t2, j), M22(i, j)));
            }
        }
    }
    
    SparseMatrix<double> P(nDof, nDof);
    P.setFromTriplets(triplets.begin(), triplets.end());
    return P;
}

VectorXd assembleLoadVector(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, const ExactSolution& sol) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    VectorXd F = VectorXd::Zero(nDof);
    
    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = w.size();
    
    // Precompute basis values at quadrature points
    std::vector<MatrixXd> phi_q(nq);
    for(int q=0; q<nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q)); // 1 x locDof
    }
    
    MatrixXd pt_mat(1, 2);  // Reuse matrix
    VectorXd Ft(locDof);
    
    for (int t = 0; t < NT; ++t) {
        Ft.setZero();
        
        // Vertices for mapping
        int i1 = mesh.elem(t, 0);
        int i2 = mesh.elem(t, 1);
        int i3 = mesh.elem(t, 2);
        const Vector2d& p1 = mesh.node.row(i1);
        const Vector2d& p2 = mesh.node.row(i2);
        const Vector2d& p3 = mesh.node.row(i3);
        double area_t = fem.area(t);
        
        for (int q = 0; q < nq; ++q) {
            const Vector3d& lam = quadL.row(q);
            Vector2d pt = lam(0)*p1 + lam(1)*p2 + lam(2)*p3;
            
            // Evaluate f
            pt_mat(0, 0) = pt(0);
            pt_mat(0, 1) = pt(1);
            VectorXd lap = sol.laplace_u_exact(pt_mat);
            double f_val = -lap(0); // f = - lap u
            
            const MatrixXd& phi = phi_q[q]; // 1 x locDof
            
            Ft.noalias() += w(q) * area_t * f_val * phi.transpose();
        }
        
        for (int i = 0; i < locDof; ++i) {
            F(elem2dof(t, i)) += Ft(i);
        }
    }
    return F;
}

static bool pointOnSegment(const Vector2d& p, const Vector2d& a, const Vector2d& b, double tol) {
    Vector2d v = b - a;
    double len = v.norm();
    double cross = v(0) * (p(1) - a(1)) - v(1) * (p(0) - a(0));
    double dotVal = (p - a).dot(p - b);
    return std::abs(cross) <= tol * std::max(1.0, len) && dotVal <= tol;
}

void interpStrongBDC(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, 
                     const ExactSolution& sol, VectorXd& c, std::vector<int>& freeDof,
                     const MatrixXd* polygonVertices) {
    int nDof = elem2dof.maxCoeff() + 1;
    c = VectorXd::Zero(nDof);
    std::vector<bool> isFixed(nDof, false);
    
    // Compute Lagrange nodes for reference element in structured order
    int locDof = fem.locDof;
    int ord = fem.ord;
    MatrixXd pnt(locDof, 3);
    if (ord == 0) {
        pnt << 1.0 / 3, 1.0 / 3, 1.0 / 3;
    } else {
        pnt << 1.0, 0.0, 0.0,
               0.0, 1.0, 0.0,
               0.0, 0.0, 1.0;
        int idx = 3;
        int nEdgeDof = ord - 1;
        for (int k = 1; k <= nEdgeDof; ++k) {
            double t = static_cast<double>(k) / ord;
            pnt.row(idx++) << 0.0, 1.0 - t, t; // edge [2,3]
        }
        for (int k = 1; k <= nEdgeDof; ++k) {
            double t = static_cast<double>(k) / ord;
            pnt.row(idx++) << t, 0.0, 1.0 - t; // edge [3,1]
        }
        for (int k = 1; k <= nEdgeDof; ++k) {
            double t = static_cast<double>(k) / ord;
            pnt.row(idx++) << 1.0 - t, t, 0.0; // edge [1,2]
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
    
    const double tol = 1e-12;
    MatrixXd pt_mat(1, 2);  // Reuse matrix
    int NT = mesh.elem.rows();
    bool usePolygon = (polygonVertices != nullptr && polygonVertices->rows() > 0);
    
    for (int t = 0; t < NT; ++t) {
        int i1 = mesh.elem(t, 0);
        int i2 = mesh.elem(t, 1);
        int i3 = mesh.elem(t, 2);
        const Vector2d& p1 = mesh.node.row(i1);
        const Vector2d& p2 = mesh.node.row(i2);
        const Vector2d& p3 = mesh.node.row(i3);
        
        for (int i = 0; i < locDof; ++i) {
            const Vector3d& lam = pnt.row(i);
            Vector2d pt = lam(0)*p1 + lam(1)*p2 + lam(2)*p3;
            
            bool onBd = false;
            if (usePolygon) {
                int M = polygonVertices->rows();
                for (int j = 0; j < M; ++j) {
                    Vector2d a = polygonVertices->row(j);
                    Vector2d b = polygonVertices->row((j + 1) % M);
                    if (pointOnSegment(pt, a, b, tol)) {
                        onBd = true; break;
                    }
                }
            } else {
                double x = pt(0), y = pt(1);
                onBd = (std::abs(x) < tol || std::abs(x-1) < tol ||
                        std::abs(y) < tol || std::abs(y-1) < tol);
            }

            if (onBd) {
                int gdof = elem2dof(t, i);
                if (!isFixed[gdof]) {  // Avoid duplicate work
                    isFixed[gdof] = true;
                    pt_mat(0, 0) = pt(0);
                    pt_mat(0, 1) = pt(1);
                    c(gdof) = sol.u_exact(pt_mat)(0);
                }
            }
        }
    }
    
    freeDof.clear();
    freeDof.reserve(nDof);  // Reserve space
    for (int i = 0; i < nDof; ++i) {
        if (!isFixed[i]) freeDof.push_back(i);
    }
}

void getH1Err(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, const VectorXd& c, 
              const ExactSolution& sol, double& errH1, double& errL2) {
    int NT = mesh.elem.rows();
    int locDof = fem.locDof;
    errH1 = 0;
    errL2 = 0;
    
    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = w.size();
    
    // Precompute basis values and gradients at quadrature points
    std::vector<MatrixXd> phi_q(nq);
    std::vector<MatrixXd> dphi_dl_q(nq);
    for(int q=0; q<nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q)); // 1 x locDof
        dphi_dl_q[q] = fem.computeBasisDlam_all(quadL.row(q));
    }
    
    VectorXd uh_coeffs(locDof);
    MatrixXd pt_mat(1, 2);
    
    for (int t = 0; t < NT; ++t) {
        // Element DoF values
        for(int i=0; i<locDof; ++i) uh_coeffs(i) = c(elem2dof(t, i));
        
        int i1 = mesh.elem(t, 0);
        int i2 = mesh.elem(t, 1);
        int i3 = mesh.elem(t, 2);
        Vector2d p1 = mesh.node.row(i1);
        Vector2d p2 = mesh.node.row(i2);
        Vector2d p3 = mesh.node.row(i3);
        
        MatrixXd Dlam = fem.Dlam[t];
        double area = fem.area(t);
        
        MatrixXd grad_phi(2, locDof);  // Reuse matrix
        VectorXd grad_u_h(2);  // Reuse vector
        
        for(int q=0; q<nq; ++q) {
            const Vector3d& lam = quadL.row(q);
            Vector2d pt = lam(0)*p1 + lam(1)*p2 + lam(2)*p3;
            
            pt_mat(0, 0) = pt(0);
            pt_mat(0, 1) = pt(1);
            double u_ex = sol.u_exact(pt_mat)(0);
            VectorXd grad_u_ex = sol.grad_u_exact(pt_mat).row(0);
            
            // Use precomputed basis values
            double u_h = phi_q[q].row(0).dot(uh_coeffs);
            
            // Use precomputed basis gradients
            // grad_phi = Dlam * dphi_dl  (2x3 * 3xlocDof = 2xlocDof)
            grad_phi.noalias() = Dlam * dphi_dl_q[q];
            grad_u_h.noalias() = grad_phi * uh_coeffs;
            
            double weight = w(q) * area;
            double diff_u = u_h - u_ex;
            errL2 += weight * diff_u * diff_u;  // More efficient than pow
            grad_u_h -= grad_u_ex;
            errH1 += weight * grad_u_h.squaredNorm();
        }
    }
    
    errH1 = std::sqrt(errL2 + errH1);
    errL2 = std::sqrt(errL2);
}
