#include "DGAssembly.h"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
double cross2d(const Vector2d& a, const Vector2d& b) {
    return a(0) * b(1) - a(1) * b(0);
}

std::pair<int, int> getEdgeTypeAndDir(int idx1, int idx2) {
    // Edge types: 0 = (0,1), 1 = (1,2), 2 = (2,0)
    if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0)) {
        return {0, (idx1 > idx2) ? 1 : 0};
    } else if ((idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1)) {
        return {1, (idx1 > idx2) ? 1 : 0};
    } else { // (2,0) or (0,2)
        return {2, (idx1 == 0 && idx2 == 2) ? 1 : 0};
    }
}

struct EdgeQuadData {
    std::vector<std::vector<std::vector<RowVectorXd>>> phi;
    std::vector<std::vector<std::vector<MatrixXd>>> dphi_dl;
    std::vector<std::vector<std::vector<MatrixXd>>> Hlam;
};

EdgeQuadData precomputeEdgeData(FEM& fem, const MatrixXd& quad1d) {
    int nq = quad1d.rows();
    EdgeQuadData data;
    data.phi.assign(3, std::vector<std::vector<RowVectorXd>>(2, std::vector<RowVectorXd>(nq)));
    data.dphi_dl.assign(3, std::vector<std::vector<MatrixXd>>(2, std::vector<MatrixXd>(nq)));
    data.Hlam.assign(3, std::vector<std::vector<MatrixXd>>(2, std::vector<MatrixXd>(nq)));

    for (int edge_type = 0; edge_type < 3; ++edge_type) {
        int idx1 = 0, idx2 = 1;
        if (edge_type == 1) { idx1 = 1; idx2 = 2; }
        else if (edge_type == 2) { idx1 = 2; idx2 = 0; }

        for (int dir = 0; dir < 2; ++dir) {
            for (int q = 0; q < nq; ++q) {
                double lam1 = quad1d(q, 0);
                double lam2 = quad1d(q, 1);

                Vector3d lam = Vector3d::Zero();
                if (dir == 0) {
                    lam(idx1) = lam1; lam(idx2) = lam2;
                } else {
                    lam(idx1) = lam2; lam(idx2) = lam1;
                }
                MatrixXd lamMat(1, 3);
                lamMat << lam.transpose();

                data.phi[edge_type][dir][q] = fem.computeBasisValue_all(lamMat).row(0);
                data.dphi_dl[edge_type][dir][q] = fem.computeBasisDlam_all(lam);
                data.Hlam[edge_type][dir][q] = fem.computeBasisHlam_all(lam);
            }
        }
    }
    return data;
}
}

SparseMatrix<double> assembleK_Poi2D(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    
    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = w.size();
    
    // Precompute reference gradients dphi/dlam at quad points
    std::vector<MatrixXd> dphi_dlam_q(nq);
    for(int q=0; q<nq; ++q) {
        dphi_dlam_q[q] = fem.computeBasisDlam_all(quadL.row(q));
    }

    int nThreads = 1;
#ifdef _OPENMP
    nThreads = omp_get_max_threads();
#endif
    std::vector<std::vector<Triplet<double>>> triplets_private(nThreads);
    int perThreadEstimate = static_cast<int>((NT / nThreads + 1) * locDof * locDof);
    for (auto& vec : triplets_private) {
        vec.reserve(perThreadEstimate);
    }

#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        auto& localTrip = triplets_private[tid];
        MatrixXd At(locDof, locDof);
        MatrixXd grad_phi(2, locDof);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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
                    localTrip.emplace_back(elem2dof(t, i), elem2dof(t, j), At(i, j));
                }
            }
        }
    }
    
    std::vector<Triplet<double>> triplets;
    size_t totalTriplets = 0;
    for (const auto& vec : triplets_private) totalTriplets += vec.size();
    triplets.reserve(totalTriplets);
    for (auto& vec : triplets_private) {
        triplets.insert(triplets.end(), vec.begin(), vec.end());
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
    
    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = w.size();
    
    // Precompute basis values at quadrature points
    std::vector<MatrixXd> phi_q(nq);
    for(int q=0; q<nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q)); // 1 x locDof
    }

    int nThreads = 1;
#ifdef _OPENMP
    nThreads = omp_get_max_threads();
#endif
    std::vector<VectorXd> F_private(nThreads, VectorXd::Zero(nDof));
    
#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        VectorXd& F_local = F_private[tid];
        MatrixXd pt_mat(1, 2);  // Thread-local matrix
        VectorXd Ft(locDof);
        
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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
                F_local(elem2dof(t, i)) += Ft(i);
            }
        }
    }
    
    VectorXd F = VectorXd::Zero(nDof);
    for (const auto& Floc : F_private) {
        F += Floc;
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

SparseMatrix<double> assembleK_Bihar2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;

    std::vector<Triplet<double>> trip;
    trip.reserve(NT * locDof * locDof);

    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<MatrixXd> Hlam_q(nq);
    for (int q = 0; q < nq; ++q) {
        Vector3d lam = quadL.row(q).transpose();
        Hlam_q[q] = fem.computeBasisHlam_all(lam); // 6 x locDof
    }

    MatrixXd At(locDof, locDof);
    for (int t = 0; t < NT; ++t) {
        At.setZero();
        double area = fem.area(t);
        const MatrixXd& Rt = fem.R[t];
        for (int q = 0; q < nq; ++q) {
            MatrixXd Hvoigt = Rt * Hlam_q[q]; // 3 x locDof
            At.noalias() += w(q) * area * (Hvoigt.transpose() * Hvoigt);
        }

        for (int i = 0; i < locDof; ++i) {
            int gi = elem2dof(t, i);
            for (int j = 0; j < locDof; ++j) {
                trip.emplace_back(gi, elem2dof(t, j), At(i, j));
            }
        }
    }

    SparseMatrix<double> A(nDof, nDof);
    A.setFromTriplets(trip.begin(), trip.end());
    return A;
}

SparseMatrix<double> assembleIP_Bihar2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                        const MatrixXi& edge, const MatrixXi& edge2side,
                                        double sigma, double beta) {
    int NE = edge.rows();
    int locDof = fem.locDof;
    int nDof = elem2dof.maxCoeff() + 1;

    std::vector<Triplet<double>> trip;
    trip.reserve(NE * 4 * locDof * locDof);

    MatrixXd quad1d;
    VectorXd w1d;
    fem.quad1d(quad1d, w1d);
    int nq = static_cast<int>(w1d.size());

    EdgeQuadData edgeData = precomputeEdgeData(fem, quad1d);
    const double sqrt2 = std::sqrt(2.0);

    for (int e = 0; e < NE; ++e) {
        int t1 = edge2side(e, 0);
        int t2 = edge2side(e, 1);

        int n1 = edge(e, 0);
        int n2 = edge(e, 1);
        Vector2d p1 = mesh.node.row(n1);
        Vector2d p2 = mesh.node.row(n2);
        Vector2d evec = p2 - p1;
        double he = evec.norm();
        Vector2d nvec(evec(1) / he, -evec(0) / he);

        bool hasLeft = (t1 != -1);
        bool hasRight = (t2 != -1);

        if (hasLeft && hasRight) {
            // Interior edge
            int idx1_1 = -1, idx1_2 = -1;
            int idx2_1 = -1, idx2_2 = -1;

            const auto& elem1 = mesh.elem.row(t1);
            const auto& elem2 = mesh.elem.row(t2);
            for (int k = 0; k < 3; ++k) {
                if (elem1(k) == n1) idx1_1 = k;
                else if (elem1(k) == n2) idx1_2 = k;
                if (elem2(k) == n1) idx2_1 = k;
                else if (elem2(k) == n2) idx2_2 = k;
            }

            auto [edge_type1, dir1] = getEdgeTypeAndDir(idx1_1, idx1_2);
            auto [edge_type2, dir2] = getEdgeTypeAndDir(idx2_1, idx2_2);

            MatrixXd locMat = MatrixXd::Zero(2 * locDof, 2 * locDof);
            double nx = nvec(0), ny = nvec(1);
            double coefXY = sqrt2 * nx * ny;

            for (int q = 0; q < nq; ++q) {
                const RowVectorXd& phi1 = edgeData.phi[edge_type1][dir1][q];
                const RowVectorXd& phi2 = edgeData.phi[edge_type2][dir2][q];
                (void)phi1; (void)phi2; // phi not used but kept for clarity

                const MatrixXd& dphi1_dl = edgeData.dphi_dl[edge_type1][dir1][q];
                const MatrixXd& dphi2_dl = edgeData.dphi_dl[edge_type2][dir2][q];

                MatrixXd g1 = fem.Dlam[t1] * dphi1_dl;
                MatrixXd g2 = fem.Dlam[t2] * dphi2_dl;
                RowVectorXd ngrad1 = nvec.transpose() * g1;
                RowVectorXd ngrad2 = nvec.transpose() * g2;

                const MatrixXd& Hlam1 = edgeData.Hlam[edge_type1][dir1][q];
                const MatrixXd& Hlam2 = edgeData.Hlam[edge_type2][dir2][q];
                MatrixXd Hvoigt1 = fem.R[t1] * Hlam1;
                MatrixXd Hvoigt2 = fem.R[t2] * Hlam2;

                RowVectorXd d2n1 = nx * nx * Hvoigt1.row(0) + ny * ny * Hvoigt1.row(1) + coefXY * Hvoigt1.row(2);
                RowVectorXd d2n2 = nx * nx * Hvoigt2.row(0) + ny * ny * Hvoigt2.row(1) + coefXY * Hvoigt2.row(2);

                RowVectorXd dn_jump(2 * locDof);
                dn_jump << ngrad1, -ngrad2;
                RowVectorXd d2n_mean(2 * locDof);
                d2n_mean << 0.5 * d2n1, 0.5 * d2n2;

                double w_he = w1d(q) * he;
                locMat.noalias() += w_he * (
                    -beta * d2n_mean.transpose() * dn_jump
                    - dn_jump.transpose() * d2n_mean
                    + (sigma / he) * dn_jump.transpose() * dn_jump
                );
            }

            VectorXi idx(2 * locDof);
            idx << elem2dof.row(t1).transpose(), elem2dof.row(t2).transpose();
            for (int i = 0; i < idx.size(); ++i) {
                for (int j = 0; j < idx.size(); ++j) {
                    trip.emplace_back(idx(i), idx(j), locMat(i, j));
                }
            }
        } else if (hasLeft || hasRight) {
            // Boundary edge
            int tid = hasLeft ? t1 : t2;
            const auto& elemRow = mesh.elem.row(tid);
            int idx_a = -1, idx_b = -1, idx_w = -1;
            for (int k = 0; k < 3; ++k) {
                int v = elemRow(k);
                if (v == n1) idx_a = k;
                else if (v == n2) idx_b = k;
                else idx_w = k;
            }

            Vector2d wpt = mesh.node.row(elemRow(idx_w));
            if (cross2d(evec, wpt - p1) < 0) nvec = -nvec;

            auto [edge_type, dir] = getEdgeTypeAndDir(idx_a, idx_b);
            MatrixXd locMat = MatrixXd::Zero(locDof, locDof);
            double nx = nvec(0), ny = nvec(1);
            double coefXY = sqrt2 * nx * ny;

            for (int q = 0; q < nq; ++q) {
                const MatrixXd& dphi_dl = edgeData.dphi_dl[edge_type][dir][q];
                MatrixXd g = fem.Dlam[tid] * dphi_dl;
                RowVectorXd ngrad = nvec.transpose() * g;

                const MatrixXd& Hlam = edgeData.Hlam[edge_type][dir][q];
                MatrixXd Hvoigt = fem.R[tid] * Hlam;
                RowVectorXd d2n = nx * nx * Hvoigt.row(0) + ny * ny * Hvoigt.row(1) + coefXY * Hvoigt.row(2);

                double w_he = w1d(q) * he;
                locMat.noalias() += w_he * (
                    -beta * d2n.transpose() * ngrad
                    - ngrad.transpose() * d2n
                    + (sigma / he) * ngrad.transpose() * ngrad
                );
            }

            const VectorXi idx = elem2dof.row(tid).transpose();
            for (int i = 0; i < locDof; ++i) {
                for (int j = 0; j < locDof; ++j) {
                    trip.emplace_back(idx(i), idx(j), locMat(i, j));
                }
            }
        }
    }

    SparseMatrix<double> P(nDof, nDof);
    P.setFromTriplets(trip.begin(), trip.end());
    return P;
}

VectorXd assembleLoadVectorBihar(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                 const ExactSolution& sol) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    VectorXd F = VectorXd::Zero(nDof);

    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<MatrixXd> phi_q(nq);
    for (int q = 0; q < nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q));
    }

    MatrixXd pt_mat(1, 2);
    VectorXd Ft(locDof);
    for (int t = 0; t < NT; ++t) {
        Ft.setZero();
        int i1 = mesh.elem(t, 0);
        int i2 = mesh.elem(t, 1);
        int i3 = mesh.elem(t, 2);
        Vector2d p1 = mesh.node.row(i1);
        Vector2d p2 = mesh.node.row(i2);
        Vector2d p3 = mesh.node.row(i3);
        double area = fem.area(t);

        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;

            pt_mat(0, 0) = pt(0);
            pt_mat(0, 1) = pt(1);
            double f_val = sol.lap_lap_u_exact(pt_mat)(0);
            const MatrixXd& phi = phi_q[q];
            Ft.noalias() += w(q) * area * f_val * phi.transpose();
        }

        for (int i = 0; i < locDof; ++i) {
            F(elem2dof(t, i)) += Ft(i);
        }
    }
    return F;
}

VectorXd modifyLoadVectorNitsche_Bihar2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                         const MatrixXi& edge, const MatrixXi& edge2side,
                                         double sigma, double beta, const ExactSolution& sol) {
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    VectorXd Fb = VectorXd::Zero(nDof);

    MatrixXd quad1d;
    VectorXd w1d;
    fem.quad1d(quad1d, w1d);
    int nq = static_cast<int>(w1d.size());
    EdgeQuadData edgeData = precomputeEdgeData(fem, quad1d);
    const double sqrt2 = std::sqrt(2.0);

    MatrixXd pt_mat(1, 2);
    for (int e = 0; e < edge.rows(); ++e) {
        int t1 = edge2side(e, 0);
        int t2 = edge2side(e, 1);
        if (!((t1 == -1) ^ (t2 == -1))) continue; // only boundary edges

        int tid = (t1 != -1) ? t1 : t2;
        int n1 = edge(e, 0);
        int n2 = edge(e, 1);
        Vector2d p1 = mesh.node.row(n1);
        Vector2d p2 = mesh.node.row(n2);
        Vector2d evec = p2 - p1;
        double he = evec.norm();
        Vector2d nvec(evec(1) / he, -evec(0) / he);

        const auto& elemRow = mesh.elem.row(tid);
        int idx_a = -1, idx_b = -1, idx_w = -1;
        for (int k = 0; k < 3; ++k) {
            int v = elemRow(k);
            if (v == n1) idx_a = k;
            else if (v == n2) idx_b = k;
            else idx_w = k;
        }

        Vector2d wpt = mesh.node.row(elemRow(idx_w));
        if (cross2d(evec, wpt - p1) < 0) nvec = -nvec;

        auto [edge_type, dir] = getEdgeTypeAndDir(idx_a, idx_b);
        double nx = nvec(0), ny = nvec(1);
        double coefXY = sqrt2 * nx * ny;

        VectorXd Ft = VectorXd::Zero(locDof);
        MatrixXd vtx(3, 2);
        vtx.row(0) = mesh.node.row(elemRow(0));
        vtx.row(1) = mesh.node.row(elemRow(1));
        vtx.row(2) = mesh.node.row(elemRow(2));

        for (int q = 0; q < nq; ++q) {
            const MatrixXd& dphi_dl = edgeData.dphi_dl[edge_type][dir][q];
            MatrixXd g = fem.Dlam[tid] * dphi_dl;
            RowVectorXd ngrad = nvec.transpose() * g;

            const MatrixXd& Hlam = edgeData.Hlam[edge_type][dir][q];
            MatrixXd Hvoigt = fem.R[tid] * Hlam;
            RowVectorXd d2n = nx * nx * Hvoigt.row(0) + ny * ny * Hvoigt.row(1) + coefXY * Hvoigt.row(2);

            // Physical point on the edge for exact gradient
            double lam1 = quad1d(q, 0);
            double lam2 = quad1d(q, 1);
            Vector3d lamVec = Vector3d::Zero();
            lamVec(idx_a) = lam1;
            lamVec(idx_b) = lam2;
            Vector2d pt = lamVec(0) * vtx.row(0).transpose()
                        + lamVec(1) * vtx.row(1).transpose()
                        + lamVec(2) * vtx.row(2).transpose();
            pt_mat(0, 0) = pt(0);
            pt_mat(0, 1) = pt(1);
            RowVectorXd gradExact = sol.grad_u_exact(pt_mat).row(0);
            double g1val = nvec.dot(gradExact.transpose());

            double w_he = w1d(q) * he;
            Ft.noalias() += w_he * ( -beta * d2n.transpose() * g1val + (sigma / he) * ngrad.transpose() * g1val );
        }

        const VectorXi idx = elem2dof.row(tid).transpose();
        for (int i = 0; i < locDof; ++i) {
            Fb(idx(i)) += Ft(i);
        }
    }
    return Fb;
}
