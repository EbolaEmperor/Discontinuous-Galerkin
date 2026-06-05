#include "AssemblyHDG.h"
#include "Quadrature.h"
#include <vector>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXi;

namespace {
double cross2d(const Vector2d& a, const Vector2d& b) {
    return a(0) * b(1) - a(1) * b(0);
}
} // namespace

SparseMatrix<double> assembleMass(VecFEM& fem, const Mesh& mesh, const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows();
    int locDof = fem.locDof;
    int nDof = elem2dof.maxCoeff() + 1;

    std::vector<Triplet<double>> trip;
    trip.reserve(NT * locDof * locDof);

    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());

    // Basis on reference element is independent of tid
    std::vector<MatrixXd> phi_q(nq);
    for (int q = 0; q < nq; ++q) {
        Vector3d lam = quadL.row(q).transpose();
        phi_q[q] = fem.computeBasisValue_all(lam); // 2 x locDof
    }

    MatrixXd At(locDof, locDof);
    for (int t = 0; t < NT; ++t) {
        At.setZero();
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) {
            double weight = w(q) * area;
            At.noalias() += weight * (phi_q[q].transpose() * phi_q[q]);
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

SparseMatrix<double> assembleDivMass(VecFEM& femSigma, FEM& femU, const Mesh& mesh,
                                     const MatrixXi& elem2dofSigma, const MatrixXi& elem2dofU) {
    int NT = mesh.elem.rows();
    int locSig = femSigma.locDof;
    int locU = femU.locDof;
    int nSig = elem2dofSigma.maxCoeff() + 1;
    int nU = elem2dofU.maxCoeff() + 1;

    std::vector<Triplet<double>> trip;
    trip.reserve(NT * locSig * locU);

    MatrixXd quadL;
    VectorXd w;
    Quadrature::quadpts(femSigma.ord + femU.ord - 1, quadL, w);
    int nq = static_cast<int>(w.size());

    // Precompute scalar basis of u on reference element
    std::vector<RowVectorXd> phi_u_q(nq);
    for (int q = 0; q < nq; ++q) {
        MatrixXd lamMat(1, 3);
        lamMat << quadL.row(q);
        phi_u_q[q] = femU.computeBasisValue_all(lamMat).row(0);
    }

    MatrixXd Bt(locSig, locU);
    for (int t = 0; t < NT; ++t) {
        Bt.setZero();
        double area = femSigma.area(t);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            RowVectorXd psi_div = femSigma.computeBasisDiv_all(t, lam);
            double weight = w(q) * area;
            Bt.noalias() += weight * (psi_div.transpose() * phi_u_q[q]);
        }

        for (int i = 0; i < locSig; ++i) {
            int gi = elem2dofSigma(t, i);
            for (int j = 0; j < locU; ++j) {
                trip.emplace_back(gi, elem2dofU(t, j), Bt(i, j));
            }
        }
    }

    SparseMatrix<double> B(nSig, nU);
    B.setFromTriplets(trip.begin(), trip.end());
    return B;
}

SparseMatrix<double> assembleGradMass(FEM& femU, VecFEM& femSigma, const Mesh& mesh,
                                      const MatrixXi& elem2dofU, const MatrixXi& elem2dofSigma) {
    int NT = mesh.elem.rows();
    int locSig = femSigma.locDof;
    int locU = femU.locDof;
    int nSig = elem2dofSigma.maxCoeff() + 1;
    int nU = elem2dofU.maxCoeff() + 1;

    std::vector<Triplet<double>> trip;
    trip.reserve(NT * locU * locSig);

    MatrixXd quadL;
    VectorXd w;
    Quadrature::quadpts(femSigma.ord + femU.ord - 1, quadL, w);
    int nq = static_cast<int>(w.size());

    // Precompute reference data for u
    std::vector<MatrixXd> grad_ref_u(nq);
    for (int q = 0; q < nq; ++q) {
        Vector3d lam = quadL.row(q).transpose();
        grad_ref_u[q] = femU.computeBasisDlam_all(lam); // 3 x locU
    }

    MatrixXd Bt(locU, locSig);
    for (int t = 0; t < NT; ++t) {
        Bt.setZero();
        double area = femU.area(t);
        const MatrixXd& Dlam = femU.Dlam[t];

        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            MatrixXd grad_phi = Dlam * grad_ref_u[q];              // 2 x locU
            MatrixXd phi_sigma = femSigma.computeBasisValue_all(lam); // 2 x locSig
            double weight = w(q) * area;
            Bt.noalias() += weight * (grad_phi.transpose() * phi_sigma);
        }

        for (int i = 0; i < locU; ++i) {
            int gi = elem2dofU(t, i);
            for (int j = 0; j < locSig; ++j) {
                trip.emplace_back(gi, elem2dofSigma(t, j), Bt(i, j));
            }
        }
    }

    SparseMatrix<double> B(nU, nSig);
    B.setFromTriplets(trip.begin(), trip.end());
    return B;
}

SparseMatrix<double> assembleIP_HDG(VecFEM& femSigma, FEM& femU, const Mesh& mesh,
                                    const MatrixXi& elem2dofSigma, const MatrixXi& elem2dofU,
                                    const MatrixXi& edge, const MatrixXi& edge2side,
                                    double alpha) {
    int NE = edge.rows();
    int locSig = femSigma.locDof;
    int locU = femU.locDof;
    int nSig = elem2dofSigma.maxCoeff() + 1;
    int nU = elem2dofU.maxCoeff() + 1;
    int nTot = nSig + nU;

    std::vector<Triplet<double>> trip;
    trip.reserve(NE * (4 * locSig * locSig + 4 * locSig * locU + 4 * locU * locU));

    MatrixXd quad1d;
    VectorXd w1d;
    Quadrature::quadpts1_my(femSigma.ord + femU.ord, quad1d, w1d);
    int nq = static_cast<int>(w1d.size());

    for (int e = 0; e < NE; ++e) {
        int t1 = edge2side(e, 0);
        int t2 = edge2side(e, 1);
        if (t1 == -1 && t2 == -1) continue;

        int n1 = edge(e, 0);
        int n2 = edge(e, 1);
        Vector2d p1 = mesh.node.row(n1).transpose();
        Vector2d p2 = mesh.node.row(n2).transpose();
        Vector2d evec = p2 - p1;
        double he = evec.norm();
        Vector2d nvec = Vector2d(evec(1), -evec(0)) / he;

        if (t1 == -1 || t2 == -1) {
            int tid = (t1 != -1) ? t1 : t2;
            Eigen::RowVectorXi tri = mesh.elem.row(tid);

            int locA = -1, locB = -1;
            for (int k = 0; k < 3; ++k) {
                if (tri(k) == n1) locA = k;
                if (tri(k) == n2) locB = k;
            }
            int wid = 3 - (locA + locB);
            Vector2d wpt = mesh.node.row(tri(wid)).transpose();
            if (cross2d(evec, wpt - p1) < 0) nvec = -nvec;

            MatrixXd Pt21 = MatrixXd::Zero(locU, locSig);
            MatrixXd Pt22 = MatrixXd::Zero(locU, locU);

            for (int q = 0; q < nq; ++q) {
                double lam1 = quad1d(q, 0);
                double lam2 = quad1d(q, 1);
                Vector3d lam = Vector3d::Zero();
                lam(locA) = lam1;
                lam(locB) = lam2;

                MatrixXd sigma_in = femSigma.computeBasisValue_all(lam); // 2 x locSig
                RowVectorXd nsigma = nvec.transpose() * sigma_in;        // 1 x locSig
                MatrixXd lamMat(1, 3);
                lamMat << lam.transpose();
                RowVectorXd u_in = femU.computeBasisValue_all(lamMat).row(0);   // 1 x locU

                double wgt = w1d(q) * he;
                Pt21.noalias() -= wgt * (u_in.transpose() * nsigma);
                Pt22.noalias() += wgt * alpha * (u_in.transpose() * u_in);
            }

            VectorXi idxSig = elem2dofSigma.row(tid).transpose();
            VectorXi idxU = elem2dofU.row(tid).transpose();
            for (int i = 0; i < locU; ++i) {
                int gi = idxU(i) + nSig;
                for (int j = 0; j < locSig; ++j) {
                    trip.emplace_back(gi, idxSig(j), Pt21(i, j));
                }
                for (int j = 0; j < locU; ++j) {
                    trip.emplace_back(gi, idxU(j) + nSig, Pt22(i, j));
                }
            }
        } else {
            Eigen::RowVectorXi tri1 = mesh.elem.row(t1);
            Eigen::RowVectorXi tri2 = mesh.elem.row(t2);

            int idx1a = -1, idx1b = -1, idx2a = -1, idx2b = -1;
            for (int k = 0; k < 3; ++k) {
                if (tri1(k) == n1) idx1a = k;
                if (tri1(k) == n2) idx1b = k;
                if (tri2(k) == n1) idx2a = k;
                if (tri2(k) == n2) idx2b = k;
            }

            MatrixXd Pt11 = MatrixXd::Zero(2 * locSig, 2 * locSig);
            MatrixXd Pt12 = MatrixXd::Zero(2 * locSig, 2 * locU);
            MatrixXd Pt21 = MatrixXd::Zero(2 * locU, 2 * locSig);
            MatrixXd Pt22 = MatrixXd::Zero(2 * locU, 2 * locU);

            for (int q = 0; q < nq; ++q) {
                double lam1 = quad1d(q, 0);
                double lam2 = quad1d(q, 1);
                Vector3d lamL = Vector3d::Zero();
                Vector3d lamR = Vector3d::Zero();
                lamL(idx1a) = lam1; lamL(idx1b) = lam2;
                lamR(idx2a) = lam1; lamR(idx2b) = lam2;

                MatrixXd sigma_left = MatrixXd::Zero(2, 2 * locSig);
                MatrixXd sigma_right = MatrixXd::Zero(2, 2 * locSig);
                sigma_left.block(0, 0, 2, locSig) = femSigma.computeBasisValue_all(lamL);
                sigma_right.block(0, locSig, 2, locSig) = femSigma.computeBasisValue_all(lamR);

                RowVectorXd sigma_jump = nvec.transpose() * (sigma_left - sigma_right); // 1 x (2*locSig)
                MatrixXd sigma_mean = 0.5 * (sigma_left + sigma_right);                 // 2 x (2*locSig)

                RowVectorXd u_left = RowVectorXd::Zero(2 * locU);
                RowVectorXd u_right = RowVectorXd::Zero(2 * locU);
                MatrixXd lamMatL(1, 3), lamMatR(1, 3);
                lamMatL << lamL.transpose();
                lamMatR << lamR.transpose();
                u_left.segment(0, locU) = femU.computeBasisValue_all(lamMatL).row(0);
                u_right.segment(locU, locU) = femU.computeBasisValue_all(lamMatR).row(0);

                MatrixXd u_jump = nvec * (u_left - u_right); // 2 x (2*locU)
                RowVectorXd u_mean = 0.5 * (u_left + u_right);           // 1 x (2*locU)

                double wgt = w1d(q) * he;
                Pt11.noalias() += (wgt / (2.0 * alpha)) * (sigma_jump.transpose() * sigma_jump);
                Pt12.noalias() -= wgt * (sigma_jump.transpose() * u_mean);
                Pt21.noalias() -= wgt * (u_jump.transpose() * sigma_mean);
                Pt22.noalias() += (wgt * alpha / 2.0) * (u_jump.transpose() * u_jump);
            }

            VectorXi idxSig(2 * locSig);
            idxSig << elem2dofSigma.row(t1).transpose(), elem2dofSigma.row(t2).transpose();
            VectorXi idxU(2 * locU);
            idxU << elem2dofU.row(t1).transpose(), elem2dofU.row(t2).transpose();

            for (int i = 0; i < 2 * locSig; ++i) {
                int gi = idxSig(i);
                for (int j = 0; j < 2 * locSig; ++j) {
                    trip.emplace_back(gi, idxSig(j), Pt11(i, j));
                }
                for (int j = 0; j < 2 * locU; ++j) {
                    trip.emplace_back(gi, idxU(j) + nSig, Pt12(i, j));
                }
            }
            for (int i = 0; i < 2 * locU; ++i) {
                int gi = idxU(i) + nSig;
                for (int j = 0; j < 2 * locSig; ++j) {
                    trip.emplace_back(gi, idxSig(j), Pt21(i, j));
                }
                for (int j = 0; j < 2 * locU; ++j) {
                    trip.emplace_back(gi, idxU(j) + nSig, Pt22(i, j));
                }
            }
        }
    }

    SparseMatrix<double> P(nTot, nTot);
    P.setFromTriplets(trip.begin(), trip.end());
    return P;
}

VectorXd assembleWeakBDC(VecFEM& femSigma, FEM& femU, const Mesh& mesh,
                         const MatrixXi& elem2dofSigma, const MatrixXi& elem2dofU,
                         const MatrixXi& edge, const MatrixXi& edge2side,
                         double alpha, const ExactSolution& sol) {
    int NE = edge.rows();
    int locSig = femSigma.locDof;
    int locU = femU.locDof;
    int nSig = elem2dofSigma.maxCoeff() + 1;
    int nU = elem2dofU.maxCoeff() + 1;
    int nTot = nSig + nU;

    VectorXd F1 = VectorXd::Zero(nSig);
    VectorXd F2 = VectorXd::Zero(nU);

    MatrixXd quad1d;
    VectorXd w1d;
    Quadrature::quadpts1_my(femSigma.ord + femU.ord, quad1d, w1d);
    int nq = static_cast<int>(w1d.size());

    for (int e = 0; e < NE; ++e) {
        int t1 = edge2side(e, 0);
        int t2 = edge2side(e, 1);
        if (t1 != -1 && t2 != -1) continue;

        int tid = (t1 != -1) ? t1 : t2;
        Eigen::RowVectorXi tri = mesh.elem.row(tid);
        VectorXi idxSig = elem2dofSigma.row(tid);
        VectorXi idxU = elem2dofU.row(tid);

        int nA = edge(e, 0);
        int nB = edge(e, 1);
        Vector2d u = mesh.node.row(nA).transpose();
        Vector2d v = mesh.node.row(nB).transpose();
        Vector2d evec = v - u;
        double he = evec.norm();
        Vector2d nvec = Vector2d(evec(1), -evec(0)) / he;

        int locA = -1, locB = -1;
        for (int k = 0; k < 3; ++k) {
            if (tri(k) == nA) locA = k;
            if (tri(k) == nB) locB = k;
        }
        int wid = 3 - (locA + locB);
        Vector2d wpt = mesh.node.row(tri(wid)).transpose();
        if (cross2d(evec, wpt - u) < 0) nvec = -nvec;

        for (int q = 0; q < nq; ++q) {
            double lam1 = quad1d(q, 0);
            double lam2 = quad1d(q, 1);
            Vector3d lam = Vector3d::Zero();
            lam(locA) = lam1; lam(locB) = lam2;

            MatrixXd sigma_in = femSigma.computeBasisValue_all(lam);
            RowVectorXd nsigma = nvec.transpose() * sigma_in; // 1 x locSig

            MatrixXd lamMat(1, 3);
            lamMat << lam.transpose();
            RowVectorXd u_in = femU.computeBasisValue_all(lamMat).row(0); // 1 x locU

            Vector2d pt = lam(0) * mesh.node.row(tri(0)).transpose()
                        + lam(1) * mesh.node.row(tri(1)).transpose()
                        + lam(2) * mesh.node.row(tri(2)).transpose();
            MatrixXd ptMat(1, 2);
            ptMat << pt.transpose();
            double gval = sol.u_exact(ptMat)(0);

            double wgt = w1d(q) * he;
            for (int i = 0; i < locSig; ++i) {
                F1(idxSig(i)) += wgt * nsigma(i) * gval;
            }
            for (int i = 0; i < locU; ++i) {
                F2(idxU(i)) += wgt * alpha * u_in(i) * gval;
            }
        }
    }

    VectorXd F(nTot);
    F.head(nSig) = F1;
    F.tail(nU) = F2;
    return F;
}

double l2ErrorScalar(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                     const VectorXd& uh, const ExactSolution& sol) {
    int NT = mesh.elem.rows();
    int locDof = fem.locDof;

    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());

    double err = 0.0;
    for (int t = 0; t < NT; ++t) {
        VectorXd uh_t(locDof);
        for (int i = 0; i < locDof; ++i) uh_t(i) = uh(elem2dof(t, i));

        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)).transpose();
        Vector2d p2 = mesh.node.row(mesh.elem(t, 1)).transpose();
        Vector2d p3 = mesh.node.row(mesh.elem(t, 2)).transpose();
        double area = fem.area(t);

        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            MatrixXd lamMat(1, 3);
            lamMat << lam.transpose();
            RowVectorXd phi = fem.computeBasisValue_all(lamMat).row(0);
            double uh_val = phi * uh_t;

            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            MatrixXd ptMat(1, 2);
            ptMat << pt.transpose();
            double u_ex = sol.u_exact(ptMat)(0);

            double diff = uh_val - u_ex;
            err += w(q) * area * diff * diff;
        }
    }

    return std::sqrt(err);
}

double l2ErrorVector(VecFEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                     const VectorXd& sigmah, const ExactSolution& sol) {
    int NT = mesh.elem.rows();
    int locDof = fem.locDof;

    MatrixXd quadL;
    VectorXd w;
    fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());

    double err = 0.0;
    for (int t = 0; t < NT; ++t) {
        VectorXd sigma_t(locDof);
        for (int i = 0; i < locDof; ++i) sigma_t(i) = sigmah(elem2dof(t, i));

        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)).transpose();
        Vector2d p2 = mesh.node.row(mesh.elem(t, 1)).transpose();
        Vector2d p3 = mesh.node.row(mesh.elem(t, 2)).transpose();
        double area = fem.area(t);

        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            MatrixXd phi_sigma = fem.computeBasisValue_all(lam); // 2 x locDof
            Vector2d sigma_val = phi_sigma * sigma_t;

            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            MatrixXd ptMat(1, 2);
            ptMat << pt.transpose();
            MatrixXd gradMat = sol.grad_u_exact(ptMat);
            Vector2d grad_ex = gradMat.row(0).transpose();

            Vector2d diff = sigma_val - grad_ex;
            err += w(q) * area * diff.squaredNorm();
        }
    }

    return std::sqrt(err);
}
