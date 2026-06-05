#include "PostProcess.h"
#include "Quadrature.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXi;

namespace {
double cross2d(const Vector2d& a, const Vector2d& b) {
    return a(0) * b(1) - a(1) * b(0);
}
} // namespace

VectorXd solveLocalPoisson(FEM& femStar, const VecFEM& femSigma, FEM& femU,
                           const Mesh& mesh, const MatrixXi& elem2dofStar,
                           const MatrixXi& elem2dofSigma, const MatrixXi& elem2dofU,
                           const VectorXd& sigmah, const VectorXd& uh,
                           const std::function<double(const Vector2d&)>& f) {
    int NT = mesh.elem.rows();
    int locStar = femStar.locDof;
    int locSig = femSigma.locDof;
    int locU = femU.locDof;
    int nStar = elem2dofStar.maxCoeff() + 1;

    VectorXd ustar = VectorXd::Zero(nStar);

    MatrixXd quad2L;
    VectorXd w2;
    femStar.quad2d(quad2L, w2);
    int nq2 = static_cast<int>(w2.size());

    MatrixXd quad1L;
    VectorXd w1;
    femStar.quad1d(quad1L, w1);
    int nq1 = static_cast<int>(w1.size());

    int edgeLoc[3][2] = {{1, 2}, {2, 0}, {0, 1}}; // zero-based

    for (int t = 0; t < NT; ++t) {
        VectorXi idxS = elem2dofStar.row(t).transpose();
        VectorXi idxSig = elem2dofSigma.row(t).transpose();
        VectorXi idxU = elem2dofU.row(t).transpose();

        VectorXd sigma_t(locSig);
        for (int i = 0; i < locSig; ++i) sigma_t(i) = sigmah(idxSig(i));
        VectorXd uh_t(locU);
        for (int i = 0; i < locU; ++i) uh_t(i) = uh(idxU(i));

        MatrixXd Kloc = MatrixXd::Zero(locStar, locStar);
        VectorXd Floc = VectorXd::Zero(locStar);
        RowVectorXd m_row = RowVectorXd::Zero(locStar);
        double m_u = 0.0;

        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)).transpose();
        Vector2d p2 = mesh.node.row(mesh.elem(t, 1)).transpose();
        Vector2d p3 = mesh.node.row(mesh.elem(t, 2)).transpose();
        double area = femStar.area(t);
        const MatrixXd& DlamStar = femStar.Dlam[t];

        for (int q = 0; q < nq2; ++q) {
            Vector3d lam = quad2L.row(q).transpose();
            MatrixXd lamMat(1, 3);
            lamMat << lam.transpose();

            RowVectorXd phiS = femStar.computeBasisValue_all(lamMat).row(0);
            MatrixXd grad_phiS = DlamStar * femStar.computeBasisDlam_all(lam); // 2 x locStar

            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            double fq = f(pt);

            double wgt = w2(q) * area;
            Kloc.noalias() += wgt * (grad_phiS.transpose() * grad_phiS);
            Floc.noalias() += wgt * (phiS.transpose() * fq);

            m_row.noalias() += phiS * wgt;
            double uval = femU.computeBasisValue_all(lamMat).row(0).dot(uh_t);
            m_u += uval * wgt;
        }

        MatrixXd vtx(3, 2);
        vtx.row(0) = p1.transpose();
        vtx.row(1) = p2.transpose();
        vtx.row(2) = p3.transpose();

        for (int ee = 0; ee < 3; ++ee) {
            int a = edgeLoc[ee][0];
            int b = edgeLoc[ee][1];
            Vector2d uxy = vtx.row(a).transpose();
            Vector2d vxy = vtx.row(b).transpose();
            Vector2d evec = vxy - uxy;
            double he = evec.norm();
            Vector2d n = Vector2d(evec(1), -evec(0)) / he;

            int wid = 3 - (a + b);
            Vector2d wpt = vtx.row(wid).transpose();
            if (cross2d(evec, wpt - uxy) < 0) n = -n;

            for (int q = 0; q < nq1; ++q) {
                Vector3d lamEdge = Vector3d::Zero();
                lamEdge(a) = quad1L(q, 0);
                lamEdge(b) = quad1L(q, 1);

                MatrixXd sigBasis = femSigma.computeBasisValue_all(lamEdge); // 2 x locSig
                RowVectorXd nsig_row = n.transpose() * sigBasis;            // 1 x locSig
                double nsig_val = nsig_row.dot(sigma_t);

                MatrixXd lamMat(1, 3);
                lamMat << lamEdge.transpose();
                RowVectorXd phiS = femStar.computeBasisValue_all(lamMat).row(0);

                Floc.noalias() += phiS.transpose() * (nsig_val * w1(q) * he);
            }
        }

        MatrixXd KKT(locStar + 1, locStar + 1);
        KKT.setZero();
        KKT.topLeftCorner(locStar, locStar) = Kloc;
        KKT.block(locStar, 0, 1, locStar) = m_row;
        KKT.block(0, locStar, locStar, 1) = m_row.transpose();

        VectorXd rhs(locStar + 1);
        rhs.head(locStar) = Floc;
        rhs(locStar) = m_u;

        VectorXd sol = KKT.fullPivHouseholderQr().solve(rhs);
        for (int i = 0; i < locStar; ++i) {
            ustar(idxS(i)) = sol(i);
        }
    }

    return ustar;
}
