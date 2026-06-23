#include "DG.h"

#include <cmath>
#include <functional>

namespace euler {

using namespace Eigen;

void makeRectMesh(Mesh& mesh, double x0, double x1, double y0, double y1, int nx, int ny) {
    int nnx = nx + 1, nny = ny + 1;
    mesh.node.resize(nnx * nny, 2);
    double dx = (x1 - x0) / nx, dy = (y1 - y0) / ny;
    auto nid = [&](int i, int j){ return j * nnx + i; };
    for (int j = 0; j < nny; ++j)
        for (int i = 0; i < nnx; ++i) {
            mesh.node(nid(i, j), 0) = x0 + i * dx;
            mesh.node(nid(i, j), 1) = y0 + j * dy;
        }
    mesh.elem.resize(2 * nx * ny, 3);
    int t = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int a = nid(i, j), b = nid(i + 1, j), c = nid(i + 1, j + 1), d = nid(i, j + 1);
            mesh.elem.row(t++) << a, b, c;
            mesh.elem.row(t++) << a, c, d;
        }
}

MatrixXd projectInitial(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                        const std::function<Vector4d(double, double)>& primField) {
    int NT = static_cast<int>(mesh.elem.rows()), locDof = fem.locDof;
    int nDof = static_cast<int>(elem2dof.maxCoeff()) + 1;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    MatrixXd Mref = MatrixXd::Zero(locDof, locDof);
    for (int q = 0; q < nq; ++q) { phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0); Mref.noalias() += w(q) * (phi[q].transpose() * phi[q]); }
    MatrixXd MrefInv = Mref.inverse();
    MatrixXd U = MatrixXd::Zero(nDof, 4);
    for (int t = 0; t < NT; ++t) {
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)), p2 = mesh.node.row(mesh.elem(t, 1)), p3 = mesh.node.row(mesh.elem(t, 2));
        Matrix<double, Dynamic, 4> be = Matrix<double, Dynamic, 4>::Zero(locDof, 4);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d Pp = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            Vector4d pr = primField(Pp.x(), Pp.y());
            Vector4d Uc = primToCons(pr(0), pr(1), pr(2), pr(3));
            be.noalias() += w(q) * (phi[q].transpose() * Uc.transpose());
        }
        Matrix<double, Dynamic, 4> Ue = MrefInv * be;
        for (int i = 0; i < locDof; ++i) U.row(elem2dof(t, i)) = Ue.row(i);
    }
    return U;
}

double l2Error(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
               const MatrixXd& U, int comp,
               const std::function<Vector4d(double, double)>& exactCons) {
    int NT = static_cast<int>(mesh.elem.rows()), locDof = fem.locDof;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
    double err = 0.0;
    for (int t = 0; t < NT; ++t) {
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0)), p2 = mesh.node.row(mesh.elem(t, 1)), p3 = mesh.node.row(mesh.elem(t, 2));
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d Pp = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            double uh = 0.0;
            for (int i = 0; i < locDof; ++i) uh += phi[q](i) * U(elem2dof(t, i), comp);
            double ue = exactCons(Pp.x(), Pp.y())(comp);
            err += w(q) * area * (uh - ue) * (uh - ue);
        }
    }
    return std::sqrt(err);
}

} // namespace euler
