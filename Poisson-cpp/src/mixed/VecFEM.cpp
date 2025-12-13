#include "VecFEM.h"
#include "Quadrature.h"

VecFEM::VecFEM(int order, const Mesh& mesh) : ord(order) {
    scalarLocDof = (ord + 1) * (ord + 2) / 2;
    locDof = 2 * scalarLocDof;
    gradbasis_my(mesh, Dlam, area);
}

void VecFEM::getDOF(const Mesh& mesh, Eigen::MatrixXi& elem2dof, int& nDof) const {
    int NT = mesh.elem.rows();
    nDof = NT * locDof;
    elem2dof.resize(NT, locDof);

    int idx = 0;
    for (int t = 0; t < NT; ++t) {
        for (int i = 0; i < locDof; ++i) {
            elem2dof(t, i) = idx++;
        }
    }
}

Eigen::MatrixXd VecFEM::computeBasisValue_all(const Eigen::Vector3d& lam) const {
    Eigen::MatrixXd lamMat(1, 3);
    lamMat << lam.transpose();
    Eigen::MatrixXd span = polyBasisHomo3D(ord, lamMat); // 1 x scalarLocDof

    Eigen::MatrixXd val = Eigen::MatrixXd::Zero(2, locDof);
    val.block(0, 0, 1, scalarLocDof) = span;
    val.block(1, scalarLocDof, 1, scalarLocDof) = span;
    return val;
}

Eigen::RowVectorXd VecFEM::computeBasisDiv_all(int tid, const Eigen::Vector3d& lam) const {
    Eigen::MatrixXd gradRef = polyBasisHomoGrad3D(ord, lam); // 3 x scalarLocDof
    Eigen::MatrixXd gradXY = Dlam[tid] * gradRef;            // 2 x scalarLocDof

    Eigen::RowVectorXd div(locDof);
    div << gradXY.row(0), gradXY.row(1);
    return div;
}

void VecFEM::quad2d(Eigen::MatrixXd& quadL, Eigen::VectorXd& w) const {
    Quadrature::quadpts2_my(2 * (ord + 1), quadL, w);
}

void VecFEM::quad1d(Eigen::MatrixXd& quadL, Eigen::VectorXd& w) const {
    Quadrature::quadpts1_my(2 * (ord + 1), quadL, w);
}
