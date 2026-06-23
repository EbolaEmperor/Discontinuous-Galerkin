#include "NeoHookeanSolid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace euler_ale {
namespace {

double orient2dNeo(const Vector2d& a, const Vector2d& b, const Vector2d& c) {
    return (b.x() - a.x()) * (c.y() - a.y()) -
           (b.y() - a.y()) * (c.x() - a.x());
}

double minReferenceEdge(const Mesh& mesh) {
    double h = std::numeric_limits<double>::max();
    for (int e = 0; e < mesh.elem.rows(); ++e) {
        for (int k = 0; k < 3; ++k) {
            Vector2d a = mesh.node.row(mesh.elem(e, k)).transpose();
            Vector2d b = mesh.node.row(mesh.elem(e, (k + 1) % 3)).transpose();
            h = std::min(h, (b - a).norm());
        }
    }
    return h == std::numeric_limits<double>::max() ? 0.0 : h;
}

} // namespace

NeoHookeanMaterial NeoHookeanMaterial::fromSolidMaterial(const SolidMaterial& material) {
    NeoHookeanMaterial out;
    out.density = material.density;
    out.thickness = material.thickness;
    out.young = material.young;
    out.poisson = material.poisson;
    out.damping = material.damping;
    return out;
}

double NeoHookeanMaterial::mu() const {
    double nu = std::clamp(poisson, -0.95, 0.495);
    return young / std::max(1e-14, 2.0 * (1.0 + nu));
}

double NeoHookeanMaterial::lambda() const {
    double nu = std::clamp(poisson, -0.95, 0.495);
    return young * nu /
           std::max(1e-14, (1.0 + nu) * (1.0 - 2.0 * nu));
}

NeoHookeanSolidModel::NeoHookeanSolidModel(const NeoHookeanMaterial& material)
    : material_(material) {}

const NeoHookeanMaterial& NeoHookeanSolidModel::material() const {
    return material_;
}

MatrixXd NeoHookeanSolidModel::internalForces(const ElasticSolid2D& solid,
                                              double* strainEnergy) const {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    const MatrixXi& elem = solid.elements();
    MatrixXd force = MatrixXd::Zero(X.rows(), 2);
    double energy = 0.0;
    const double mu = material_.mu();
    const double lambda = material_.lambda();

    for (int e = 0; e < elem.rows(); ++e) {
        int id[3] = {elem(e, 0), elem(e, 1), elem(e, 2)};
        Vector2d Xp[3] = {
            X.row(id[0]).transpose(),
            X.row(id[1]).transpose(),
            X.row(id[2]).transpose()
        };
        Vector2d xp[3] = {
            x.row(id[0]).transpose(),
            x.row(id[1]).transpose(),
            x.row(id[2]).transpose()
        };
        double twiceA = orient2dNeo(Xp[0], Xp[1], Xp[2]);
        double area0 = 0.5 * std::abs(twiceA);
        if (area0 <= 1e-16) continue;

        double b[3] = {Xp[1].y() - Xp[2].y(), Xp[2].y() - Xp[0].y(),
                       Xp[0].y() - Xp[1].y()};
        double c[3] = {Xp[2].x() - Xp[1].x(), Xp[0].x() - Xp[2].x(),
                       Xp[1].x() - Xp[0].x()};
        Vector2d gradN[3];
        for (int i = 0; i < 3; ++i) {
            gradN[i] = Vector2d(b[i], c[i]) / (2.0 * area0);
        }

        Matrix2d F = Matrix2d::Zero();
        for (int i = 0; i < 3; ++i) {
            F += xp[i] * gradN[i].transpose();
        }
        double Jraw = F.determinant();
        double J = ((Jraw > 1e-10) && std::isfinite(Jraw)) ? Jraw : 1e-10;
        Matrix2d FinvT;
        FinvT << F(1, 1), -F(1, 0),
                 -F(0, 1), F(0, 0);
        FinvT /= J;
        double logJ = std::log(J);
        Matrix2d P = mu * (F - FinvT) + lambda * logJ * FinvT;

        double I1 = (F.transpose() * F).trace();
        double W = 0.5 * mu * (I1 - 2.0 - 2.0 * logJ) +
                   0.5 * lambda * logJ * logJ;
        energy += material_.thickness * area0 * W;

        for (int i = 0; i < 3; ++i) {
            Vector2d nodal = material_.thickness * area0 * (P * gradN[i]);
            force.row(id[i]) += nodal.transpose();
        }
    }
    if (strainEnergy) *strainEnergy = energy;
    return force;
}

double NeoHookeanSolidModel::strainEnergy(const ElasticSolid2D& solid) const {
    double energy = 0.0;
    internalForces(solid, &energy);
    return energy;
}

double NeoHookeanSolidModel::kineticEnergy(const ElasticSolid2D& solid) const {
    const MatrixXd& v = solid.velocities();
    const VectorXd& m = solid.lumpedMasses();
    double energy = 0.0;
    for (int i = 0; i < v.rows(); ++i) {
        energy += 0.5 * m(i) * v.row(i).squaredNorm();
    }
    return energy;
}

double NeoHookeanSolidModel::stableTimeStep(const ElasticSolid2D& solid, double cfl) const {
    Mesh mesh;
    mesh.node = solid.referenceNodes();
    mesh.elem = solid.elements();
    double h = minReferenceEdge(mesh);
    if (!(h > 0.0)) return 1e300;
    double wave = std::sqrt((material_.lambda() + 2.0 * material_.mu()) /
                            std::max(1e-14, material_.density));
    return cfl * h / std::max(wave, 1e-14);
}

void NeoHookeanSolidModel::advanceExplicit(ElasticSolid2D& solid, double dt) const {
    if (solid.numNodes() == 0) return;
    MatrixXd fint = internalForces(solid);
    MatrixXd x = solid.currentNodes();
    MatrixXd v = solid.velocities();
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& fext = solid.externalForces();
    const VectorXd& mass = solid.lumpedMasses();
    const VectorXi& fixed = solid.fixedMask();

    for (int i = 0; i < x.rows(); ++i) {
        if (fixed(i)) {
            x.row(i) = X.row(i);
            v.row(i).setZero();
            continue;
        }
        Vector2d fdamp = material_.damping * mass(i) * v.row(i).transpose();
        Vector2d a = (fext.row(i).transpose() - fint.row(i).transpose() - fdamp) /
                     std::max(mass(i), 1e-14);
        v.row(i) += (dt * a).transpose();
        x.row(i) += (dt * v.row(i).transpose()).transpose();
    }
    solid.setState(x, v);
}

} // namespace euler_ale
