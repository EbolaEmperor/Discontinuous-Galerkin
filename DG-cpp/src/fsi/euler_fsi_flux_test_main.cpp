#include "EulerALE1D.h"
#include "EulerALE2D.h"
#include "EulerDG.h"
#include "FEM.h"
#include "FSIDiagnostics.h"
#include "Mesh.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace euler_fsi;

namespace {

bool near(double a, double b, double tol, const char* label) {
    double e = std::abs(a - b);
    if (e > tol) {
        std::cerr << label << " mismatch: got " << a << " expected " << b
                  << " abs_err=" << e << "\n";
        return false;
    }
    return true;
}

bool nearState(const State1D& a, const State1D& b, double tol, const char* label) {
    bool ok = true;
    ok = near(a.rho, b.rho, tol, label) && ok;
    ok = near(a.mom, b.mom, tol, label) && ok;
    ok = near(a.ene, b.ene, tol, label) && ok;
    return ok;
}

bool nearState(const State2D& a, const State2D& b, double tol, const char* label) {
    bool ok = true;
    ok = near(a.rho, b.rho, tol, label) && ok;
    ok = near(a.mx, b.mx, tol, label) && ok;
    ok = near(a.my, b.my, tol, label) && ok;
    ok = near(a.ene, b.ene, tol, label) && ok;
    return ok;
}

bool nearVector(const Eigen::Vector4d& a, const Eigen::Vector4d& b,
                double tol, const char* label) {
    bool ok = true;
    for (int i = 0; i < 4; ++i) {
        std::string component = std::string(label) + "[" + std::to_string(i) + "]";
        ok = near(a(i), b(i), tol, component.c_str()) && ok;
    }
    return ok;
}

} // namespace

int main() {
    bool ok = true;

    State1D u1 = primToCons(Prim1D{1.2, 0.7, 2.3});
    double w1 = 0.15;
    ok = nearState(aleRusanovFlux(u1, u1, w1), flux(u1) - w1 * u1, 1e-13,
                   "1D ALE equal-state flux") && ok;
    State1D wall1 = movingWallFlux(u1, -0.22);
    double p1 = pressure(u1);
    ok = near(wall1.rho, 0.0, 1e-14, "1D moving wall mass") && ok;
    ok = near(wall1.mom, p1, 1e-14, "1D moving wall momentum") && ok;
    ok = near(wall1.ene, p1 * -0.22, 1e-14, "1D moving wall energy") && ok;

    State2D u2 = primToCons(Prim2D{0.9, 0.4, -0.2, 1.7});
    double wx = -0.08, wy = 0.11;
    ok = nearState(aleRusanovFluxX(u2, u2, wx), fluxX(u2) - wx * u2, 1e-13,
                   "2D ALE-X equal-state flux") && ok;
    ok = nearState(aleRusanovFluxY(u2, u2, wy), fluxY(u2) - wy * u2, 1e-13,
                   "2D ALE-Y equal-state flux") && ok;
    double p2 = pressure(u2);
    ok = nearState(movingWallFluxX(u2, 0.31), State2D{0.0, p2, 0.0, p2 * 0.31},
                   1e-14, "2D moving X-wall flux") && ok;
    ok = nearState(movingWallFluxY(u2, -0.27), State2D{0.0, 0.0, p2, p2 * -0.27},
                   1e-14, "2D moving Y-wall flux") && ok;

    Eigen::Vector4d ug = euler::primToCons(1.1, 0.3, -0.4, 2.2);
    ok = nearVector(euler::aleRusanov(ug, ug, 1.0, 0.0, wx),
                    euler::normalFlux(ug, 1.0, 0.0) - wx * ug,
                    1e-13, "DG ALE x equal-state flux") && ok;
    ok = nearVector(euler::aleRusanov(ug, ug, 0.0, 1.0, wy),
                    euler::normalFlux(ug, 0.0, 1.0) - wy * ug,
                    1e-13, "DG ALE y equal-state flux") && ok;
    double nx = 0.6;
    double ny = 0.8;
    double wallU = 0.25;
    double wallV = -0.15;
    double pg = euler::pressure(ug);
    Eigen::Vector4d expectedWall(0.0, pg * nx, pg * ny, pg * (wallU * nx + wallV * ny));
    ok = nearVector(euler::movingWallFlux(ug, nx, ny, wallU, wallV),
                    expectedWall, 1e-14, "DG moving-wall ALE flux") && ok;
    auto wallFlux = euler::movingWallBoundaryFlux(
        [wallU, wallV](double, double, double, int) {
            return std::array<double, 2>{wallU, wallV};
        });
    ok = nearVector(wallFlux(0.2, 0.3, 0.4, ug, nx, ny, 7),
                    expectedWall, 1e-14, "DG moving-wall boundary callback") && ok;
    Eigen::Vector4d ghost = euler::movingSlipWallExterior(ug, nx, ny, wallU, wallV);
    double relIn = ((ug(1) / ug(0)) - wallU) * nx + ((ug(2) / ug(0)) - wallV) * ny;
    double relGhost = ((ghost(1) / ghost(0)) - wallU) * nx + ((ghost(2) / ghost(0)) - wallV) * ny;
    ok = near(relGhost, -relIn, 1e-14, "DG moving-wall ghost relative normal") && ok;
    ok = near(euler::pressure(ghost), pg, 1e-14, "DG moving-wall ghost pressure") && ok;

    {
        Mesh mesh;
        euler::makeRectMesh(mesh, 0.0, 1.0, 0.0, 1.0, 2, 2);
        FEM fem(1, mesh, false);
        Eigen::MatrixXi elem2dof;
        int nDof = 0;
        fem.getDOF(mesh, elem2dof, nDof);
        Eigen::MatrixXi edge, edge2side;
        mesh.getEdge2Side(edge, edge2side);
        Eigen::VectorXi tag = Eigen::VectorXi::Zero(edge.rows());
        euler::EulerConfig cfg;
        cfg.use_av = false;
        cfg.use_positivity = false;
        cfg.use_hllc = false;
        euler::EulerDG dg(fem, mesh, elem2dof, edge, edge2side, tag, cfg);
        Eigen::Vector4d primitive(1.0, 0.25, -0.1, 1.0);
        dg.setState(euler::projectInitial(fem, mesh, elem2dof,
            [primitive](double, double) { return primitive; }));
        Eigen::MatrixXd U0 = dg.state();
        euler::BoundaryFluxFn transmissiveFlux(
            [](double, double, double, const Eigen::Vector4d& Uin,
               double nx, double ny, int) {
                return euler::normalFlux(Uin, nx, ny);
            });
        Eigen::MatrixXd R;
        dg.inviscidResidualWithBoundaryFlux(dg.state(), 0.0, transmissiveFlux, R);
        ok = near(R.norm(), 0.0, 1e-11, "DG direct boundary-flux residual constant state") && ok;
        ok = dg.stepWithBoundaryFlux(1e-3, 1e-3, transmissiveFlux) && ok;
        ok = near((dg.state() - U0).norm(), 0.0, 1e-11,
                  "DG direct boundary-flux step constant state") && ok;
    }

    SpringPiston piston;
    piston.x = 1.0;
    piston.v = 0.0;
    piston.mass = 2.0;
    piston.stiffness = 50.0;
    piston.damping = 4.0;
    piston.xRest = 1.0;
    piston.area = 2.0;
    piston.externalPressure = 1.0;
    FSIEnergyBudget budget;
    budget.reset(10.0, piston);
    double xOld = piston.x, vOld = piston.v;
    piston.x = 1.1;
    piston.v = 0.2;
    FSIEnergySnapshot snap = budget.advance(0.1, 9.4, piston, 3.0, xOld, vOld);
    ok = near(snap.gasWallWork, 0.6, 1e-14, "FSI gas wall work") && ok;
    ok = near(snap.dampingLoss, 0.008, 1e-14, "FSI damping loss") && ok;
    ok = near(snap.gasWorkDrift, 0.0, 1e-14, "FSI gas work drift") && ok;
    double expectedCoupled = 9.4 + piston.kineticEnergy() + piston.springEnergy()
                           + piston.externalPressurePotential() + snap.dampingLoss;
    double initialCoupled = 10.0 + 2.0;
    ok = near(snap.coupledEnergyDrift, expectedCoupled - initialCoupled,
              1e-14, "FSI coupled energy drift") && ok;

    std::cout << "Verification " << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 2;
}
