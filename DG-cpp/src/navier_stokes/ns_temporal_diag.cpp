// TEMPORARY diagnostic: characterise the *temporal* order of the DG-NS splitting
// integrator on the Taylor-Green vortex, decoupled from the spatial error by
// measuring against a fine-dt reference on the SAME mesh (so the spatial error
// cancels and only the time-splitting error remains):
//   error(dt) = || u_h(dt) - u_h(dt_ref) ||_{L2},   dt_ref << all dt.
// Also isolates the pressure boundary condition (high-order Neumann vs exact
// Dirichlet) and reports || div u_h ||_{L2} vs dt (the splitting boundary layer).
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "FEM.h"
#include "Mesh.h"
#include "NavierStokes.h"

using namespace Eigen;
using namespace std;
using namespace ns;

struct Exact {
    double nu;
    double e2(double t) const { return std::exp(-2.0 * nu * t); }
    double e4(double t) const { return std::exp(-4.0 * nu * t); }
    Vector2d vel(double x, double y, double t) const {
        double e = e2(t);
        return {-std::cos(x) * std::sin(y) * e, std::sin(x) * std::cos(y) * e};
    }
    Vector2d acc(double x, double y, double t) const { return -2.0 * nu * vel(x, y, t); }
    double pres(double x, double y, double t) const {
        return -0.25 * (std::cos(2 * x) + std::cos(2 * y)) * e4(t);
    }
    Vector2d presGrad(double x, double y, double t) const {  // exact grad p
        double e = e4(t);
        return {0.5 * std::sin(2 * x) * e, 0.5 * std::sin(2 * y) * e};
    }
};

// BC mode:  0 = shipped (p-Dirichlet on right wall, high-order Neumann elsewhere)
//           1 = exact pressure Dirichlet on ALL walls (removes the Neumann splitting BC)
static BCData makeBC(const Mesh& mesh, const MatrixXi& edge, const MatrixXi& edge2side,
                     const Exact& ex, int pmode, bool exactNeu = false) {
    int NE = edge.rows();
    BCData bc;
    bc.bcU = VectorXi::Zero(NE); bc.bcV = VectorXi::Zero(NE); bc.bcP = VectorXi::Zero(NE);
    const double tol = 1e-9;
    for (int e = 0; e < NE; ++e) {
        if (edge2side(e, 0) != -1 && edge2side(e, 1) != -1) continue;
        Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1)));
        bc.bcU(e) = BCN_DIRICHLET;
        bc.bcV(e) = BCN_DIRICHLET;
        if (pmode == 1)                       bc.bcP(e) = BCN_DIRICHLET;   // exact p everywhere
        else if (std::abs(m.x() - 1.0) < tol) bc.bcP(e) = BCN_DIRICHLET;
        else                                  bc.bcP(e) = BCN_NEUMANN_HO;
    }
    bc.velDir = [ex](double x, double y, double t) { return ex.vel(x, y, t); };
    bc.velAcc = [ex](double x, double y, double t) { return ex.acc(x, y, t); };
    bc.presDir = [ex](double x, double y, double t) { return ex.pres(x, y, t); };
    if (exactNeu) bc.presGrad = [ex](double x, double y, double t) { return ex.presGrad(x, y, t); };
    return bc;
}

static void projectExact(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof, const Exact& ex,
                         double t, VectorXd& u0, VectorXd& v0, VectorXd& p0) {
    int nDof = elem2dof.maxCoeff() + 1, NT = mesh.elem.rows(), locDof = fem.locDof;
    SparseMatrix<double> M = assembleScalarMassDG(fem, mesh, elem2dof);
    SimplicialLDLT<SparseMatrix<double>> luM(M);
    VectorXd bu = VectorXd::Zero(nDof), bv = VectorXd::Zero(nDof), bp = VectorXd::Zero(nDof);
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
    for (int e = 0; e < NT; ++e) {
        Vector2d p1 = mesh.node.row(mesh.elem(e, 0));
        Vector2d p2 = mesh.node.row(mesh.elem(e, 1));
        Vector2d p3 = mesh.node.row(mesh.elem(e, 2));
        double area = fem.area(e);
        VectorXd lu = VectorXd::Zero(locDof), lv = lu, lp = lu;
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d P = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            Vector2d ve = ex.vel(P.x(), P.y(), t);
            double pe = ex.pres(P.x(), P.y(), t);
            double wa = w(q) * area;
            lu.noalias() += wa * ve.x() * phi[q].transpose();
            lv.noalias() += wa * ve.y() * phi[q].transpose();
            lp.noalias() += wa * pe * phi[q].transpose();
        }
        for (int i = 0; i < locDof; ++i) { bu(elem2dof(e, i)) += lu(i); bv(elem2dof(e, i)) += lv(i); bp(elem2dof(e, i)) += lp(i); }
    }
    u0 = luM.solve(bu); v0 = luM.solve(bv); p0 = luM.solve(bp);
}

static double velL2Diff(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
                        const VectorXd& ua, const VectorXd& va,
                        const VectorXd& ub, const VectorXd& vb) {
    int NT = mesh.elem.rows(), locDof = fem.locDof;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
    double err = 0.0;
    VectorXd due(locDof), dve(locDof);
    for (int e = 0; e < NT; ++e) {
        for (int i = 0; i < locDof; ++i) { due(i) = ua(elem2dof(e, i)) - ub(elem2dof(e, i));
                                           dve(i) = va(elem2dof(e, i)) - vb(elem2dof(e, i)); }
        double area = fem.area(e);
        for (int q = 0; q < nq; ++q) {
            double du = phi[q].dot(due), dv = phi[q].dot(dve);
            err += w(q) * area * (du * du + dv * dv);
        }
    }
    return std::sqrt(err);
}

struct Ctx {
    FEM* fem; Mesh* mesh; const MatrixXi* elem2dof; const MatrixXi* edge; const MatrixXi* edge2side;
    SparseMatrix<double> M, Gx, Gy; SimplicialLDLT<SparseMatrix<double>> luM;
};

// L2 norm of the discrete divergence field  div = M^{-1}(Gx u + Gy v).
static double divL2(Ctx& c, const VectorXd& u, const VectorXd& v) {
    VectorXd d = c.luM.solve(c.Gx * u + c.Gy * v);
    return std::sqrt(d.dot(c.M * d));
}

static void runCase(Ctx& c, int N, double dt, double T, double nu, int pmode, bool rot,
                    bool exactNeu, VectorXd& uOut, VectorXd& vOut, bool conv = true,
                    double gradDiv = 0.0, int pressureMode = NSPRESSURE_PROJECTION) {
    Exact ex{nu};
    BCData bc = makeBC(*c.mesh, *c.edge, *c.edge2side, ex, pmode, exactNeu);
    double sigma = 8.0 * (c.fem->ord + 1) * (c.fem->ord + 1);
    NSIntegrator integ(*c.fem, *c.mesh, *c.elem2dof, *c.edge, *c.edge2side, bc, nu, dt,
                       sigma, 1.0, gradDiv, pressureMode);
    integ.rotationalPressureBC = rot;
    integ.includeConvection = conv;
    VectorXd u0, v0, p0;
    projectExact(*c.fem, *c.mesh, *c.elem2dof, ex, 0.0, u0, v0, p0);
    integ.setInitial(u0, v0, p0);
    int nsteps = static_cast<int>(std::lround(T / dt));
    for (int s = 1; s <= nsteps; ++s) integ.step(s * dt);
    uOut = integ.u(); vOut = integ.v();
}

static void study(int ord, int N, double T, double nu, int pmode, bool rot, bool exactNeu,
                  const char* label, bool conv = true, double gradDiv = 0.0,
                  int pressureMode = NSPRESSURE_PROJECTION) {
    Ctx c;
    static Mesh mesh; mesh.getMesh(1.0 / N);
    static FEM fem(ord, mesh);
    static MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    static MatrixXi edge, edge2side; mesh.getEdge2Side(edge, edge2side);
    c.fem = &fem; c.mesh = &mesh; c.elem2dof = &elem2dof; c.edge = &edge; c.edge2side = &edge2side;
    c.M = assembleScalarMassDG(fem, mesh, elem2dof);
    assembleWeakGrad(fem, mesh, elem2dof, edge, edge2side, c.Gx, c.Gy);
    c.luM.compute(c.M);

    double dtRef = T / 3200.0;
    VectorXd uRef, vRef;
    runCase(c, N, dtRef, T, nu, pmode, rot, exactNeu, uRef, vRef, conv, gradDiv, pressureMode);

    cout << "\n=== " << label << " | P" << ord << " N=" << N << " T=" << T
         << " nu=" << nu << " (dt_ref=" << scientific << setprecision(2) << dtRef << ") ===\n";
    cout << "   dt         ||u(dt)-u_ref||   rate     ||div u(dt)||\n";

    std::vector<double> dts;
    for (int k = 0; k < 6; ++k) dts.push_back(T / (25.0 * std::pow(2.0, k)));  // 1/25 .. 1/800
    double prev = -1, prevDt = 0;
    for (double dt : dts) {
        VectorXd uu, vv;
        runCase(c, N, dt, T, nu, pmode, rot, exactNeu, uu, vv, conv, gradDiv, pressureMode);
        double d = velL2Diff(fem, mesh, elem2dof, uu, vv, uRef, vRef);
        double dv = divL2(c, uu, vv);
        cout << "  " << scientific << setprecision(2) << dt << "    " << setprecision(3) << d;
        if (prev > 0) cout << "        " << fixed << setprecision(2)
                           << std::log(prev / d) / std::log(prevDt / dt);
        else          cout << "            ";
        cout << "      " << scientific << setprecision(3) << dv << "\n";
        prev = d; prevDt = dt;
    }
}

int main() {
    cout << "DG-NS TEMPORAL order diagnostic (reference-based; spatial error cancels)\n";
    double nu = 0.05;
    // Same BC layout (high-order pressure Neumann), Laplacian vs rotational viscous term.
    study(3, 16, 0.40, nu, 0, false, false, "HO-Neumann, Laplacian (OLD) ");
    study(3, 16, 0.40, nu, 0, true,  false, "HO-Neumann, rotational      ");
    // Decisive: HO-Neumann LAYOUT but feed the EXACT dp/dn (isolates reconstruction).
    study(3, 16, 0.40, nu, 0, true,  true,  "HO-Neumann, EXACT dp/dn     ");
    // Sanity: exact pressure Dirichlet everywhere (removes the pressure-BC error).
    study(3, 16, 0.40, nu, 1, true,  false, "exact-p-Dirichlet (sanity)  ");
    // New route: pressure from the consistent PPE plus a coupled grad-div velocity
    // solve.  This targets the pressure-gradient/divergence boundary layer directly.
    study(3, 16, 0.40, nu, 0, true,  false, "Direct PPE + grad-div(1.0)  ",
          true, 1.0, NSPRESSURE_DIRECT_PPE);
    study(3, 16, 0.40, nu, 0, true,  false, "Direct PPE + grad-div(10)   ",
          true, 10.0, NSPRESSURE_DIRECT_PPE);
    // H6: isolate the convection treatment.  Drop convection (unsteady Stokes) and
    // keep each pressure-BC layout; if the rate becomes a clean 2 the cap is in the
    // convection splitting; if it still degrades the cap is in the pressure/viscous
    // splitting and is INDEPENDENT of convection.
    study(3, 16, 0.40, nu, 1, true,  false, "STOKES exact-p-Dirichlet    ", false);
    study(3, 16, 0.40, nu, 0, true,  true,  "STOKES HO-Neu EXACT dp/dn   ", false);
    study(3, 16, 0.40, nu, 0, true,  false, "STOKES HO-Neumann rotational", false);
    return 0;
}
