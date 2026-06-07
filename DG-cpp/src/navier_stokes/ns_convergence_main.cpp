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

// ---------------------------------------------------------------------------
// Space-time convergence study for the DG Navier-Stokes solver, using the
// Taylor-Green decaying vortex (an exact solution of the incompressible NS
// equations) on the unit square:
//     u = -cos(x) sin(y) e^{-2 nu t},  v = sin(x) cos(y) e^{-2 nu t},
//     p = -1/4 (cos 2x + cos 2y) e^{-4 nu t}.
// Velocity is imposed (Dirichlet/Nitsche) on all four walls; pressure is pinned
// by a Dirichlet value on the right wall and the consistent high-order Neumann
// condition on the other three.  Expected rates: ||u-u_h||_{L2} ~ h^{k+1} in
// space and ~ dt^2 in time.
// ---------------------------------------------------------------------------

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
};

// L2 error of the DG velocity field against the exact solution at time t.
static double velL2Error(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
                         const VectorXd& uh, const VectorXd& vh, const Exact& ex, double t) {
    int NT = mesh.elem.rows(), locDof = fem.locDof;
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
    double err = 0.0;
    VectorXd ue(locDof), ve(locDof);
    for (int e = 0; e < NT; ++e) {
        for (int i = 0; i < locDof; ++i) { ue(i) = uh(elem2dof(e, i)); ve(i) = vh(elem2dof(e, i)); }
        Vector2d p1 = mesh.node.row(mesh.elem(e, 0));
        Vector2d p2 = mesh.node.row(mesh.elem(e, 1));
        Vector2d p3 = mesh.node.row(mesh.elem(e, 2));
        double area = fem.area(e);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d P = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            double Uh = phi[q].dot(ue), Vh = phi[q].dot(ve);
            Vector2d ue_x = ex.vel(P.x(), P.y(), t);
            double du = Uh - ue_x.x(), dv = Vh - ue_x.y();
            err += w(q) * area * (du * du + dv * dv);
        }
    }
    return std::sqrt(err);
}

// Build the boundary-condition data for the unit-square Taylor-Green test.
static BCData makeBC(const Mesh& mesh, const MatrixXi& edge, const MatrixXi& edge2side,
                     const Exact& ex) {
    int NE = edge.rows();
    BCData bc;
    bc.bcU = VectorXi::Zero(NE); bc.bcV = VectorXi::Zero(NE); bc.bcP = VectorXi::Zero(NE);
    const double tol = 1e-9;
    for (int e = 0; e < NE; ++e) {
        if (edge2side(e, 0) != -1 && edge2side(e, 1) != -1) continue;  // interior
        Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1)));
        bc.bcU(e) = BCN_DIRICHLET;            // velocity Dirichlet on all walls
        bc.bcV(e) = BCN_DIRICHLET;
        if (std::abs(m.x() - 1.0) < tol) bc.bcP(e) = BCN_DIRICHLET;    // pin pressure on right wall
        else                             bc.bcP(e) = BCN_NEUMANN_HO;   // high-order Neumann elsewhere
    }
    bc.velDir = [ex](double x, double y, double t) { return ex.vel(x, y, t); };
    bc.velAcc = [ex](double x, double y, double t) { return ex.acc(x, y, t); };
    bc.presDir = [ex](double x, double y, double t) { return ex.pres(x, y, t); };
    return bc;
}

// L2-project the exact field onto the DG space at time t (initial condition).
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

// L2 norm of the velocity difference between two DG fields on the same mesh.
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

// Run a case; optionally return the final velocity field (deterministic mesh/DOF
// numbering, so fields from different dt on the same (ord,N) are comparable).
static double runCase(int ord, int N, double dt, double T, double nu,
                      VectorXd* uOut = nullptr, VectorXd* vOut = nullptr,
                      double gradDiv = 0.0, int pressureMode = NSPRESSURE_PROJECTION,
                      double ppeDivDamping = 0.0) {
    Exact ex{nu};
    Mesh mesh; mesh.getMesh(1.0 / N);
    FEM fem(ord, mesh);
    MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    MatrixXi edge, edge2side; mesh.getEdge2Side(edge, edge2side);

    BCData bc = makeBC(mesh, edge, edge2side, ex);
    double sigma = 8.0 * (ord + 1) * (ord + 1);
    NSIntegrator integ(fem, mesh, elem2dof, edge, edge2side, bc, nu, dt, sigma, 1.0,
                       gradDiv, pressureMode, ppeDivDamping);

    VectorXd u0, v0, p0;
    projectExact(fem, mesh, elem2dof, ex, 0.0, u0, v0, p0);
    integ.setInitial(u0, v0, p0);

    int nsteps = static_cast<int>(std::lround(T / dt));
    for (int s = 1; s <= nsteps; ++s)
        if (!integ.step(s * dt)) { std::cerr << "  solve failed at step " << s << "\n"; return -1; }

    if (uOut) *uOut = integ.u();
    if (vOut) *vOut = integ.v();
    return velL2Error(fem, mesh, elem2dof, integ.u(), integ.v(), ex, nsteps * dt);
}

int main() {
    cout << "DG Navier-Stokes -- Taylor-Green space-time convergence\n";
    cout << "  u=-cos x sin y e^{-2nu t}, v=sin x cos y e^{-2nu t}, p=-1/4(cos2x+cos2y)e^{-4nu t}\n";
    double nu = 0.05;

    // ---- Spatial: tiny dt so the temporal error stays below the spatial error;
    //      N capped per order so the spatial error stays above the dt^2 floor. ----
    cout << "\n[Spatial]  dt=2.5e-4 fixed, T=0.02, expect ||u-u_h||_L2 ~ h^(k+1)\n";
    std::vector<std::vector<int>> NsByOrd = {{}, {8, 16, 32, 64}, {4, 8, 16, 32}, {4, 8, 16}};
    for (int ord = 1; ord <= 3; ++ord) {
        double prev = -1; int prevN = 0;
        for (int N : NsByOrd[ord]) {
            double err = runCase(ord, N, 2.5e-4, 0.02, nu);
            cout << "  P" << ord << "  N=" << setw(3) << N
                 << "  L2=" << scientific << setprecision(3) << err;
            if (prev > 0) cout << "   rate=" << fixed << setprecision(2)
                               << std::log(prev / err) / std::log((double)N / prevN);
            cout << "\n";
            prev = err; prevN = N;
        }
    }

    // ---- Temporal: self-convergence (Richardson) on a fixed mesh, so the spatial
    //      error cancels exactly; ||u(dt)-u(dt/2)|| should halve at rate 2^2. ----
    cout << "\n[Temporal] P3, N=16 fixed, T=0.40, DirectPPE+PPE-div-damping=10 self-convergence; expect rate 2\n";
    int to = 3, tN = 16;
    Mesh tmesh; tmesh.getMesh(1.0 / tN);
    FEM tfem(to, tmesh);
    MatrixXi tdof; int tnDof; tfem.getDOF(tmesh, tdof, tnDof);
    std::vector<double> dts = {8e-3, 4e-3, 2e-3, 1e-3};
    std::vector<VectorXd> us, vs;
    for (double dt : dts) {
        VectorXd uu, vv; runCase(to, tN, dt, 0.40, nu, &uu, &vv, 0.0, NSPRESSURE_DIRECT_PPE, 10.0);
        us.push_back(uu); vs.push_back(vv);
    }
    double prevDiff = -1;
    for (size_t i = 0; i + 1 < dts.size(); ++i) {
        double d = velL2Diff(tfem, tmesh, tdof, us[i], vs[i], us[i + 1], vs[i + 1]);
        cout << "  ||u(" << scientific << setprecision(1) << dts[i] << ")-u("
             << dts[i + 1] << ")||=" << setprecision(3) << d;
        if (prevDiff > 0) cout << "   rate=" << fixed << setprecision(2) << std::log2(prevDiff / d);
        cout << "\n";
        prevDiff = d;
    }
    return 0;
}
