// ===========================================================================
// Space-time convergence study for the DG + 2nd-order IMEX(ARS(2,2,2)) Euler
// solver, using the isentropic (Shu) vortex -- an exact, infinitely smooth
// solution of the 2-D compressible Euler equations.  A stationary vortex on a
// uniform mean flow (rho_inf=1, u_inf=v_inf=1, p_inf=1, beta=5) is convected
// without change of shape, so the exact solution at time t is the t=0 vortex
// rigidly translated by (u_inf,v_inf)*t.  The exact solution is imposed as the
// Rusanov exterior (ghost) state on every boundary edge -- a consistent,
// high-order boundary condition equivalent to periodic BC.
//
//   * Spatial order:  fixed (tiny) dt, refine the mesh -> ||rho-rho_h||_L2 ~ h^{k+1}.
//   * Temporal order: fixed mesh, halve dt -> successive-difference (Richardson)
//                     rate -> 2 (the spatial error cancels in the differences).
//
// Artificial viscosity and the positivity limiter are OFF here (the flow is
// smooth; limiters would cap the asymptotic order).
// ===========================================================================
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "EulerDG.h"
#include "FEM.h"
#include "Mesh.h"

using namespace Eigen;
using namespace std;
using namespace euler;

struct Vortex {
    double uinf = 1.0, vinf = 1.0, beta = 5.0, x0 = 5.0, y0 = 5.0;
    Vector4d prim(double x, double y, double t) const {
        const double g = GAMMA, PI = M_PI;
        double xb = x - x0 - uinf * t, yb = y - y0 - vinf * t;
        double r2 = xb * xb + yb * yb;
        double f  = std::exp(0.5 * (1.0 - r2));
        double u  = uinf - (beta / (2 * PI)) * yb * f;
        double v  = vinf + (beta / (2 * PI)) * xb * f;
        double T  = 1.0 - (g - 1.0) * beta * beta / (8.0 * g * PI * PI) * std::exp(1.0 - r2);
        double rho = std::pow(T, 1.0 / (g - 1.0));
        double p   = std::pow(rho, g);
        return Vector4d(rho, u, v, p);
    }
    Vector4d cons(double x, double y, double t) const {
        Vector4d q = prim(x, y, t);
        return primToCons(q(0), q(1), q(2), q(3));
    }
};

// run the vortex to time T on an nx x ny mesh of [0,L]^2, return the DG state.
static MatrixXd runVortex(int ord, int nx, double L, double T, double dt,
                          const Vortex& vx, FEM*& femOut, Mesh*& meshOut, MatrixXi*& e2dOut) {
    Mesh* mesh = new Mesh();
    makeRectMesh(*mesh, 0, L, 0, L, nx, nx);
    FEM* fem = new FEM(ord, *mesh);
    MatrixXi* e2d = new MatrixXi; int nDof; fem->getDOF(*mesh, *e2d, nDof);
    MatrixXi edge, e2s; mesh->getEdge2Side(edge, e2s);
    VectorXi tag = VectorXi::Zero(edge.rows());
    EulerConfig cfg; cfg.use_av = false; cfg.use_positivity = false;
    EulerDG* dg = new EulerDG(*fem, *mesh, *e2d, edge, e2s, tag, cfg);
    dg->setState(projectInitial(*fem, *mesh, *e2d, [&](double x, double y){ return vx.prim(x, y, 0.0); }));
    ExteriorStateFn bc = [&](double x, double y, double t, const Vector4d&, double, double, int){ return vx.cons(x, y, t); };
    int nsteps = std::max(1, (int)std::lround(T / dt));
    double h = dt; // keep dt fixed (T may not be an exact multiple; use nsteps*dt ~ T)
    for (int s = 1; s <= nsteps; ++s) dg->step(h, s * h, bc);
    MatrixXd U = dg->state();
    femOut = fem; meshOut = mesh; e2dOut = e2d;
    delete dg;
    return U;
}

int main() {
    cout << "=========================================================================\n";
    cout << " DG + 2nd-order IMEX (ARS(2,2,2)) Euler -- isentropic-vortex convergence\n";
    cout << "=========================================================================\n";
    const double L = 10.0;
    Vortex vx;

    // ----------------------- spatial order (fixed tiny dt) -----------------------
    cout << "\n[1] SPATIAL ORDER   (T=1, fixed dt=5e-4 so temporal error is negligible)\n";
    const double Tspace = 1.0, dtSpace = 5e-4;
    for (int ord : {1, 2, 3}) {
        cout << "  dP" << ord << ":\n";
        double prevErr = 0, prevH = 0;
        for (int nx : {16, 24, 32, 48}) {
            FEM* fem; Mesh* mesh; MatrixXi* e2d;
            MatrixXd U = runVortex(ord, nx, L, Tspace, dtSpace, vx, fem, mesh, e2d);
            double err = l2Error(*fem, *mesh, *e2d, U, 0, [&](double x, double y){ return vx.cons(x, y, Tspace); });
            double h = L / nx;
            double rate = (prevErr > 0) ? std::log(prevErr / err) / std::log(prevH / h) : 0.0;
            cout << "     nx=" << setw(3) << nx << "  h=" << fixed << setprecision(4) << h
                 << "  ||rho-rho_h||_L2=" << scientific << setprecision(3) << err
                 << "  rate=" << fixed << setprecision(2) << (prevErr > 0 ? rate : 0.0) << "\n";
            prevErr = err; prevH = h;
            delete fem; delete mesh; delete e2d;
        }
        cout << "     (expected rate -> " << ord + 1 << ")\n";
    }

    // ----------------------- temporal order (Richardson) -----------------------
    cout << "\n[2] TEMPORAL ORDER  (fixed dP3 mesh nx=24; successive-difference / Richardson)\n";
    const int ordT = 3, nxT = 24;
    const double Ttime = 1.0;
    std::vector<double> dts = {1.0/100, 1.0/200, 1.0/400, 1.0/800, 1.0/1600};
    std::vector<MatrixXd> sols;
    for (double dt : dts) {
        FEM* fem; Mesh* mesh; MatrixXi* e2d;
        MatrixXd U = runVortex(ordT, nxT, L, Ttime, dt, vx, fem, mesh, e2d);
        sols.push_back(U);
        delete fem; delete mesh; delete e2d;
    }
    // build a throwaway space to measure differences (same mesh for all)
    {
        Mesh mesh; makeRectMesh(mesh, 0, L, 0, L, nxT, nxT);
        FEM fem(ordT, mesh); MatrixXi e2d; int nDof; fem.getDOF(mesh, e2d, nDof);
        auto diffL2 = [&](const MatrixXd& A, const MatrixXd& B){
            return l2Error(fem, mesh, e2d, A, 0, [&](double x, double y){
                // exact arg unused; we hijack l2Error by passing B as 'exact' through a closure
                (void)x; (void)y; return Vector4d::Zero(); }); };
        (void)diffL2;
        // proper successive differences via direct L2 of (sols[i]-sols[i+1]) density
        auto densDiff = [&](const MatrixXd& A, const MatrixXd& B){
            int NT = mesh.elem.rows(), locDof = fem.locDof;
            MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
            int nq = w.size(); std::vector<RowVectorXd> phi(nq);
            for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
            double e = 0;
            for (int t = 0; t < NT; ++t) { double area = fem.area(t);
                for (int q = 0; q < nq; ++q) { double da = 0;
                    for (int i = 0; i < locDof; ++i) da += phi[q](i) * (A(e2d(t,i),0) - B(e2d(t,i),0));
                    e += w(q) * area * da * da; } }
            return std::sqrt(e); };
        std::vector<double> d;
        for (size_t i = 0; i + 1 < sols.size(); ++i) d.push_back(densDiff(sols[i], sols[i+1]));
        double prev = 0;
        for (size_t i = 0; i < d.size(); ++i) {
            double rate = (i > 0) ? std::log(prev / d[i]) / std::log(2.0) : 0.0;
            cout << "     dt=" << scientific << setprecision(2) << dts[i] << " -> " << dts[i+1]
                 << "   ||rho(dt)-rho(dt/2)||_L2=" << setprecision(3) << d[i]
                 << "   rate=" << fixed << setprecision(2) << (i > 0 ? rate : 0.0) << "\n";
            prev = d[i];
        }
        cout << "     (expected rate -> 2)\n";
    }
    cout << "\nDone.\n";
    return 0;
}
