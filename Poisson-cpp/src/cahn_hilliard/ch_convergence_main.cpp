// Cahn-Hilliard space-time convergence test (no video).
//
// Uses the method of manufactured solutions (MMS): the exact solution
//   c_e(x,y,t) = cos(pi x) cos(pi y) exp(-t)
// satisfies the no-flux BCs, and a forcing term g = d_t c_e - M Lap(mu_e) is added
// to the c-equation (see ExactSolutionCH).  We then measure the error of the
// discrete solution at the final time and fit convergence rates.
//
//   * SPATIAL sweep : fix a tiny time step and 2nd-order (SBDF2) stepping so the
//                     temporal error sits far below the spatial floor; refine h.
//                     Expect L2 rate -> k+1 and H1 rate -> k for P_k.
//   * TEMPORAL sweep: fix a fine high-order space; refine tau for time_order 1
//                     and 2.  The order is read from successive-difference
//                     Richardson  log2( ||c_tau - c_{tau/2}|| / ||c_{tau/2} - c_{tau/4}|| ),
//                     which cancels the (common) spatial error; the direct error
//                     vs c_e is printed as a cross-check.

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "Mesh.h"
#include "FEM.h"
#include "Quadrature.h"
#include "DGAssembly.h"
#include "CahnHilliard.h"
#include "ExactSolutionCH.h"

using namespace Eigen;
using namespace std;

// L2-projection of c_e(.,0) onto the DG space:  Mm c0 = \int c_e(.,0) phi.
static VectorXd projectIC(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                          const SparseMatrix<double>& Mm, const ExactSolutionCH& sol) {
    VectorXd load = assembleScalarLoad(fem, mesh, elem2dof,
                                       [&](double x, double y) { return sol.c(x, y, 0.0); });
    SimplicialLDLT<SparseMatrix<double>> ldlt;
    ldlt.compute(Mm);
    return ldlt.solve(load);
}

// L2 and H1 error of the DG field c against c_e(.,t), high-order quadrature.
static void errorCH(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                    const VectorXd& c, const ExactSolutionCH& sol, double t,
                    double& eL2, double& eH1) {
    MatrixXd quadL;
    VectorXd w;
    Quadrature::quadpts2_my(14, quadL, w);   // exact-ish for the smooth trig solution
    int nq = static_cast<int>(w.size());
    std::vector<MatrixXd> phi_q(nq), dphi_q(nq);
    for (int q = 0; q < nq; ++q) {
        phi_q[q]  = fem.computeBasisValue_all(quadL.row(q));
        dphi_q[q] = fem.computeBasisDlam_all(quadL.row(q));
    }
    int locDof = fem.locDof;
    double l2 = 0.0, h1semi = 0.0;
    VectorXd ce(locDof);
    for (int t_e = 0; t_e < mesh.elem.rows(); ++t_e) {
        for (int i = 0; i < locDof; ++i) ce(i) = c(elem2dof(t_e, i));
        Vector2d p1 = mesh.node.row(mesh.elem(t_e, 0));
        Vector2d p2 = mesh.node.row(mesh.elem(t_e, 1));
        Vector2d p3 = mesh.node.row(mesh.elem(t_e, 2));
        const MatrixXd& Dlam = fem.Dlam[t_e];
        double area = fem.area(t_e);
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            double uh = phi_q[q].row(0).dot(ce);
            double ue = sol.c(pt(0), pt(1), t);
            double d  = uh - ue;
            MatrixXd gphi = Dlam * dphi_q[q]; // 2 x locDof
            Vector2d guh  = gphi * ce;
            double gx, gy; sol.grad(pt(0), pt(1), t, gx, gy);
            double dgx = guh(0) - gx, dgy = guh(1) - gy;
            double wq = w(q) * area;
            l2     += wq * d * d;
            h1semi += wq * (dgx * dgx + dgy * dgy);
        }
    }
    eL2 = std::sqrt(l2);
    eH1 = std::sqrt(l2 + h1semi);
}

// Integrate the MMS problem to time T = nsteps*tau on a prebuilt space; return c^N.
// The separable source g = a(t) g1 + a(t)^3 g3 is applied via two precomputed loads.
static VectorXd integrateMMS(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
                             const SparseMatrix<double>& Mm, const SparseMatrix<double>& A,
                             const VectorXd& c0, const ExactSolutionCH& sol,
                             double tau, int nsteps, double mob, double eps2, double S, int order) {
    VectorXd L1 = assembleScalarLoad(fem, mesh, elem2dof, [&](double x, double y) { return sol.g1(x, y); });
    VectorXd L3 = assembleScalarLoad(fem, mesh, elem2dof, [&](double x, double y) { return sol.g3(x, y); });
    CHIntegrator integr(fem, mesh, elem2dof, Mm, A, tau, mob, eps2, S, order);
    integr.setInitial(c0);
    for (int s = 1; s <= nsteps; ++s) {
        double a = std::exp(-(s * tau));          // a(t^{n+1})
        VectorXd src = a * L1 + (a * a * a) * L3;  // \int g(.,t^{n+1}) phi
        if (!integr.step(src)) { std::cerr << "  solve failed at step " << s << "\n"; break; }
    }
    return integr.field();
}

// ||a - b||_{Mm} = sqrt((a-b)^T Mm (a-b))
static double normMm(const SparseMatrix<double>& Mm, const VectorXd& a, const VectorXd& b) {
    VectorXd d = a - b;
    return std::sqrt(d.dot(Mm * d));
}

int main() {
    cout << "Cahn-Hilliard space-time convergence test (MMS, no video)\n";
    cout << "  exact: c_e = cos(pi x) cos(pi y) exp(-t),  no-flux BC\n";
    cout << scientific << setprecision(3);

    const double eps  = 0.1;
    const double eps2 = eps * eps;
    const double mob  = 1.0;
    // MMS solution has |c_e| <= 1, so L = max|3c^2-1| = 2 and S >= L/2 = 1 suffices.
    // Using the minimal S=1 keeps the (S-proportional) first-order error constant small
    // so the order-1 rate reaches its asymptotic value 1 within the tested tau range.
    const double S    = 1.0;

    // =====================================================================
    // SPATIAL convergence: SBDF2, tiny tau so the O(tau^2) time error is far
    // below the spatial floor; refine h.  Expect L2 ~ h^{k+1}, H1 ~ h^k.
    // =====================================================================
    {
        const double tau = 1e-5;
        const int nsteps = 1000;            // T = 0.01 (tau tiny so temporal error << spatial floor)
        const double T   = nsteps * tau;
        cout << "\n================ SPATIAL convergence (SBDF2, tau=" << tau
             << ", T=" << T << ", eps=" << eps << ") ================\n";
        cout << "  expect L2 rate -> k+1, H1 rate -> k for P_k\n";

        struct Case { int ord; std::vector<int> Ns; };
        std::vector<Case> cases = {
            {1, {8, 16, 32, 64}},
            {2, {8, 16, 32}},
            {3, {8, 16, 32}},
        };

        for (const auto& cs : cases) {
            int ord = cs.ord;
            double sigma = 3.0 * ord * (ord + 1);
            cout << "\n  P" << ord << "  (sigma=" << defaultfloat << sigma << scientific << "):\n";
            std::vector<double> eL2s, eH1s;
            for (int N : cs.Ns) {
                Mesh mesh; mesh.getMesh(1.0 / N);
                FEM fem(ord, mesh);
                MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
                MatrixXi edge, e2s; mesh.getEdge2Side(edge, e2s);
                SparseMatrix<double> Mm = assembleMass_DG2D(fem, mesh, elem2dof);
                SparseMatrix<double> A  = assembleK_Poi2D(fem, mesh, elem2dof)
                                        + assembleIP_Poi2D(fem, mesh, elem2dof, edge, e2s, sigma, 1.0);
                ExactSolutionCH sol(mob, eps2);
                VectorXd c0 = projectIC(fem, mesh, elem2dof, Mm, sol);
                VectorXd cf = integrateMMS(fem, mesh, elem2dof, Mm, A, c0, sol, tau, nsteps, mob, eps2, S, 2);
                double eL2, eH1;
                errorCH(fem, mesh, elem2dof, cf, sol, T, eL2, eH1);
                eL2s.push_back(eL2); eH1s.push_back(eH1);
                cout << "    N=" << setw(3) << defaultfloat << N << scientific
                     << "  nDof=" << setw(8) << nDof
                     << "  L2=" << eL2 << "  H1=" << eH1;
                int i = static_cast<int>(eL2s.size()) - 1;
                if (i > 0) {
                    double rL2 = std::log2(eL2s[i - 1] / eL2s[i]);
                    double rH1 = std::log2(eH1s[i - 1] / eH1s[i]);
                    cout << "   rate: L2=" << fixed << setprecision(2) << rL2
                         << " H1=" << rH1 << scientific;
                }
                cout << "\n";
            }
        }
    }

    // =====================================================================
    // TEMPORAL convergence: fixed fine space; refine tau for order 1 and 2.
    // Primary metric = successive-difference Richardson (cancels spatial error).
    // =====================================================================
    {
        const int ord = 3, N = 32;   // fine high-order space: low spatial floor (~7e-8 in L2)
        const double sigma = 3.0 * ord * (ord + 1);
        const double T = 0.5;
        const double tau0 = 0.01;
        const int nlev = 6;   // tau = tau0 / 2^k

        cout << "\n================ TEMPORAL convergence (P" << ord << ", N=" << N
             << ", eps=" << eps << ", T=" << T << ") ================\n";
        cout << "  expect Richardson rate -> 1 (time_order=1), -> 2 (time_order=2 SBDF2)\n";

        Mesh mesh; mesh.getMesh(1.0 / N);
        FEM fem(ord, mesh);
        MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
        MatrixXi edge, e2s; mesh.getEdge2Side(edge, e2s);
        SparseMatrix<double> Mm = assembleMass_DG2D(fem, mesh, elem2dof);
        SparseMatrix<double> A  = assembleK_Poi2D(fem, mesh, elem2dof)
                                + assembleIP_Poi2D(fem, mesh, elem2dof, edge, e2s, sigma, 1.0);
        ExactSolutionCH sol(mob, eps2);
        VectorXd c0 = projectIC(fem, mesh, elem2dof, Mm, sol);
        cout << "  (nDof=" << nDof << ", spatial floor fixed)\n";

        for (int order = 1; order <= 2; ++order) {
            cout << "\n  time_order=" << order << (order == 2 ? " (SBDF2):\n" : " (stabilised BE):\n");
            std::vector<double> taus, eL2s;
            std::vector<VectorXd> finals;
            for (int k = 0; k < nlev; ++k) {
                double tau = tau0 / std::pow(2.0, k);
                int nsteps = static_cast<int>(std::lround(T / tau));
                VectorXd cf = integrateMMS(fem, mesh, elem2dof, Mm, A, c0, sol, tau, nsteps, mob, eps2, S, order);
                double eL2, eH1; errorCH(fem, mesh, elem2dof, cf, sol, T, eL2, eH1);
                taus.push_back(tau); eL2s.push_back(eL2); finals.push_back(cf);
                cout << "    tau=" << tau << " (" << setw(4) << defaultfloat << nsteps
                     << " steps)" << scientific << "  ||c-c_e||_L2=" << eL2;
                if (k > 0) {
                    double rdir = std::log2(eL2s[k - 1] / eL2s[k]);
                    cout << "  dir-rate=" << fixed << setprecision(2) << rdir << scientific;
                }
                cout << "\n";
            }
            // Successive-difference Richardson (spatial error cancels).
            std::vector<double> diffs;
            for (int k = 0; k + 1 < nlev; ++k) diffs.push_back(normMm(Mm, finals[k], finals[k + 1]));
            cout << "    Richardson (||c_tau - c_{tau/2}||_Mm):\n";
            for (int k = 0; k + 1 < static_cast<int>(diffs.size()); ++k) {
                double r = std::log2(diffs[k] / diffs[k + 1]);
                cout << "      tau=" << taus[k] << " -> " << taus[k + 1]
                     << "   diff=" << diffs[k] << "   rate=" << fixed << setprecision(2) << r << scientific << "\n";
            }
        }
    }

    cout << defaultfloat << "\nDone.\n";
    return 0;
}
