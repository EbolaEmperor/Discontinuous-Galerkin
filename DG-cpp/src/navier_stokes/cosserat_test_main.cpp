#include "CosseratFilament.h"
#include "IBCoupler.h"
#include "MeshGen.h"
#include "FEM.h"
#include "Mesh.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <vector>

using namespace ns;

// Estimate the dominant period of the tip-deflection signal y(t) by averaging
// the spacing between zero-crossings of y - mean(y), using the LATTER HALF only
// (so the initial ramp-up doesn't bias the mean).
static double estimatePeriodFromSeries(const std::vector<double>& t, const std::vector<double>& y)
{
    int n = (int)t.size();
    if (n < 8) return 0.0;
    int half = n / 2;
    double mean = 0.0;
    for (int i = half; i < n; ++i) mean += y[i];
    mean /= (n - half);
    std::vector<double> crossings;
    for (int i = half + 1; i < n; ++i) {
        double a = y[i - 1] - mean, b = y[i] - mean;
        if (a < 0 && b >= 0) {
            double frac = -a / (b - a);
            crossings.push_back(t[i - 1] + frac * (t[i] - t[i - 1]));
        }
    }
    if (crossings.size() < 2) return 0.0;
    return (crossings.back() - crossings.front()) / (crossings.size() - 1);
}

// Cantilever Euler-Bernoulli beam analytic angular frequency, mode n=1,2,3:
//   omega_n = (beta_n L)^2 * sqrt(EI / (rhoLine L^4))
// where (beta_n L) are roots of  cos(bL) cosh(bL) + 1 = 0.
static double cantileverPeriod(double EI, double rhoLine, double L, int mode = 1)
{
    static const double bL[3] = { 1.8751040687119611,
                                  4.6940911329741742,
                                  7.8547574382376126 };
    double w = (bL[mode - 1] * bL[mode - 1]) * std::sqrt(EI / (rhoLine * L * L * L * L));
    return 2.0 * M_PI / w;
}

int main()
{
    // ============== Phase 1: Cantilever-rod period vs Euler-Bernoulli =====
    // Cantilever: clamped at left, free at right.  Pluck the tip with a small
    // transverse velocity and measure the resulting oscillation period.
    const double L       = 1.0;
    const double rhoLine = 1.0;       // mass per unit length
    const double EI      = 1.0;       // bending stiffness
    const double EA      = 1e6;       // very stiff axially -> nearly inextensible
    const int    N       = 60;        // segments
    const double T_end   = 30.0;      // long enough to capture many cycles
    const double dt      = 1e-3;

    CosseratFilament rod;
    rod.initStraight(N, 0.0, 0.0, L, 0.0, rhoLine, EA, EI);
    // Tiny structural damping to bleed the higher bending modes that would
    // otherwise ring forever and clutter the zero-crossing period estimate.
    // Damping ratio  zeta = (gamma) / (2 omega1) ;  with gamma=0.05 and
    // omega1 ~ 3.51, zeta ~ 0.7%, period shifts by < 0.01%.
    rod.dampStr = 0.05;
    rod.clampRoot(true);   // clamp node 0 (and node 1 -> tangent locked)

    // Initial condition: pure mode-1 shape (well-approximated by the static
    // cantilever deflection under a small tip force,  u(x) = a*(L x^2/2 - x^3/6)
    // with tip excursion a*L^3/3).  No initial velocity -> the response is a
    // clean cosine in mode 1, easy to time-resolve.
    const double tipKick = 0.005;     // small enough for linear theory
    const double a = tipKick * 3.0 / (L * L * L);
    for (int i = 0; i <= N; ++i) {
        double x = rod.X0(i, 0);
        rod.X(i, 1) = a * (L * x * x / 2.0 - x * x * x / 6.0);
    }
    rod.V.setZero();
    rod.A.setZero();

    int nsteps = (int)std::lround(T_end / dt);
    std::vector<double> tHist, yTip, energy;
    tHist.reserve(nsteps + 1);
    yTip.reserve(nsteps + 1);
    energy.reserve(nsteps + 1);
    tHist.push_back(0.0);
    yTip.push_back(rod.X(N, 1) - rod.X0(N, 1));
    energy.push_back(rod.kineticEnergy() + rod.potentialEnergy());

    for (int s = 1; s <= nsteps; ++s) {
        if (!rod.step(dt)) {
            std::cerr << "step failed at s=" << s << "\n";
            return 1;
        }
        double t = s * dt;
        tHist.push_back(t);
        yTip.push_back(rod.X(N, 1) - rod.X0(N, 1));
        energy.push_back(rod.kineticEnergy() + rod.potentialEnergy());
    }

    double Tnum = estimatePeriodFromSeries(tHist, yTip);
    double Tana = cantileverPeriod(EI, rhoLine, L, 1);

    // Energy drift over the final 80% of the run.
    int n0 = (int)(0.2 * energy.size());
    double Emin = energy[n0], Emax = energy[n0];
    for (int i = n0; i < (int)energy.size(); ++i) {
        Emin = std::min(Emin, energy[i]);
        Emax = std::max(Emax, energy[i]);
    }
    double Eref = 0.5 * (Emin + Emax);
    double drift = (Emax - Emin) / std::max(1e-30, std::abs(Eref));

    double rel = std::abs(Tnum - Tana) / Tana;
    double maxStrain = rod.maxEdgeStrain();

    std::printf("\nCantilever-rod sanity test\n");
    std::printf("  N = %d, EI = %.3g, EA = %.3g, rhoLine = %.3g, L = %.3g\n",
                N, EI, EA, rhoLine, L);
    std::printf("  dt = %.3g, T_end = %.3g (%d steps)\n", dt, T_end, nsteps);
    std::printf("  measured period  T_num = %.6f\n", Tnum);
    std::printf("  analytic period  T_ana = %.6f\n", Tana);
    std::printf("  relative error          = %.3e\n", rel);
    std::printf("  energy drift (latter 80%%) = %.3e\n", drift);
    std::printf("  max edge strain over run  = %.3e\n", maxStrain);

    // Convergence on N: also report what (β_n L)² implies for the second
    // and third mode if the discrete result happens to be tracking those.
    std::printf("\n  ref periods: T1=%.4f T2=%.4f T3=%.4f\n",
                cantileverPeriod(EI, rhoLine, L, 1),
                cantileverPeriod(EI, rhoLine, L, 2),
                cantileverPeriod(EI, rhoLine, L, 3));

    // Dump a coarse view of the tip signal so we can eyeball the waveform.
    std::printf("  tip y(t) sample:\n");
    int every = std::max(1, (int)tHist.size() / 40);
    for (int i = 0; i < (int)tHist.size(); i += every)
        std::printf("    t=%6.3f  y=%+9.5f  E=%+10.5e\n", tHist[i], yTip[i], energy[i]);

    // The energy here decays monotonically because we use a small structural
    // damping to suppress higher-mode beating (it would otherwise ruin the
    // zero-crossing period estimate); the relevant numerical-quality check is
    // therefore (a) tip-period match, and (b) bounded axial strain (no
    // unphysical stretching).  We log E(t) for inspection but don't gate on it.
    bool ok1 = (rel < 0.05) && (maxStrain < 1e-2);
    std::printf("\nPhase 1 (cantilever): %s\n", ok1 ? "PASS" : "FAIL");

    // ============== Phase 2: IB transfer-kernel adjoint identity ==========
    // Build a small DG mesh on the existing cylinder geometry (a rectangle
    // minus a disk), project a smooth analytic flow onto the DG space, place
    // a rod across the wake, give it random per-node forces, and check
    //     <F, I(u)>_rod  ==  <S(F), u>_dof    to machine precision,
    // where I = meshToRod (sample), S = rodToMesh (scatter).
    std::printf("\nPhase 2 (IB transfer-kernel adjoint identity)\n");

    Mesh mesh;
    CylinderGeom geom{-3.0, 6.0, -3.0, 3.0,  0.0, 0.0,  0.5};
    generateCylinderMesh(mesh, geom, 0.30, 6.0, 0.20);
    int ord = 2;
    FEM fem(ord, mesh, /*withHessian=*/false);
    Eigen::MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    std::printf("  mesh: %lld nodes, %lld triangles, dP%d, nDof=%d\n",
                (long long)mesh.node.rows(), (long long)mesh.elem.rows(), ord, nDof);

    MeshLocator loc; loc.build(mesh);

    // u(x,y) = a + b*x + c*y,  v(x,y) = d + e*x + f*y  -- linear, exactly
    // representable by dP_k for any k>=1.  Project to the DG space by simply
    // sampling the basis at the dof positions: the basis evaluation lambda
    // -> phi turns out to coincide with the polynomial evaluation in the
    // monomial barycentric form, so we get the same result whether we
    // initialise dofs by analytic interpolation or by L2 projection.  Easier
    // path: build dof values from triangle barycentre + interpolate.
    auto exactU = [](double x, double y){ return  0.7 + 0.3 * x - 0.5 * y; };
    auto exactV = [](double x, double y){ return -0.2 + 0.4 * x + 0.6 * y; };

    Eigen::VectorXd uDof(nDof), vDof(nDof);
    {
        // Use the monomial-coefficient form of the field: any dof-vector
        // satisfying  field(lam) = sum_j coef[j] * lam_1^a * lam_2^b * lam_3^c
        // gives the right answer when re-evaluated through the basis.  Easiest:
        // set dof values by interpolation at the *Lagrange* nodes, but the
        // codebase uses a hierarchical basis.  Workaround: project by solving
        // a small per-element linear system for the dof values that match the
        // analytic field at locDof barycentric points.
        int locDof = fem.locDof;
        // Generate locDof distinct barycentric sample points (the order-k
        // simplex Lagrange-nodes are the classical choice).  We just take a
        // simple one-parameter family: the (i,j,k=ord-i-j)/ord lattice.
        std::vector<Eigen::Vector3d> samples;
        for (int i = 0; i <= ord; ++i)
            for (int j = 0; j + i <= ord; ++j) {
                int k = ord - i - j;
                samples.emplace_back((double)i / ord, (double)j / ord, (double)k / ord);
            }
        if ((int)samples.size() != locDof) {
            std::cerr << "Lagrange-lattice size mismatch: " << samples.size()
                      << " vs locDof=" << locDof << "\n";
            return 1;
        }
        Eigen::MatrixXd Phi(locDof, locDof);
        for (int q = 0; q < locDof; ++q)
            Phi.row(q) = fem.computeBasisValue_all(samples[q].transpose()).row(0);
        Eigen::PartialPivLU<Eigen::MatrixXd> lu(Phi);
        for (int t = 0; t < (int)mesh.elem.rows(); ++t) {
            Eigen::Vector2d P1 = mesh.node.row(mesh.elem(t, 0));
            Eigen::Vector2d P2 = mesh.node.row(mesh.elem(t, 1));
            Eigen::Vector2d P3 = mesh.node.row(mesh.elem(t, 2));
            Eigen::VectorXd bU(locDof), bV(locDof);
            for (int q = 0; q < locDof; ++q) {
                Eigen::Vector3d L = samples[q];
                Eigen::Vector2d P = L(0) * P1 + L(1) * P2 + L(2) * P3;
                bU(q) = exactU(P.x(), P.y());
                bV(q) = exactV(P.x(), P.y());
            }
            Eigen::VectorXd ue = lu.solve(bU), ve = lu.solve(bV);
            for (int i = 0; i < locDof; ++i) {
                uDof(elem2dof(t, i)) = ue(i);
                vDof(elem2dof(t, i)) = ve(i);
            }
        }
    }

    // Place a rod just downstream of the cylinder, well inside the mesh.
    CosseratFilament rod2;
    int Nseg = 32;
    rod2.initStraight(Nseg, 0.55, 0.0, 2.5, 0.0, 1.0, 1.0, 1.0);
    // Nudge the shape so it isn't exactly horizontal -- exercises both
    // components of u and v.  Apply a small sinusoidal lateral displacement.
    for (int i = 0; i <= Nseg; ++i) {
        double s = (double)i / Nseg;
        rod2.X(i, 1) += 0.20 * std::sin(2.0 * M_PI * s);
    }
    // Random per-node forces.
    std::mt19937 rng(20250101);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    Eigen::MatrixXd F(Nseg + 1, 2);
    for (int k = 0; k <= Nseg; ++k) { F(k, 0) = U(rng); F(k, 1) = U(rng); }

    Eigen::VectorXd ds = rodArclengthWeights(rod2);
    std::vector<int> hint;

    // LHS  = sum_k F_k . u(X_k) * ds_k
    auto sample = meshToRod(fem, mesh, elem2dof, uDof, vDof, loc, rod2, hint);
    int nAlive = 0; for (int k = 0; k <= Nseg; ++k) nAlive += sample.alive[k];
    if (nAlive < Nseg + 1) {
        std::printf("  WARNING: %d / %d rod nodes outside the mesh -- "
                    "the adjoint test still holds for the survivors.\n",
                    nAlive, Nseg + 1);
    }
    double lhs = 0.0;
    for (int k = 0; k <= Nseg; ++k) {
        if (!sample.alive[k]) continue;
        lhs += (F(k, 0) * sample.uv(k, 0) + F(k, 1) * sample.uv(k, 1)) * ds(k);
    }

    // RHS  = sum_j (S F)_j . u_dof_j   (per-component dot product)
    Eigen::VectorXd loadX = Eigen::VectorXd::Zero(nDof);
    Eigen::VectorXd loadY = Eigen::VectorXd::Zero(nDof);
    rodToMesh(fem, mesh, elem2dof, loc, rod2, F, ds, hint, loadX, loadY);
    double rhs = loadX.dot(uDof) + loadY.dot(vDof);

    double err = std::abs(lhs - rhs) / std::max(1e-30, std::abs(lhs));
    std::printf("  LHS = <F, I(u)>_rod  = %+.12e\n", lhs);
    std::printf("  RHS = <S(F), u>_dof  = %+.12e\n", rhs);
    std::printf("  relative error       = %.3e\n", err);
    bool ok2 = (err < 1e-12);
    std::printf("\nPhase 2 (IB adjoint): %s\n", ok2 ? "PASS" : "FAIL");

    bool ok = ok1 && ok2;
    std::printf("\nResult: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
