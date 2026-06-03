#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <filesystem>

#include <Eigen/Sparse>

#include "Mesh.h"
#include "FEM.h"
#include "DGAssembly.h"   // assembleK_Poi2D, assembleIP_Poi2D
#include "CahnHilliard.h"
#include "Json.h"         // tiny dependency-free JSON config reader

using namespace Eigen;
using namespace std;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // ======================= Parameters =======================
    // Parameters live in a JSON config file so they can be changed WITHOUT
    // recompiling. Resolution order:
    //   1) a path given on the command line:  ./cahn_hilliard my_config.json
    //   2) else "ch_config.json" in the current directory, if it exists
    //   3) else the built-in defaults below.
    // The values below are also the per-key fallback for anything omitted in the JSON.
    int      ord        = 1;                  // polynomial degree k (P_k)
    int      N          = 64;                 // cells per side; h = 1/N
    double   eps        = 0.02;               // interface-width parameter
    double   mob        = 1.0;                // mobility M
    double   S          = 2.0;                // stabilisation (S >= L/2; safe for |c| <= 1.29)
    double   beta       = 1.0;                // SIPG (symmetric); keep 1
    double   sigmaIn    = -1.0;               // SIPG penalty; <=0 / "auto" / omitted => 3*ord*(ord+1)
    double   tau        = 2.5e-5;             // time step
    int      n_steps    = 16000;              // total steps (T = n_steps*tau)
    int      save_every = 80;                 // frame cadence (-> ~200 frames)
    double   cbar       = 0.0;                // mean concentration
    double   amp        = 0.05;               // IC white-noise amplitude
    unsigned seed       = 20240601u;          // RNG seed (reproducible)
    int      Npix       = 256;                // frame resolution (MUST be even for yuv420p)
    double   cmin       = -1.0, cmax = 1.0;   // colormap range when normalize=false
    bool     normalize  = true;               // per-frame rescale to [-1,1] using max|c| (vivid early frames)
    int      time_order = 1;                   // 1 = stabilised backward-Euler (robust), 2 = SBDF2 (2nd order)
    // ---- adaptive time stepping (movie) ----
    bool     adaptive   = false;              // true: choose tau per step from the solution change
    double   t_end      = 0.0;                // adaptive: stop at this physical time (0 -> use n_steps*tau)
    double   tau_min    = 1e-6;               // adaptive: smallest tau (used during the fast early phase)
    double   tau_max    = 2e-4;               // adaptive: largest tau (used during slow coarsening)
    double   rel_change = 1.5e-3;             // adaptive: target relative change per step (sets the cadence)
    int      max_steps  = 400000;             // adaptive: safety cap on number of steps
    string   framesDir  = "ch_frames";        // output directory for PPM frames

    // ---- resolve config path ----
    string cfgPath;
    if (argc > 1)                       cfgPath = argv[1];
    else if (fs::exists("ch_config.json")) cfgPath = "ch_config.json";

    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) {
            cout << "ERROR: cannot open config file '" << cfgPath << "'\n";
            return 1;
        }
        cfgjson::Json cfg;
        try {
            cfg = cfgjson::parse(text);
        } catch (const std::exception& e) {
            cout << "ERROR parsing '" << cfgPath << "': " << e.what() << "\n";
            return 1;
        }
        if (cfg.type != cfgjson::Json::Object) {
            cout << "ERROR: config root of '" << cfgPath << "' must be a JSON object\n";
            return 1;
        }
        ord        = cfg.getInt("ord", ord);
        N          = cfg.getInt("N", N);
        eps        = cfg.getNumber("eps", eps);
        mob        = cfg.getNumber("mob", mob);
        S          = cfg.getNumber("S", S);
        beta       = cfg.getNumber("beta", beta);
        if (cfg.hasNumber("sigma")) sigmaIn = cfg.getNumber("sigma", sigmaIn); // else "auto"/omitted -> formula
        tau        = cfg.getNumber("tau", tau);
        n_steps    = cfg.getInt("n_steps", n_steps);
        save_every = cfg.getInt("save_every", save_every);
        cbar       = cfg.getNumber("cbar", cbar);
        amp        = cfg.getNumber("amp", amp);
        seed       = static_cast<unsigned>(cfg.getLong("seed", static_cast<long>(seed)));
        Npix       = cfg.getInt("Npix", Npix);
        cmin       = cfg.getNumber("cmin", cmin);
        cmax       = cfg.getNumber("cmax", cmax);
        normalize  = cfg.getBool("normalize", normalize);
        time_order = cfg.getInt("time_order", time_order);
        adaptive   = cfg.getBool("adaptive", adaptive);
        t_end      = cfg.getNumber("t_end", t_end);
        tau_min    = cfg.getNumber("tau_min", tau_min);
        tau_max    = cfg.getNumber("tau_max", tau_max);
        rel_change = cfg.getNumber("rel_change", rel_change);
        max_steps  = cfg.getInt("max_steps", max_steps);
        framesDir  = cfg.getString("frames_dir", framesDir);
    }
    if (t_end <= 0.0) t_end = n_steps * tau;   // default horizon = fixed-mode T

    double sigma = (sigmaIn > 0.0) ? sigmaIn : 3.0 * ord * (ord + 1);

    // ---- validate ----
    bool ok = true;
    if (ord < 1)         { cout << "config error: ord must be >= 1\n";        ok = false; }
    if (N < 1)           { cout << "config error: N must be >= 1\n";          ok = false; }
    if (eps <= 0.0)      { cout << "config error: eps must be > 0\n";         ok = false; }
    if (mob <= 0.0)      { cout << "config error: mob must be > 0\n";         ok = false; }
    if (tau <= 0.0)      { cout << "config error: tau must be > 0\n";         ok = false; }
    if (n_steps < 0)     { cout << "config error: n_steps must be >= 0\n";    ok = false; }
    if (save_every <= 0) { cout << "config error: save_every must be >= 1\n"; ok = false; }
    if (Npix < 2)        { cout << "config error: Npix must be >= 2\n";       ok = false; }
    if (time_order != 1 && time_order != 2) { cout << "config error: time_order must be 1 or 2\n"; ok = false; }
    if (adaptive) {
        if (tau_min <= 0.0 || tau_max < tau_min) { cout << "config error: need 0 < tau_min <= tau_max\n"; ok = false; }
        if (rel_change <= 0.0) { cout << "config error: rel_change must be > 0\n"; ok = false; }
        if (t_end <= 0.0)      { cout << "config error: t_end must be > 0 (adaptive)\n"; ok = false; }
        if (max_steps < 1)     { cout << "config error: max_steps must be >= 1\n"; ok = false; }
    }
    if (!ok) return 1;
    if (Npix % 2 != 0) {
        cout << "note: Npix must be even for yuv420p; rounding " << Npix << " -> " << Npix + 1 << "\n";
        ++Npix;
    }

    const double eps2 = eps * eps;
    const double h    = 1.0 / N;

    cout << "Cahn-Hilliard mixed-DG / stabilised IMEX (C++)\n";
    cout << "  PDE: dc/dt = M*Lap(mu),  mu = c^3 - c - eps^2 Lap(c),  no-flux BC\n";
    cout << "  config: " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  ord=" << ord << "  N=" << N << "  h=" << h << "  eps=" << eps
         << "  M=" << mob << "  S=" << S << "  sigma=" << sigma
         << "  tau=" << tau << "\n";
    if (adaptive)
        cout << "  mode=ADAPTIVE  t_end=" << t_end << "  save_every=" << save_every << " steps";
    else
        cout << "  mode=fixed  n_steps=" << n_steps << " (T=" << n_steps * tau << ")  save_every=" << save_every;
    cout << "  Npix=" << Npix << "  amp=" << amp << "  normalize=" << (normalize ? "true" : "false") << "\n";

    // ======================= Mesh + DG space =======================
    Mesh mesh;
    mesh.getMesh(h);                    // [0,1]^2 structured triangulation
    FEM fem(ord, mesh);
    MatrixXi elem2dof;
    int nDof;
    fem.getDOF(mesh, elem2dof, nDof);   // discontinuous (per-element) numbering
    MatrixXi edge, edge2side;
    mesh.getEdge2Side(edge, edge2side);
    cout << "  NT=" << mesh.elem.rows() << "  nDof=" << nDof
         << "  (block system 2*nDof=" << 2 * nDof << ")\n";

    // ======================= Operators (assembled once) =======================
    auto t_asm0 = chrono::high_resolution_clock::now();
    SparseMatrix<double> Mm = assembleMass_DG2D(fem, mesh, elem2dof);
    SparseMatrix<double> K  = assembleK_Poi2D(fem, mesh, elem2dof);
    SparseMatrix<double> P  = assembleIP_Poi2D(fem, mesh, elem2dof, edge, edge2side, sigma, beta);
    SparseMatrix<double> A  = K + P;    // symmetric SIPG operator for -Laplacian, A*1 = 0

    // Mass-conservation prerequisite: the constant must be in the nullspace of A.
    {
        VectorXd ones = VectorXd::Ones(nDof);
        double r = (A * ones).cwiseAbs().maxCoeff();
        cout << "  ||A*1||_inf = " << scientific << r << defaultfloat
             << " (must be ~0 for exact mass conservation)\n";
        if (r > 1e-8)
            cout << "  WARNING: A*1 != 0 -> discrete mass will drift!\n";
    }

    // ======================= Initial condition + output helpers =======================
    VectorXd c0 = initSpinodal(nDof, cbar, amp, seed);
    const double mass0 = computeMassCH(Mm, c0);

    fs::create_directories(framesDir);
    for (const auto& e : fs::directory_iterator(framesDir))   // clear stale frames
        if (e.path().extension() == ".ppm") fs::remove(e.path());

    auto writeFrame = [&](const VectorXd& cf, int frameIdx) {
        char fn[512];
        snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", framesDir.c_str(), frameIdx);
        double lo = cmin, hi = cmax;
        if (normalize) {
            // Symmetric per-frame rescale: map [-max|c|, +max|c|] onto the colormap so
            // c=0 stays the neutral mid colour. Makes the faint early field show full
            // contrast (the floor avoids a degenerate all-zero frame).
            double m = std::max(1e-3, cf.cwiseAbs().maxCoeff());
            lo = -m; hi = m;
        }
        writeFramePPM(fn, mesh, elem2dof, cf, Npix, lo, hi);
    };
    auto report = [&](const VectorXd& cf, int step, double t, double tau_now, int frameIdx) {
        double E    = computeEnergyCH(fem, mesh, elem2dof, cf, A, eps2);
        double mass = computeMassCH(Mm, cf);
        double maxc = cf.cwiseAbs().maxCoeff();
        cout << "  step " << setw(6) << step
             << "  t=" << fixed << setprecision(5) << t
             << "  tau=" << scientific << setprecision(2) << tau_now
             << "  frame=" << setw(4) << frameIdx
             << "  E=" << fixed << setprecision(6) << E
             << "  drift=" << scientific << setprecision(2) << (mass - mass0)
             << "  max|c|=" << fixed << setprecision(4) << maxc << defaultfloat;
        if (3.0 * maxc * maxc - 1.0 > 2.0 * S) cout << "  [WARN E-stab]";
        cout << "\n";
    };

    // ======================= Time loop =======================
    int frame = 0;
    auto t_run0 = chrono::high_resolution_clock::now();

    if (!adaptive) {
        // ---- fixed-step, order 1 (stabilised BE) or 2 (SBDF2); J factorised once ----
        CHIntegrator integr(fem, mesh, elem2dof, Mm, A, tau, mob, eps2, S, time_order);
        integr.setInitial(c0);
        cout << "  time_order=" << time_order << (time_order == 2 ? " (SBDF2)" : " (stabilised BE)")
             << ";  assemble + factorise: "
             << chrono::duration<double>(chrono::high_resolution_clock::now() - t_asm0).count() << "s\n";
        cout << "\nEvolving (fixed tau; mass drift ~ machine eps; energy non-increasing for order 1):\n";
        writeFrame(integr.field(), frame); report(integr.field(), 0, 0.0, tau, frame); ++frame;
        for (int step = 1; step <= n_steps; ++step) {
            if (!integr.step()) { cout << "ERROR: solve failed at step " << step << "\n"; return 1; }
            if (step % save_every == 0 || step == n_steps) {
                writeFrame(integr.field(), frame); report(integr.field(), step, step * tau, tau, frame); ++frame;
            }
        }
    } else {
        // ---- adaptive backward-Euler: small tau when fast, large when slow ----
        if (time_order == 2)
            cout << "  note: adaptive mode uses the 1st-order (backward-Euler) scheme; time_order=2 ignored.\n";
        CHAdaptiveStepper stepper(fem, mesh, elem2dof, Mm, A, mob, eps2, S, tau_min, tau_max);
        stepper.setInitial(c0);
        cout << "  adaptive BE: tau in [" << scientific << setprecision(1) << tau_min << ", " << tau_max
             << "]  rel_change=" << rel_change << ";  operators assembled in "
             << fixed << setprecision(3) << chrono::duration<double>(chrono::high_resolution_clock::now() - t_asm0).count()
             << "s (factorisation is lazy per tau-level)\n" << defaultfloat;
        cout << "\nEvolving (adaptive tau grows as the dynamics slow; mass drift ~ machine eps):\n";
        writeFrame(stepper.field(), frame); report(stepper.field(), 0, 0.0, tau_min, frame); ++frame;
        double t = 0.0, tauDesired = tau_min;
        int step = 0;
        while (t < t_end && step < max_steps) {
            if (!stepper.step(tauDesired)) { cout << "ERROR: solve failed at step " << step << "\n"; return 1; }
            ++step;
            double aTau = stepper.lastTau();
            double r    = stepper.lastRelChange();
            t += aTau;
            // controller: keep the per-step relative change near rel_change (clamp the per-step ratio)
            double factor = std::min(2.0, std::max(0.5, 0.9 * rel_change / std::max(r, 1e-12)));
            tauDesired = std::min(tau_max, std::max(tau_min, aTau * factor));
            if (step % save_every == 0 || t >= t_end || step >= max_steps) {
                writeFrame(stepper.field(), frame); report(stepper.field(), step, t, aTau, frame); ++frame;
            }
        }
        cout << "  (" << step << " steps; " << stepper.numFactorizations()
             << " factorisations across tau levels)\n";
    }
    auto t_run1 = chrono::high_resolution_clock::now();

    cout << "\nDone. " << frame << " frames in '" << framesDir << "/'. "
         << "Time-stepping wall: " << chrono::duration<double>(t_run1 - t_run0).count() << "s\n";
    cout << "\nAssemble the movie with ffmpeg:\n"
         << "  ffmpeg -y -framerate 25 -i " << framesDir << "/frame_%05d.ppm"
         << " -c:v libx264 -pix_fmt yuv420p -crf 18 cahn_hilliard.mp4\n";
    return 0;
}
