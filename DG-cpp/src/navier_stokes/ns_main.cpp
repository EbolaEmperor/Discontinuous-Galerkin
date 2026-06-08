#include <chrono>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "FEM.h"
#include "Mesh.h"
#include "MeshGen.h"
#include "NavierStokes.h"
#include "Json.h"            // tiny dependency-free JSON config reader (shared with CH)

using namespace Eigen;
using namespace std;
using namespace ns;
namespace fs = std::filesystem;

// ===========================================================================
// Classic 2-D flow past a circular cylinder (the von Karman vortex street),
// solved with the DG / high-order-splitting Navier-Stokes integrator on an
// automatically generated, graded, body-fitted triangular mesh.
//
//   * uniform free-stream inflow on the left,   u = (Uinf, 0)
//   * "do-nothing" outflow on the right         (Neumann u, p = 0)
//   * slip (symmetry) side walls                (v = 0, free tangential u)
//   * no-slip cylinder                          (u = 0)
//
//   Re = Uinf * D / nu,   D = 2*radius.  Default Re = 100 -> periodic shedding.
//
// Output: a vorticity movie (PPM frames -> ffmpeg) plus a drag/lift/Strouhal
// time-history.  All parameters live in a JSON config (no recompile).
// ===========================================================================

int main(int argc, char** argv) {
    // ----------------------- parameters (JSON-overridable) -----------------------
    int    ord    = 2;                 // polynomial degree k (dP_k)
    double radius = 0.5;               // cylinder radius (diameter D = 1)
    double cx = 0.0, cy = 0.0;         // cylinder centre
    double xa = -5.0, xb = 18.0;       // domain in x  (5D upstream, 18D downstream)
    double ya = -6.0, yb = 6.0;        // domain in y  (slip walls)
    double h      = 0.07;              // target element size on the cylinder
    double farRatio = 9.0;             // far-field size / h
    double grade  = 0.22;              // mesh grading rate
    double Re     = 100.0;             // Reynolds number
    double Uinf   = 1.0;               // free-stream speed
    double cfl    = 0.5;               // sets dt = cfl*h/((2k+1)Uinf) when dt<=0
    double dt     = 0.0;               // time step (<=0 -> from CFL)
    double t_end  = 90.0;              // final time
    int    save_every = 0;             // frame cadence in steps (0 -> ~auto)
    int    n_frames = 240;             // target number of frames (used if save_every<=0)
    double sigmaFac = 8.0;             // SIPG/Nitsche penalty = sigmaFac*(k+1)^2
    double gradDiv = 0.0;              // optional coupled grad-div stabilisation strength
    double ppeDivDamping = 10.0;       // cheap PPE divergence damping (no coupled velocity solve)
    string pressureModeName = "direct_ppe"; // "direct_ppe" | "projection"
    double rampTime = 0.5;             // smooth inflow ramp; >0 gives compatible zero initial data
    double perturb = 0.30;             // initial v-perturbation amplitude (triggers shedding)
    // rendering
    int    Wpix = 1100;                // frame width in pixels
    double rxa = -3.0, rxb = 16.0;     // render window
    double rya = -4.0, ryb = 4.0;
    double vortClip = 4.0;             // vorticity colour range +/- vortClip
    // Output channels.  Each entry is one of  "vorticity" | "speed" | "flow"
    // and produces an independent frame directory + ffmpeg command line.  If
    // the user supplies neither `outputs` nor the legacy `field`/`render_flow`
    // keys, default to vorticity + flow (backwards-compatible "lots of pretty").
    std::vector<string> outputs;       // filled from JSON (or legacy keys) below
    string framesDirVort  = "ns_frames";        // legacy default for vorticity / speed
    string framesDirFlow  = "ns_flow_frames";
    int    nParticles     = 1500;
    unsigned int particleSeed = 12345u;
    int    trailLen       = 24;       // ribbon length per particle
    int    trailStride    = 1;        // sample every k-th step (k>1 stretches the streak)
    double bgDim          = 0.12;     // 0=black bg, 1=full viridis
    // Legacy (still honoured if `outputs` is missing): pick exactly one colour
    // field via `field`, optionally also enable the streakline movie via
    // `render_flow`.
    string legacyField    = "";        // "" -> not set in config
    int    legacyRenderFlow = -1;      // -1 -> not set in config

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("ns_config.json")) cfgPath = "ns_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        ord = cfg.getInt("ord", ord);
        radius = cfg.getNumber("radius", radius);
        cx = cfg.getNumber("cx", cx); cy = cfg.getNumber("cy", cy);
        xa = cfg.getNumber("xa", xa); xb = cfg.getNumber("xb", xb);
        ya = cfg.getNumber("ya", ya); yb = cfg.getNumber("yb", yb);
        h = cfg.getNumber("h", h);
        farRatio = cfg.getNumber("far_ratio", farRatio);
        grade = cfg.getNumber("grade", grade);
        Re = cfg.getNumber("Re", Re);
        Uinf = cfg.getNumber("Uinf", Uinf);
        cfl = cfg.getNumber("cfl", cfl);
        dt = cfg.getNumber("dt", dt);
        t_end = cfg.getNumber("t_end", t_end);
        save_every = cfg.getInt("save_every", save_every);
        n_frames = cfg.getInt("n_frames", n_frames);
        sigmaFac = cfg.getNumber("sigma_fac", sigmaFac);
        gradDiv = cfg.getNumber("grad_div", gradDiv);
        ppeDivDamping = cfg.getNumber("ppe_div_damping", ppeDivDamping);
        pressureModeName = cfg.getString("pressure_mode", pressureModeName);
        rampTime = cfg.getNumber("ramp_time", rampTime);
        perturb = cfg.getNumber("perturb", perturb);
        Wpix = cfg.getInt("Wpix", Wpix);
        rxa = cfg.getNumber("render_xa", rxa); rxb = cfg.getNumber("render_xb", rxb);
        rya = cfg.getNumber("render_ya", rya); ryb = cfg.getNumber("render_yb", ryb);
        vortClip = cfg.getNumber("vort_clip", vortClip);
        // ---- output selection ----
        // Preferred form: an array, e.g.  "outputs": ["vorticity", "flow"]
        // The order is also the order in which frame dirs / ffmpeg lines are emitted.
        if (cfg.contains("outputs") && cfg.obj.at("outputs").type == cfgjson::Json::Array) {
            for (const auto& it : cfg.obj.at("outputs").arr)
                if (it.type == cfgjson::Json::String) outputs.push_back(it.str);
        }
        // Legacy form, only consulted when `outputs` is missing.
        legacyField      = cfg.getString("field", legacyField);
        legacyRenderFlow = cfg.contains("render_flow")
                           ? (cfg.getInt("render_flow", 0) != 0 ? 1 : 0)
                           : -1;
        framesDirVort = cfg.getString("frames_dir", framesDirVort);
        framesDirFlow = cfg.getString("flow_frames_dir", framesDirFlow);
        nParticles    = cfg.getInt("n_particles", nParticles);
        particleSeed  = (unsigned int)cfg.getInt("particle_seed", (int)particleSeed);
        trailLen      = cfg.getInt("trail_len", trailLen);
        trailStride   = cfg.getInt("trail_stride", trailStride);
        bgDim         = cfg.getNumber("bg_dim", bgDim);
    }
    // Resolve `outputs` from legacy keys if needed.
    if (outputs.empty()) {
        if (!legacyField.empty()) outputs.push_back(legacyField);
        else                      outputs.push_back("vorticity");
        // legacy default kept the streakline movie on; only suppress if user
        // explicitly set render_flow=0
        if (legacyRenderFlow != 0) outputs.push_back("flow");
    }
    // Normalise + validate.
    for (auto& s : outputs) {
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s != "vorticity" && s != "speed" && s != "flow") {
            cout << "ERROR: outputs entry '" << s
                 << "' must be 'vorticity', 'speed', or 'flow'\n";
            return 1;
        }
    }

    const double D = 2.0 * radius;
    const double nu = Uinf * D / Re;
    if (dt <= 0.0) dt = cfl * h / ((2.0 * ord + 1.0) * Uinf);
    int nsteps = std::max(1, (int)std::lround(t_end / dt));
    if (save_every <= 0) save_every = std::max(1, nsteps / std::max(1, n_frames));
    double sigma = sigmaFac * (ord + 1) * (ord + 1);
    int pressureMode = (pressureModeName == "projection") ? NSPRESSURE_PROJECTION : NSPRESSURE_DIRECT_PPE;

    cout << "DG incompressible Navier-Stokes -- flow past a cylinder (von Karman street)\n";
    cout << "  config:  " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  domain:  [" << xa << "," << xb << "] x [" << ya << "," << yb << "]"
         << "  cylinder r=" << radius << " at (" << cx << "," << cy << ")\n";
    cout << "  Re=" << Re << "  Uinf=" << Uinf << "  D=" << D << "  nu=" << nu
         << "  dP" << ord << "  h=" << h << "\n";
    cout << "  dt=" << dt << "  t_end=" << t_end << "  steps=" << nsteps
         << "  save_every=" << save_every << "  sigma=" << sigma << "\n";
    cout << "  pressure_mode=" << (pressureMode == NSPRESSURE_DIRECT_PPE ? "direct_ppe" : "projection")
         << "  grad_div=" << gradDiv << "  ppe_div_damping=" << ppeDivDamping
         << "  ramp_time=" << rampTime << "\n";

    // ----------------------- mesh + DG space -----------------------
    CylinderGeom geom{xa, xb, ya, yb, cx, cy, radius};
    Mesh mesh;
    generateCylinderMesh(mesh, geom, h, farRatio, grade);
    FEM fem(ord, mesh);
    MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    MatrixXi edge, edge2side; mesh.getEdge2Side(edge, edge2side);
    VectorXi tag = classifyEdges(mesh, edge, edge2side, geom);
    cout << "  nDof=" << nDof << " (x3 fields = " << 3 * nDof << ")\n";

    // ----------------------- boundary conditions from the edge tags -----------------------
    int NE = edge.rows();
    BCData bc;
    bc.bcU = VectorXi::Zero(NE); bc.bcV = VectorXi::Zero(NE); bc.bcP = VectorXi::Zero(NE);
    VectorXi isCyl = VectorXi::Zero(NE);
    for (int e = 0; e < NE; ++e) {
        switch (tag(e)) {
            case BD_INFLOW: bc.bcU(e)=BCN_DIRICHLET; bc.bcV(e)=BCN_DIRICHLET; bc.bcP(e)=BCN_NEUMANN_HO;   break;
            case BD_CYL:    bc.bcU(e)=BCN_DIRICHLET; bc.bcV(e)=BCN_DIRICHLET; bc.bcP(e)=BCN_NEUMANN_HO; isCyl(e)=1; break;
            case BD_WALL:   bc.bcU(e)=BCN_NEUMANN;   bc.bcV(e)=BCN_DIRICHLET; bc.bcP(e)=BCN_NEUMANN_ZERO; break;
            case BD_OUTFLOW:bc.bcU(e)=BCN_NEUMANN;   bc.bcV(e)=BCN_NEUMANN;   bc.bcP(e)=BCN_NEUMANN_HO;   break;
            default: break;  // interior
        }
    }
    auto ramp = [](double t, double tr) {
        if (tr <= 0.0 || t >= tr) return 1.0;
        double s = std::max(0.0, t / tr);
        return s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
    };
    auto rampDer = [](double t, double tr) {
        if (tr <= 0.0 || t <= 0.0 || t >= tr) return 0.0;
        double s = t / tr;
        return (30.0 * s * s - 60.0 * s * s * s + 30.0 * s * s * s * s) / tr;
    };
    double rcyl = radius, ccx = cx, ccy = cy, htol = 0.3 * h, Uin = Uinf, tr = rampTime;
    bc.velDir = [=](double x, double y, double t) -> Vector2d {
        if (std::hypot(x - ccx, y - ccy) < rcyl + htol) return Vector2d(0.0, 0.0); // no-slip cylinder
        return Vector2d(Uin * ramp(t, tr), 0.0);                                   // inflow (wall: v=0 used)
    };
    bc.velAcc = [=](double, double, double t) {
        return Vector2d(Uin * rampDer(t, tr), 0.0);
    };
    bc.presDir = [](double, double, double) { return 0.0; };                        // used only if pressure Dirichlet is selected

    NSIntegrator integ(fem, mesh, elem2dof, edge, edge2side, bc, nu, dt, sigma, 1.0,
                       gradDiv, pressureMode, ppeDivDamping);

    // ----------------------- initial condition (uniform flow + wake kick) -----------------------
    SparseMatrix<double> Msc = assembleScalarMassDG(fem, mesh, elem2dof);
    SimplicialLDLT<SparseMatrix<double>> luMsc(Msc);
    auto project = [&](const std::function<double(double,double)>& f) -> VectorXd {
        int NT = mesh.elem.rows(), locDof = fem.locDof;
        MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
        int nq = (int)w.size();
        std::vector<RowVectorXd> phi(nq);
        for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
        VectorXd b = VectorXd::Zero(nDof);
        for (int e = 0; e < NT; ++e) {
            Vector2d p1 = mesh.node.row(mesh.elem(e,0)), p2 = mesh.node.row(mesh.elem(e,1)), p3 = mesh.node.row(mesh.elem(e,2));
            double area = fem.area(e); VectorXd le = VectorXd::Zero(locDof);
            for (int q = 0; q < nq; ++q) {
                Vector3d lam = quadL.row(q).transpose();
                Vector2d P = lam(0)*p1 + lam(1)*p2 + lam(2)*p3;
                le.noalias() += (w(q)*area*f(P.x(),P.y())) * phi[q].transpose();
            }
            for (int i = 0; i < locDof; ++i) b(elem2dof(e,i)) += le(i);
        }
        return luMsc.solve(b);
    };
    VectorXd u0, v0;
    if (rampTime > 0.0) {
        const double xc = cx + 3.0 * D;
        const double sig2 = D * D;
        const double amp = 0.10 * perturb * Uinf * D;
        u0 = project([&](double x, double y){
            double X = x - xc, Y = y - cy;
            double psi = amp * std::exp(-(X * X + 4.0 * Y * Y) / sig2);
            return psi * (-8.0 * Y / sig2);                 // d psi / dy
        });
        v0 = project([&](double x, double y){
            double X = x - xc, Y = y - cy;
            double psi = amp * std::exp(-(X * X + 4.0 * Y * Y) / sig2);
            return psi * (2.0 * X / sig2);                  // -d psi / dx
        });
    } else {
        // Legacy instantaneous-start initial condition.
        u0 = project([&](double, double){ return Uinf; });
        v0 = project([&](double x, double y){
            double xi = x - cx;
            if (xi <= 0) return 0.0;
            return perturb * Uinf * std::sin(M_PI * xi / (4.0 * D))
                   * std::exp(-(xi*xi + 4.0*(y-cy)*(y-cy)) / (9.0 * D * D));
        });
    }
    integ.setInitial(u0, v0);

    // ----------------------- output setup -----------------------
    int Hpix = std::max(2, (int)std::lround(Wpix * (ryb - rya) / (rxb - rxa)));
    if (Hpix % 2) ++Hpix; if (Wpix % 2) ++Wpix;
    auto inDomain = [=](double x, double y) {
        return std::hypot(x - cx, y - cy) > radius * 1.001 &&
               x > xa - 1e-9 && x < xb + 1e-9 && y > ya - 1e-9 && y < yb + 1e-9;
    };
    // Per-output channel: a kind, an ouput dir, and (if a streakline channel)
    // its own particle population.  Each channel gets its own ffmpeg invocation.
    struct OutChannel {
        string kind;          // "vorticity" | "speed" | "flow"
        string dir;
        string mp4;           // suggested output filename
        std::unique_ptr<ParticleTracer> tracer;   // only set when kind=="flow"
    };
    std::vector<OutChannel> chans;
    chans.reserve(outputs.size());
    bool needLocator = false;
    int flowIdx = 0;
    for (const auto& kind : outputs) {
        OutChannel ch; ch.kind = kind;
        if (kind == "vorticity") { ch.dir = framesDirVort; ch.mp4 = "cylinder_vortex.mp4"; }
        else if (kind == "speed"){ ch.dir = framesDirVort; ch.mp4 = "cylinder_speed.mp4"; }
        else /* flow */          {
            // multiple flow channels would clobber each other -- name them apart
            ch.dir = framesDirFlow + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx)));
            ch.mp4 = "cylinder_flow"   + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx))) + ".mp4";
            ++flowIdx;
            needLocator = true;
            ch.tracer = std::make_unique<ParticleTracer>();
            ch.tracer->N = std::max(1, nParticles);
            ch.tracer->xa = rxa; ch.tracer->xb = rxb; ch.tracer->ya = rya; ch.tracer->yb = ryb;
            ch.tracer->inflowX = std::max(rxa, xa) + 0.02 * (rxb - rxa);
            ch.tracer->inDomain = inDomain;
            ch.tracer->rng = particleSeed ? particleSeed : 12345u;
            ch.tracer->trailLen    = std::max(1, trailLen);
            ch.tracer->trailStride = std::max(1, trailStride);
            ch.tracer->reset();
        }
        fs::create_directories(ch.dir);
        for (const auto& e : fs::directory_iterator(ch.dir))
            if (e.path().extension() == ".ppm") fs::remove(e.path());
        chans.push_back(std::move(ch));
    }
    MeshLocator locator;
    if (needLocator) locator.build(mesh);

    auto writeAllFrames = [&](int frameIdx) {
        char fn[512];
        for (auto& ch : chans) {
            snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", ch.dir.c_str(), frameIdx);
            if (ch.kind == "vorticity") {
                writeFieldPPM(fn, fem, mesh, elem2dof, integ.vorticity(), Wpix, Hpix,
                              rxa, rxb, rya, ryb, -vortClip, vortClip, CM_COOLWARM, inDomain);
            } else if (ch.kind == "speed") {
                writeFieldPPM(fn, fem, mesh, elem2dof, integ.speed(), Wpix, Hpix,
                              rxa, rxb, rya, ryb, 0.0, 1.6 * Uinf, CM_VIRIDIS, inDomain);
            } else { // flow
                writeParticlesPPM(fn, fem, mesh, elem2dof, integ.speed(), Wpix, Hpix,
                                  rxa, rxb, rya, ryb, 0.0, 1.6 * Uinf, *ch.tracer, inDomain, bgDim);
            }
        }
    };
    auto advanceAllTracers = [&](double dtStep) {
        for (auto& ch : chans)
            if (ch.tracer)
                ch.tracer->advance(fem, mesh, elem2dof, integ.u(), integ.v(), locator, dtStep);
    };

    ofstream hist("ns_forces.csv");
    hist << "t,CD,CL\n";
    std::vector<double> tCL, CLval, CDval;     // for Strouhal / mean-drag estimation

    // ----------------------- time loop -----------------------
    cout << "\nEvolving (explicit LF convection, BDF2/EX2, "
         << (pressureMode == NSPRESSURE_DIRECT_PPE ? "direct PPE pressure" : "projection pressure")
         << "; matrices factorised once)...\n";
    auto t0wall = chrono::high_resolution_clock::now();
    int frame = 0;
    writeAllFrames(frame);
    ++frame;
    const double qdyn = 0.5 * Uinf * Uinf * D;     // dynamic pressure * D (rho = 1)
    for (int s = 1; s <= nsteps; ++s) {
        if (!integ.step(s * dt)) { cout << "ERROR: solve failed at step " << s << "\n"; return 1; }
        advanceAllTracers(dt);
        double Fx, Fy; integ.cylinderForce(isCyl, Fx, Fy);
        double CD = Fx / qdyn, CL = Fy / qdyn;
        double t = s * dt;
        hist << t << "," << CD << "," << CL << "\n";
        tCL.push_back(t); CLval.push_back(CL); CDval.push_back(CD);
        if (s % save_every == 0 || s == nsteps) {
            writeAllFrames(frame);
            ++frame;
            double el = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
            cout << "  step " << setw(6) << s << "/" << nsteps << "  t=" << fixed << setprecision(2) << t
                 << "  CD=" << setprecision(3) << CD << "  CL=" << setw(7) << CL
                 << "  frame=" << setw(4) << (frame-1) << "  (" << setprecision(1) << el << "s)\n";
        }
    }
    hist.close();

    // ----------------------- Strouhal from lift zero-crossings (latter half) -----------------------
    std::vector<double> crossings;
    for (size_t i = 1; i < CLval.size(); ++i)
        if (tCL[i] > 0.5 * t_end && CLval[i-1] < 0 && CLval[i] >= 0) {
            double frac = -CLval[i-1] / (CLval[i] - CLval[i-1]);
            crossings.push_back(tCL[i-1] + frac * (tCL[i] - tCL[i-1]));
        }
    double St = 0.0, period = 0.0;
    if (crossings.size() >= 2) {
        period = (crossings.back() - crossings.front()) / (crossings.size() - 1);
        St = D / (Uinf * period);
    }
    // mean drag and lift amplitude over the statistically-stationary latter half
    double sumCD = 0, maxCL = -1e30, minCL = 1e30; int cnt = 0;
    for (size_t i = 0; i < CLval.size(); ++i) if (tCL[i] > 0.5 * t_end) {
        sumCD += CDval[i]; maxCL = std::max(maxCL, CLval[i]); minCL = std::min(minCL, CLval[i]); ++cnt;
    }
    double meanCD = cnt ? sumCD / cnt : 0.0;
    double wall = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
    cout << "\nDone. " << frame << " frames in";
    for (size_t i = 0; i < chans.size(); ++i)
        cout << (i ? ", " : " ") << "'" << chans[i].dir << "/'";
    cout << ".  wall=" << fixed << setprecision(1) << wall << "s\n";
    cout << "  mean CD = " << setprecision(3) << meanCD;
    if (St > 0) cout << ",  Strouhal St = " << setprecision(4) << St
                     << "  (period " << setprecision(3) << period << ", CL in ["
                     << minCL << ", " << maxCL << "])\n";
    else cout << "  (not enough lift oscillations for a Strouhal estimate; increase t_end)\n";
    cout << "  force history -> ns_forces.csv\n";
    cout << "\nAssemble the movies:\n";
    for (const auto& ch : chans)
        cout << "  ffmpeg -y -framerate 25 -i " << ch.dir << "/frame_%05d.ppm"
             << " -c:v libx264 -pix_fmt yuv420p -crf 18 " << ch.mp4 << "\n";
    return 0;
}
