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
#include "IBCoupler.h"               // overlayPolylineOnPPM
#include "Spoon.h"
#include "Json.h"

using namespace Eigen;
using namespace std;
using namespace ns;
namespace fs = std::filesystem;

// ===========================================================================
// "A spoon stirs a bowl of soup."  A closed circular bowl (no-slip wall all
// the way round) is filled with quiescent fluid.  A rigid spoon blade,
// submerged near the rim, is swept rigidly about the bowl centre through ONE
// arc-stroke and then withdrawn.  Being broadside to its own motion the blade
// acts as a paddle: each tip sheds a shear layer that rolls up into a vortex,
// and the two counter-rotating vortices pair into a self-propelled dipole that
// shoots across the bowl after the stroke and eventually collides with the
// curved wall.  ("Vortices on both sides of the spoon" -- the impulsively
// started/stopped flat-plate dipole, confined.)
//
// The blade motion is PRESCRIBED and imposed on the DG incompressible
// Navier-Stokes flow by a direct-forcing immersed boundary
//   F = (alpha/dt) * gain(t) * (V_blade - u)
// scattered at a marker cloud over the blade footprint (Spoon::applyForcing).
// The flow solver is the same high-order splitting NS integrator as the
// cylinder/von-Karman demo; only the domain (closed disk, all-Neumann pressure
// regularised by pinning one rim edge) and the forcing differ.
// ===========================================================================

namespace {

// Integrate a scalar DG field and its square over the mesh:  I = \int f,
// I2 = \int f^2  (used for total circulation / enstrophy diagnostics).
void integrateField(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                    const VectorXd& f, double& I, double& I2) {
    MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w);
    int nq = (int)w.size(), NT = mesh.elem.rows(), locDof = fem.locDof;
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = fem.computeBasisValue_all(quadL.row(q)).row(0);
    I = 0.0; I2 = 0.0;
    for (int e = 0; e < NT; ++e) {
        double area = fem.area(e);
        for (int q = 0; q < nq; ++q) {
            double fv = 0.0;
            for (int i = 0; i < locDof; ++i) fv += phi[q](i) * f(elem2dof(e, i));
            I  += w(q) * area * fv;
            I2 += w(q) * area * fv * fv;
        }
    }
}

} // anon

int main(int argc, char** argv) {
    // ----------------------- flow / domain params ----------------------------
    int    ord    = 2;
    double R      = 1.0;            // bowl radius
    double h_     = 0.025;          // target element size in the stirring core
    double farRatio = 1.9;          // coarse/fine size ratio away from the refine band
    double grade  = 0.2;            // size growth rate with distance from the band
    double refineInMargin  = 0.13;  // adaptive band extends this far inward of the blade
    double refineOutMargin = 0.25;  // ... and this far outward (the shed-dipole corridor)
    double Re     = 800.0;          // Reynolds number on blade length & stroke speed
    double cfl    = 0.40;
    double dt     = 0.0;
    double t_end  = 12.0;
    int    save_every = 0;
    int    n_frames = 400;
    double sigmaFac = 8.0;
    double gradDiv = 0.0;
    double ppeDivDamping = 30.0;
    string pressureModeName = "direct_ppe";
    // semi-implicit immersed-boundary constraint
    double ib_eps   = 1e-6;   // Tikhonov regularisation of the Schur matrix (0 = exact)
    int    ib_maxcg = 300;    // max CG iterations for the marker-force Schur solve
    double ib_tol   = 1e-6;   // CG relative tolerance
    int    ib_subiters = 0;   // EXPERIMENTAL: (constraint <-> projection) rounds
    int    ib_end_project = 0; // 1 = end with projection (div-free), 0 = end with constraint
    double ib_kernel_fac = 0.0; // mollified-transfer kernel radius, in units of h (0 = pointwise)
    double filter_strength = 0.0;   // modal top-mode damping in [0,1) (0 = off)
    double filter_sensor_lo = 0.02; // hi-mode energy fraction: below -> no filtering
    double filter_sensor_hi = 0.08; // hi-mode energy fraction: at/above -> full strength
    double av_beta = 0.0;           // artificial-viscosity coeff in units of h^2 (0 = off)
    double av_sensor_lo = 0.15;     // checkerboard sensor: below -> no AV
    double av_sensor_hi = 0.50;     // checkerboard sensor: at/above -> full AV blend
    // After the blade withdraws, free decay has no thin shear layers to preserve, so
    // apply a small UNGATED global AV floor that continuously damps residual grid-
    // scale content the sensor misses (a no-op on the smooth resolved vortices).
    double av_global_decay = 0.0;   // global AV blend floor in free decay (0 = off)
    double hv_fac = 0.0;            // global hyperviscosity coeff in units of h^4 (0 = off)
    // rendering
    int    Wpix = 900;
    double margin = 0.04;           // render window margin (fraction of R)
    double vortClip = 12.0;
    double divClip  = 5.0;
    std::vector<string> outputs;
    string framesDirVort  = "ns_spoon_frames";
    string framesDirFlow  = "ns_spoon_flow_frames";
    int    nParticles     = 2600;
    unsigned int particleSeed = 12345u;
    int    trailLen       = 26;
    int    trailStride    = 2;
    double bgDim          = 0.14;
    bool   drawRim        = true;

    // ----------------------- spoon params ------------------------------------
    Spoon spoon;   // defaults live in Spoon.h

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("spoon_config.json")) cfgPath = "spoon_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        ord = cfg.getInt("ord", ord);
        R = cfg.getNumber("R", R);
        h_ = cfg.getNumber("h", h_);
        farRatio = cfg.getNumber("far_ratio", farRatio);
        grade = cfg.getNumber("grade", grade);
        refineInMargin = cfg.getNumber("refine_in_margin", refineInMargin);
        refineOutMargin = cfg.getNumber("refine_out_margin", refineOutMargin);
        Re = cfg.getNumber("Re", Re);
        cfl = cfg.getNumber("cfl", cfl);
        dt = cfg.getNumber("dt", dt);
        t_end = cfg.getNumber("t_end", t_end);
        save_every = cfg.getInt("save_every", save_every);
        n_frames = cfg.getInt("n_frames", n_frames);
        sigmaFac = cfg.getNumber("sigma_fac", sigmaFac);
        gradDiv = cfg.getNumber("grad_div", gradDiv);
        ppeDivDamping = cfg.getNumber("ppe_div_damping", ppeDivDamping);
        pressureModeName = cfg.getString("pressure_mode", pressureModeName);
        ib_eps   = cfg.getNumber("ib_eps", ib_eps);
        ib_maxcg = cfg.getInt("ib_maxcg", ib_maxcg);
        ib_tol   = cfg.getNumber("ib_tol", ib_tol);
        ib_subiters = cfg.getInt("ib_subiters", ib_subiters);
        ib_end_project = cfg.getInt("ib_end_project", ib_end_project);
        ib_kernel_fac = cfg.getNumber("ib_kernel_fac", ib_kernel_fac);
        filter_strength = cfg.getNumber("filter_strength", filter_strength);
        filter_sensor_lo = cfg.getNumber("filter_sensor_lo", filter_sensor_lo);
        filter_sensor_hi = cfg.getNumber("filter_sensor_hi", filter_sensor_hi);
        av_beta = cfg.getNumber("av_beta", av_beta);
        av_sensor_lo = cfg.getNumber("av_sensor_lo", av_sensor_lo);
        av_sensor_hi = cfg.getNumber("av_sensor_hi", av_sensor_hi);
        av_global_decay = cfg.getNumber("av_global_decay", av_global_decay);
        hv_fac = cfg.getNumber("hv_fac", hv_fac);
        Wpix = cfg.getInt("Wpix", Wpix);
        margin = cfg.getNumber("render_margin", margin);
        vortClip = cfg.getNumber("vort_clip", vortClip);
        divClip  = cfg.getNumber("div_clip", divClip);
        if (cfg.contains("outputs") && cfg.obj.at("outputs").type == cfgjson::Json::Array) {
            for (const auto& it : cfg.obj.at("outputs").arr)
                if (it.type == cfgjson::Json::String) outputs.push_back(it.str);
        }
        framesDirVort = cfg.getString("frames_dir", framesDirVort);
        framesDirFlow = cfg.getString("flow_frames_dir", framesDirFlow);
        nParticles    = cfg.getInt("n_particles", nParticles);
        particleSeed  = (unsigned int)cfg.getInt("particle_seed", (int)particleSeed);
        trailLen      = cfg.getInt("trail_len", trailLen);
        trailStride   = cfg.getInt("trail_stride", trailStride);
        bgDim         = cfg.getNumber("bg_dim", bgDim);
        drawRim       = cfg.getInt("draw_rim", drawRim ? 1 : 0) != 0;
        // spoon geometry / schedule / forcing
        spoon.aRad    = cfg.getNumber("spoon_aRad", spoon.aRad);
        spoon.bThk    = cfg.getNumber("spoon_bThk", spoon.bThk);
        spoon.dmark   = cfg.getNumber("spoon_dmark", spoon.dmark);
        {
            string sh = cfg.getString("spoon_shape", spoon.crescent ? "crescent" : "ellipse");
            for (auto& c : sh) c = (char)std::tolower((unsigned char)c);
            spoon.crescent = (sh != "ellipse");
        }
        spoon.crescentR2frac = cfg.getNumber("spoon_crescent_r2frac", spoon.crescentR2frac);
        spoon.crescentOffset = cfg.getNumber("spoon_crescent_offset", spoon.crescentOffset);
        spoon.px      = cfg.getNumber("spoon_px", spoon.px);
        spoon.py      = cfg.getNumber("spoon_py", spoon.py);
        spoon.rPivot  = cfg.getNumber("spoon_rPivot", spoon.rPivot);
        spoon.bearing0= cfg.getNumber("spoon_bearing0", spoon.bearing0);
        spoon.tEnter  = cfg.getNumber("spoon_tEnter", spoon.tEnter);
        spoon.tStroke = cfg.getNumber("spoon_tStroke", spoon.tStroke);
        spoon.tRamp   = cfg.getNumber("spoon_tRamp", spoon.tRamp);
        spoon.tLift   = cfg.getNumber("spoon_tLift", spoon.tLift);
        spoon.sweep   = cfg.getNumber("spoon_sweep", spoon.sweep);
        spoon.ibAlpha = cfg.getNumber("spoon_ibAlpha", spoon.ibAlpha);
        spoon.ibForceCap = cfg.getNumber("spoon_ibForceCap", spoon.ibForceCap);
    }
    if (outputs.empty()) outputs = {"vorticity", "flow"};
    for (auto& s : outputs) {
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s != "vorticity" && s != "speed" && s != "flow" && s != "divergence") {
            cout << "ERROR: outputs entry '" << s << "' must be 'vorticity', 'speed', 'flow', or 'divergence'\n";
            return 1;
        }
    }

    spoon.build();

    // Characteristic scales: blade length L = 2*aRad, stroke speed U = omax*rPivot.
    const double L_blade = 2.0 * spoon.aRad;
    const double omax  = spoon.peakOmega();              // peak of the parabolic stroke
    const double U_ref = std::abs(omax) * spoon.rPivot;
    const double strokeDist = std::abs(spoon.sweep) * spoon.rPivot;   // blade-centre travel
    const double U_tip = std::abs(omax) * (spoon.rPivot + spoon.aRad);
    const double nu = U_ref * L_blade / Re;
    const double Umax = std::max(1.0, U_tip);   // for CFL
    if (dt <= 0.0) dt = cfl * h_ / ((2.0 * ord + 1.0) * Umax);
    int nsteps = std::max(1, (int)std::lround(t_end / dt));
    if (save_every <= 0) save_every = std::max(1, nsteps / std::max(1, n_frames));
    double sigma = sigmaFac * (ord + 1) * (ord + 1);
    int pressureMode = (pressureModeName == "projection") ? NSPRESSURE_PROJECTION : NSPRESSURE_DIRECT_PPE;

    cout << "DG NS -- spoon stirring a bowl of soup\n";
    cout << "  config:  " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  bowl:    R=" << R << "  dP" << ord << "  h=" << h_ << "  far_ratio=" << farRatio << "\n";
    cout << "  spoon:   shape=" << (spoon.crescent ? "crescent" : "ellipse")
         << "  L=" << L_blade << " thick=" << 2*spoon.bThk
         << "  pivot=(" << spoon.px << "," << spoon.py << ") rPivot=" << spoon.rPivot
         << "  sweep=" << spoon.sweep << " rad  bearing0=" << spoon.bearing0 << "\n";
    cout << "  stroke:  tEnter=" << spoon.tEnter << " tStroke=" << spoon.tStroke
         << " tLift=" << spoon.tLift << "  (quadratic v-profile)"
         << "  travel=" << strokeDist << " (=" << strokeDist / L_blade << " blade-lengths)\n";
    cout << "  speed:   omega_max=" << omax << "  U_ref=" << U_ref << " U_tip=" << U_tip << "\n";
    cout << "  Re=" << Re << "  nu=" << nu
         << "  IB=semi-implicit constraint (eps=" << ib_eps << ")  markers=" << spoon.M() << "\n";
    cout << "  dt=" << dt << "  t_end=" << t_end << "  steps=" << nsteps
         << "  save_every=" << save_every << "  ppe_div_damping=" << ppeDivDamping << "\n";

    // ----------------------- mesh + DG space ---------------------------------
    // Adaptive refinement: a fine annular band straddling the spoon's swept path
    // (blade radial extent +/- margins, the outward margin covering where the shed
    // dipole travels), coarsening into the quiescent core and rim.
    BowlGeom geom{0.0, 0.0, R};
    double bandLo = std::max(0.0, spoon.rPivot - spoon.aRad - refineInMargin);
    double bandHi = std::min(R,   spoon.rPivot + spoon.aRad + refineOutMargin);
    cout << "  mesh:    adaptive band r in [" << bandLo << "," << bandHi << "]"
         << "  h_fine=" << h_ << " h_coarse=" << h_ * farRatio << "\n";
    Mesh mesh;
    generateBowlMesh(mesh, geom, h_, farRatio, grade, bandLo, bandHi);
    FEM fem(ord, mesh);
    MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    MatrixXi edge, edge2side; mesh.getEdge2Side(edge, edge2side);
    VectorXi tag = classifyBowlEdges(mesh, edge, edge2side, geom);
    cout << "  nDof=" << nDof << "\n";

    // ----------------------- BCs: closed no-slip bowl ------------------------
    // All boundary edges: no-slip velocity (Dirichlet u=v=0) + high-order
    // pressure Neumann.  The all-Neumann PPE is singular (pressure defined up to
    // a constant); pin ONE rim edge to Dirichlet p=0 to fix the gauge.  The
    // pinned edge is chosen at the bottom of the bowl, far from the stroke.
    int NE = edge.rows();
    BCData bc;
    bc.bcU = VectorXi::Zero(NE); bc.bcV = VectorXi::Zero(NE); bc.bcP = VectorXi::Zero(NE);
    VectorXi isWall = VectorXi::Zero(NE);
    int pinEdge = -1; double pinY = 1e300;
    for (int e = 0; e < NE; ++e) {
        if (tag(e) != BD_BOWL_WALL) continue;
        bc.bcU(e) = BCN_DIRICHLET; bc.bcV(e) = BCN_DIRICHLET; bc.bcP(e) = BCN_NEUMANN_HO;
        isWall(e) = 1;
        Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1)));
        if (m.y() < pinY) { pinY = m.y(); pinEdge = e; }
    }
    if (pinEdge >= 0) bc.bcP(pinEdge) = BCN_DIRICHLET;
    cout << "  wall edges=" << isWall.sum() << "  pressure pin edge=" << pinEdge
         << " at y=" << pinY << "\n";

    bc.velDir = [](double, double, double) { return Vector2d(0.0, 0.0); };
    bc.velAcc = [](double, double, double) { return Vector2d(0.0, 0.0); };
    bc.presDir = [](double, double, double) { return 0.0; };

    NSIntegrator integ(fem, mesh, elem2dof, edge, edge2side, bc, nu, dt, sigma, 1.0,
                       gradDiv, pressureMode, ppeDivDamping);
    integ.setInitial(VectorXd::Zero(nDof), VectorXd::Zero(nDof));   // quiescent soup
    integ.ibSubIters = std::max(0, ib_subiters);
    integ.ibEndWithProjection = (ib_end_project != 0);
    if (integ.ibSubIters > 0)
        cout << "  IB sub-iterations: " << integ.ibSubIters
             << "  (end with " << (integ.ibEndWithProjection ? "projection" : "constraint") << ")\n";
    if (ib_kernel_fac > 0.0)
        cout << "  IB mollified transfer: Wendland C2, delta=" << ib_kernel_fac << "*h = "
             << ib_kernel_fac * h_ << "\n";
    integ.filterStrength = std::max(0.0, filter_strength);
    integ.filterSensorLo = filter_sensor_lo;
    integ.filterSensorHi = filter_sensor_hi;
    if (integ.filterStrength > 0.0)
        cout << "  modal filter: strength=" << integ.filterStrength
             << " sensor=[" << integ.filterSensorLo << "," << integ.filterSensorHi << "]"
             << " (top-degree modes, sensor-gated)\n";
    integ.avBeta = std::max(0.0, av_beta);
    integ.avSensorLo = av_sensor_lo;
    integ.avSensorHi = av_sensor_hi;
    if (integ.avBeta > 0.0)
        cout << "  artificial viscosity: beta=" << integ.avBeta << "*h^2"
             << " sensor=[" << integ.avSensorLo << "," << integ.avSensorHi << "]"
             << " (inter-element checkerboard, sensor-gated)\n";
    integ.hvFac = std::max(0.0, hv_fac);
    if (integ.hvFac > 0.0)
        cout << "  hyperviscosity: beta=" << integ.hvFac << "*h^4 (global biharmonic grid-scale filter)\n";

    // ----------------------- IB locator + load buffers -----------------------
    MeshLocator locatorIB; locatorIB.build(mesh);

    // ----------------------- output channels ---------------------------------
    double win = R * (1.0 + margin);
    double rxa = -win, rxb = win, rya = -win, ryb = win;
    int Hpix = std::max(2, (int)std::lround(Wpix * (ryb - rya) / (rxb - rxa)));
    if (Hpix % 2) ++Hpix; if (Wpix % 2) ++Wpix;
    auto inDomain = [=](double x, double y) {
        return std::hypot(x, y) < R * 0.999;
    };
    struct OutChannel {
        string kind, dir, mp4;
        std::unique_ptr<ParticleTracer> tracer;
    };
    std::vector<OutChannel> chans;
    chans.reserve(outputs.size());
    bool needLocator = false;
    int flowIdx = 0;
    for (const auto& kind : outputs) {
        OutChannel ch; ch.kind = kind;
        if (kind == "vorticity") { ch.dir = framesDirVort; ch.mp4 = "spoon_vortex.mp4"; }
        else if (kind == "speed"){ ch.dir = framesDirVort; ch.mp4 = "spoon_speed.mp4"; }
        else if (kind == "divergence") { ch.dir = framesDirVort + "_div"; ch.mp4 = "spoon_div.mp4"; }
        else /* flow */          {
            ch.dir = framesDirFlow + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx)));
            ch.mp4 = "spoon_flow"   + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx))) + ".mp4";
            ++flowIdx;
            needLocator = true;
            ch.tracer = std::make_unique<ParticleTracer>();
            ch.tracer->N = std::max(1, nParticles);
            ch.tracer->xa = rxa; ch.tracer->xb = rxb; ch.tracer->ya = rya; ch.tracer->yb = ryb;
            ch.tracer->inflowX = rxa + 0.5 * (rxb - rxa);   // closed bowl: reseed near centre
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
    MeshLocator locatorTracer;
    if (needLocator) locatorTracer.build(mesh);

    // Bowl rim circle (closed polyline) for overlay.
    MatrixXd rim(73, 2);
    for (int i = 0; i < 73; ++i) {
        double th = 2.0 * M_PI * i / 72.0;
        rim(i, 0) = R * std::cos(th); rim(i, 1) = R * std::sin(th);
    }

    double speedClip = 1.2 * U_tip;
    auto writeAllFrames = [&](int frameIdx) {
        char fn[512];
        MatrixXd bladePoly = spoon.outline(56);
        // For vorticity/speed rendering: mask the blade interior (don't render field there).
        auto inDomainMasked = [&](double x, double y) {
            if (std::hypot(x, y) >= R * 0.999) return false;
            if (spoon.submerged() && spoon.contains(x, y)) return false;
            return true;
        };
        for (auto& ch : chans) {
            snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", ch.dir.c_str(), frameIdx);
            if (ch.kind == "vorticity") {
                writeFieldPPM(fn, fem, mesh, elem2dof, integ.vorticity(), Wpix, Hpix,
                              rxa, rxb, rya, ryb, -vortClip, vortClip, CM_COOLWARM, inDomainMasked);
            } else if (ch.kind == "speed") {
                writeFieldPPM(fn, fem, mesh, elem2dof, integ.speed(), Wpix, Hpix,
                              rxa, rxb, rya, ryb, 0.0, speedClip, CM_VIRIDIS, inDomainMasked);
            } else if (ch.kind == "divergence") {
                writeFieldPPM(fn, fem, mesh, elem2dof, integ.divergence(), Wpix, Hpix,
                              rxa, rxb, rya, ryb, -divClip, divClip, CM_COOLWARM, inDomainMasked);
            } else {
                writeParticlesPPM(fn, fem, mesh, elem2dof, integ.speed(), Wpix, Hpix,
                                  rxa, rxb, rya, ryb, 0.0, speedClip, *ch.tracer, inDomain, bgDim);
            }
            int rimRGB  = (ch.kind == "vorticity") ? 0x303030 : 0x606060;
            int bladeRGB = (ch.kind == "vorticity") ? 0x111111 : 0xFFFFFF;
            if (drawRim) overlayPolylineOnPPM(string(fn), rim, rxa, rxb, rya, ryb, 2, rimRGB);
            if (spoon.submerged())
                overlayPolylineOnPPM(string(fn), bladePoly, rxa, rxb, rya, ryb, 2, bladeRGB);
        }
    };
    auto advanceAllTracers = [&](double dtStep) {
        for (auto& ch : chans)
            if (ch.tracer)
                ch.tracer->advance(fem, mesh, elem2dof, integ.u(), integ.v(), locatorTracer, dtStep);
    };

    // ----------------------- diagnostics file --------------------------------
    ofstream hist("spoon_diagnostics.csv");
    hist << "t,phi,omega,gain,ib_cg_iters,ib_resid,KE,enstrophy,circulation\n";

    // Returns (kinetic energy, enstrophy = int omega^2, total circulation = int omega).
    auto flowDiagnostics = [&]() {
        double cu, cu2, cv, cv2;
        integrateField(fem, mesh, elem2dof, integ.u(), cu, cu2);
        integrateField(fem, mesh, elem2dof, integ.v(), cv, cv2);
        double cw, cw2;
        integrateField(fem, mesh, elem2dof, integ.vorticity(), cw, cw2);
        return std::make_tuple(0.5 * (cu2 + cv2), cw2, cw);
    };

    // ----------------------- time loop ---------------------------------------
    cout << "\nStirring the soup...\n";
    auto t0wall = chrono::high_resolution_clock::now();
    int frame = 0;
    spoon.update(0.0);
    writeAllFrames(frame); ++frame;

    // Adaptive time-stepping parameters.
    const double cfl_max    = cfl * 1.4;           // trigger dt reduction (aggressive)
    const double cfl_min    = cfl * 0.3;           // trigger dt growth
    const double dt_max     = dt;                  // never exceed the initial CFL-based dt
    const double dt_min     = dt * 0.05;           // floor (avoid near-zero dt)
    const double hCFL       = h_ / (2.0 * ord + 1.0);  // CFL length scale
    double dt_cur = dt;
    double t_cur  = 0.0;
    double t_next_frame = t_end / n_frames;        // next time to write a frame
    int    totalSteps = 0;
    bool   avBumped = false;                        // free-decay AV switch done?

    // Compute max speed in the fluid domain, EXCLUDING elements whose centroid
    // lies inside the immersed boundary.  This prevents the IB-interior artefacts
    // (which are always large) from driving the time step to zero.
    int NT = mesh.elem.rows();
    int locDof = fem.locDof;
    auto maxFluidSpeed = [&]() {
        const VectorXd& uu = integ.u();
        const VectorXd& vv = integ.v();
        double smax = 0.0;
        for (int e = 0; e < NT; ++e) {
            // Skip elements inside the spoon blade.
            if (spoon.submerged()) {
                double cx = 0, cy = 0;
                for (int k = 0; k < 3; ++k) {
                    cx += mesh.node(mesh.elem(e, k), 0);
                    cy += mesh.node(mesh.elem(e, k), 1);
                }
                cx /= 3.0; cy /= 3.0;
                if (spoon.contains(cx, cy)) continue;
            }
            // Max speed over DOFs in this element.
            for (int j = 0; j < locDof; ++j) {
                int dof = elem2dof(e, j);
                double sp = std::hypot(uu(dof), vv(dof));
                if (sp > smax) smax = sp;
            }
        }
        return smax;
    };

    while (t_cur < t_end - 1e-14) {
        // --- Pre-step CFL check: use last step's max speed to limit dt BEFORE stepping ---
        if (totalSteps > 0) {
            double Upre = maxFluidSpeed();
            double cfl_pre = Upre * dt_cur / hCFL;
            if (cfl_pre > cfl_max) {
                double dt_need = cfl * hCFL / std::max(Upre, 1e-12);
                dt_cur = std::max(dt_min, std::min(dt_cur, dt_need));
                integ.setDt(dt_cur);
            }
        }

        // Clamp dt so we don't overshoot t_end.
        if (t_cur + dt_cur > t_end) dt_cur = t_end - t_cur;
        double t_new = t_cur + dt_cur;
        ++totalSteps;

        // 1) update the prescribed spoon pose at the new time level.
        spoon.update(t_new);
        // 2) arm the SEMI-IMPLICIT immersed-boundary no-slip constraint.
        if (spoon.active(t_new))
            integ.setIBConstraint(spoon.X, spoon.V, locatorIB, ib_eps, ib_maxcg, ib_tol,
                                  ib_kernel_fac * h_);
        else
            integ.clearIBConstraint();
        // 2b) once the blade has withdrawn, enable the ungated global AV floor
        //     (free decay has no thin shear layers to preserve) to mop up the
        //     faint residual grid-scale checkerboard the sensor misses.
        if (av_beta > 0.0 && av_global_decay > 0.0 && !avBumped && !spoon.active(t_new)) {
            integ.avGlobalFloor = av_global_decay;
            avBumped = true;
            cout << fixed << setprecision(3)
                 << "  [t=" << t_new << "] blade withdrawn -> free-decay global AV floor="
                 << integ.avGlobalFloor << "\n" << flush;
        }
        // 3) advance the fluid.
        if (!integ.step(t_new)) {
            cout << "ERROR: NS solve failed at step " << totalSteps
                 << " t=" << t_new << "\n"; return 1;
        }
        // 4) tracers advance on the new field.
        if (needLocator) advanceAllTracers(dt_cur);

        t_cur = t_new;

        // diagnostics
        auto [KE, enstrophy, circ] = flowDiagnostics();
        if (!std::isfinite(KE) || KE > 1e6) {
            cout << "ERROR: blow-up detected at step " << totalSteps
                 << " t=" << t_cur << " KE=" << KE << "\n";
            return 1;
        }
        int   ibIt  = spoon.active(t_cur) ? integ.ibIters()    : 0;
        double ibRes = spoon.active(t_cur) ? integ.ibResidual() : 0.0;
        hist << fixed << setprecision(6)
             << t_cur << "," << spoon.phi << "," << spoon.omega << "," << spoon.gain << ","
             << ibIt << "," << ibRes << "," << KE << "," << enstrophy << "," << circ << "\n";

        // --- Adaptive CFL control (exclude IB interior) ---
        double Ufluid = maxFluidSpeed();
        double cfl_now = Ufluid * dt_cur / hCFL;
        bool dtChanged = false;
        if (cfl_now > cfl_max && dt_cur > dt_min * 1.01) {
            // Shrink: halve dt (but never below dt_min).
            dt_cur = std::max(dt_min, dt_cur * 0.5);
            integ.setDt(dt_cur);
            dtChanged = true;
        } else if (cfl_now < cfl_min && dt_cur < dt_max * 0.99) {
            // Grow: increase by 25% (conservative growth to avoid oscillation).
            dt_cur = std::min(dt_max, dt_cur * 1.25);
            integ.setDt(dt_cur);
            dtChanged = true;
        }

        // --- Frame output at uniform time intervals ---
        if (t_cur >= t_next_frame - 1e-14 || t_cur >= t_end - 1e-14) {
            writeAllFrames(frame); ++frame;
            double el = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
            cout << "  step " << setw(6) << totalSteps
                 << "  t=" << fixed << setprecision(3) << t_cur
                 << "  dt=" << scientific << setprecision(2) << dt_cur
                 << "  CFL=" << fixed << setprecision(2) << cfl_now
                 << "  phi=" << setprecision(2) << spoon.phi
                 << "  ibCG=" << ibIt << " res=" << scientific << setprecision(1) << ibRes
                 << "  KE=" << setprecision(2) << KE
                 << fixed << setprecision(1) << "  (" << el << "s)"
                 << (dtChanged ? "  *dt" : "") << "\n" << flush;
            t_next_frame += t_end / n_frames;
        }
    }
    hist.close();

    double wall = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
    cout << "\nDone. " << frame << " frames, " << totalSteps << " steps in";
    for (size_t i = 0; i < chans.size(); ++i)
        cout << (i ? ", " : " ") << "'" << chans[i].dir << "/'";
    cout << ".  wall=" << fixed << setprecision(1) << wall << "s\n";
    cout << "  diagnostics -> spoon_diagnostics.csv\n";
    cout << "\nAssemble the movies:\n";
    for (const auto& ch : chans)
        cout << "  ffmpeg -y -framerate 25 -i " << ch.dir << "/frame_%05d.ppm"
             << " -c:v libx264 -pix_fmt yuv420p -crf 18 " << ch.mp4 << "\n";
    return 0;
}
