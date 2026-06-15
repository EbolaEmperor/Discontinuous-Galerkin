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
#include "CosseratFilament.h"
#include "IBCoupler.h"
#include "Json.h"

using namespace Eigen;
using namespace std;
using namespace ns;
namespace fs = std::filesystem;

// ===========================================================================
// FSI demo: cylinder + flexible filament tail.
//
// Flow (DG dP_k, BDF2/EX2 + direct PPE, identical to navier_stokes_cylinder)
// past a circular cylinder, with a 1-D Cosserat rod attached at the rear
// stagnation point.  The rod is advected by the fluid through an immersed-
// boundary coupling (no remeshing).
//
//   - mesh→rod : sample (u,v) at every rod node (DG basis as kernel)
//   - rod→mesh : scatter direct-forcing constraint  F = alpha*(u_target - u_h(X))
//                back onto the DG nodal load.  alpha = m / dt^2 in the spirit
//                of Uhlmann 2005 / Bhalla 2013.
//   - structure: implicit Newmark-beta on the discrete elastic rod
//
// The pipeline is:
//   1. predict rod from t^n with the t^n fluid load on it,
//   2. compute alpha*(rod velocity - sampled fluid velocity) at every node,
//      scatter as body force,
//   3. NSIntegrator::stepWithBodyForce -> u^{n+1}, p^{n+1}.
// ===========================================================================

int main(int argc, char** argv) {
    // ----------------------- flow parameters (same defaults as ns_main) -----
    int    ord    = 2;
    double radius = 0.5;
    double cx_ = 0.0, cy_ = 0.0;
    double xa = -5.0, xb = 18.0;
    double ya = -6.0, yb = 6.0;
    double h      = 0.07;
    double farRatio = 9.0;
    double grade  = 0.22;
    double Re     = 100.0;
    double Uinf   = 1.0;
    double cfl    = 0.5;
    double dt     = 0.0;
    double t_end  = 90.0;
    int    save_every = 0;
    int    n_frames = 240;
    double sigmaFac = 8.0;
    double gradDiv = 0.0;
    double ppeDivDamping = 10.0;
    string pressureModeName = "direct_ppe";
    double rampTime = 0.5;
    double perturb = 0.30;
    // rendering
    int    Wpix = 1100;
    double rxa = -3.0, rxb = 16.0;
    double rya = -4.0, ryb = 4.0;
    double vortClip = 4.0;
    std::vector<string> outputs;
    string framesDirVort  = "out/ns_filament_frames";
    string framesDirFlow  = "out/ns_filament_flow_frames";
    int    nParticles     = 1500;
    unsigned int particleSeed = 12345u;
    int    trailLen       = 24;
    int    trailStride    = 1;
    double bgDim          = 0.12;

    // ----------------------- filament parameters -----------------------------
    // Geometry & material (everything non-dimensional with rho_f=U=D=1):
    int    fil_N      = 80;        // segments
    double fil_L      = 2.0;       // length / D
    double fil_thick  = 0.05;      // h_s / D  (purely cosmetic for rendering)
    double fil_mass   = 0.1;       // m* = rho_s h / (rho_f D)  -> rhoLine = m*
    double K_B        = 0.01;      // EI / (rho_f Uinf^2 D^3)   -> EI  = K_B
    double K_S        = 1e4;       // EA / (rho_f Uinf^2 D)     -> EA  = K_S
    double fil_damp   = 0.05;      // structural Rayleigh-mass damping
    double fil_kick   = 0.01;      // initial transverse perturbation amplitude
    // Newton sub-stepping inside Cosserat::step (default 12, override here)
    // not exposed -- the rod stepper picks it.

    // IB coupling parameters.
    //   coupling = "oneway"  -- fluid drags rod, rod does NOT push back.  Drag
    //              follows the simple linear-Stokes form  F = -c_drag (V-u).
    //              Unconditionally stable, demo-quality "flag in the wake".
    //              c_drag has dimension of (force per unit length per velocity),
    //              non-dimensionalised by rho_f * U * D.
    //   coupling = "twoway"  -- partitioned direct-forcing IB with Newton's
    //              third law (Uhlmann 2005).  Physically self-consistent but
    //              suffers added-mass instability for light tails (m* < 1)
    //              without semi-implicit treatment; bring rod heavier or shrink
    //              ib_alpha in that regime.
    string coupling = "oneway";
    double c_drag   = 8.0;
    double ib_alpha = 0.5;
    int    ib_subiter = 0;

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("ns_filament_config.json")) cfgPath = "ns_filament_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        // flow
        ord = cfg.getInt("ord", ord);
        radius = cfg.getNumber("radius", radius);
        cx_ = cfg.getNumber("cx", cx_); cy_ = cfg.getNumber("cy", cy_);
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
        // filament
        fil_N      = cfg.getInt("fil_N", fil_N);
        fil_L      = cfg.getNumber("fil_L", fil_L);
        fil_thick  = cfg.getNumber("fil_thick", fil_thick);
        fil_mass   = cfg.getNumber("fil_mass", fil_mass);
        K_B        = cfg.getNumber("K_B", K_B);
        K_S        = cfg.getNumber("K_S", K_S);
        fil_damp   = cfg.getNumber("fil_damp", fil_damp);
        fil_kick   = cfg.getNumber("fil_kick", fil_kick);
        ib_alpha   = cfg.getNumber("ib_alpha", ib_alpha);
        ib_subiter = cfg.getInt("ib_subiter", ib_subiter);
        coupling   = cfg.getString("coupling", coupling);
        c_drag     = cfg.getNumber("c_drag", c_drag);
    }
    if (outputs.empty()) outputs = {"vorticity", "flow"};
    for (auto& s : outputs) {
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s != "vorticity" && s != "speed" && s != "flow") {
            cout << "ERROR: outputs entry '" << s << "' must be 'vorticity', 'speed', or 'flow'\n";
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

    cout << "DG NS + flexible-filament FSI (cylinder + Cosserat tail)\n";
    cout << "  config:  " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  domain:  [" << xa << "," << xb << "] x [" << ya << "," << yb << "]"
         << "  cylinder r=" << radius << " at (" << cx_ << "," << cy_ << ")\n";
    cout << "  Re=" << Re << "  Uinf=" << Uinf << "  D=" << D << "  nu=" << nu
         << "  dP" << ord << "  h=" << h << "\n";
    cout << "  dt=" << dt << "  t_end=" << t_end << "  steps=" << nsteps
         << "  save_every=" << save_every << "\n";
    cout << "  filament: L=" << fil_L << "  N=" << fil_N
         << "  m*=" << fil_mass << "  K_B=" << K_B << "  K_S=" << K_S
         << "  damp=" << fil_damp << "\n";
    cout << "  coupling=" << coupling
         << (coupling == "oneway" ? "  c_drag=" : "  ib_alpha=")
         << (coupling == "oneway" ? c_drag : ib_alpha)
         << (coupling == "twoway" ? "  ib_subiter=" : "")
         << (coupling == "twoway" ? std::to_string(ib_subiter) : string()) << "\n";

    // ----------------------- mesh + DG space -----------------------
    CylinderGeom geom{xa, xb, ya, yb, cx_, cy_, radius};
    Mesh mesh;
    generateCylinderMesh(mesh, geom, h, farRatio, grade);
    FEM fem(ord, mesh);
    MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    MatrixXi edge, edge2side; mesh.getEdge2Side(edge, edge2side);
    VectorXi tag = classifyEdges(mesh, edge, edge2side, geom);
    cout << "  nDof=" << nDof << " (x3 fields = " << 3 * nDof << ")\n";

    // ----------------------- BCs -----------------------
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
            default: break;
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
    double rcyl = radius, ccx = cx_, ccy = cy_, htol = 0.3 * h, Uin = Uinf, tr = rampTime;
    bc.velDir = [=](double x, double y, double t) -> Vector2d {
        if (std::hypot(x - ccx, y - ccy) < rcyl + htol) return Vector2d(0.0, 0.0);
        return Vector2d(Uin * ramp(t, tr), 0.0);
    };
    bc.velAcc = [=](double, double, double t) {
        return Vector2d(Uin * rampDer(t, tr), 0.0);
    };
    bc.presDir = [](double, double, double) { return 0.0; };

    NSIntegrator integ(fem, mesh, elem2dof, edge, edge2side, bc, nu, dt, sigma, 1.0,
                       gradDiv, pressureMode, ppeDivDamping);

    // ----------------------- IC (uniform flow, near-zero with ramp) -----
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
        const double xc = cx_ + 3.0 * D;
        const double sig2 = D * D;
        const double amp = 0.10 * perturb * Uinf * D;
        u0 = project([&](double x, double y){
            double X = x - xc, Y = y - cy_;
            double psi = amp * std::exp(-(X * X + 4.0 * Y * Y) / sig2);
            return psi * (-8.0 * Y / sig2);
        });
        v0 = project([&](double x, double y){
            double X = x - xc, Y = y - cy_;
            double psi = amp * std::exp(-(X * X + 4.0 * Y * Y) / sig2);
            return psi * (2.0 * X / sig2);
        });
    } else {
        u0 = project([&](double, double){ return Uinf; });
        v0 = project([&](double x, double y){
            double xi = x - cx_;
            if (xi <= 0) return 0.0;
            return perturb * Uinf * std::sin(M_PI * xi / (4.0 * D))
                   * std::exp(-(xi*xi + 4.0*(y-cy_)*(y-cy_)) / (9.0 * D * D));
        });
    }
    integ.setInitial(u0, v0);

    // ----------------------- filament setup -----------------------
    // Anchor at the cylinder's rear stagnation point on the body surface.
    // Embed the first segment slightly inside the cylinder so node 0 is in
    // the no-slip region (the IB load there is naturally absorbed); start
    // node 2 at  cx + radius  so the rod begins on the body surface.
    double xRoot = cx_ + radius;       // body surface at theta=0
    double yRoot = cy_;
    CosseratFilament rod;
    rod.initStraight(fil_N, xRoot, yRoot, xRoot + fil_L, yRoot, fil_mass, K_S, K_B);
    rod.dampStr = fil_damp;
    rod.clampRoot(true);                // clamp x,y of node 0 AND node 1
    // Initial transverse kick (a half-sine deflection) to seed the symmetry-
    // breaking that turns the wake on.  Small enough not to slap the no-slip
    // region.
    for (int i = 0; i <= fil_N; ++i) {
        double s = (double)i / fil_N;
        rod.X(i, 1) += fil_kick * std::sin(M_PI * s);
    }
    // Rest configuration is the original straight rod; we DON'T re-bake
    // X0 from the kicked shape -- the elastic restoring forces should pull
    // the rod back to straight.
    rod.V.setZero(); rod.A.setZero();

    MeshLocator locatorIB; locatorIB.build(mesh);
    std::vector<int> ibHint;

    // ----------------------- output channels (same machinery as ns_main) -
    int Hpix = std::max(2, (int)std::lround(Wpix * (ryb - rya) / (rxb - rxa)));
    if (Hpix % 2) ++Hpix; if (Wpix % 2) ++Wpix;
    auto inDomain = [=](double x, double y) {
        return std::hypot(x - cx_, y - cy_) > radius * 1.001 &&
               x > xa - 1e-9 && x < xb + 1e-9 && y > ya - 1e-9 && y < yb + 1e-9;
    };
    struct OutChannel {
        string kind;
        string dir;
        string mp4;
        std::unique_ptr<ParticleTracer> tracer;
    };
    std::vector<OutChannel> chans;
    chans.reserve(outputs.size());
    bool needLocator = false;
    int flowIdx = 0;
    for (const auto& kind : outputs) {
        OutChannel ch; ch.kind = kind;
        if (kind == "vorticity") { ch.dir = framesDirVort; ch.mp4 = "filament_vortex.mp4"; }
        else if (kind == "speed"){ ch.dir = framesDirVort; ch.mp4 = "filament_speed.mp4"; }
        else /* flow */          {
            ch.dir = framesDirFlow + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx)));
            ch.mp4 = "filament_flow"   + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx))) + ".mp4";
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
    MeshLocator locatorTracer;
    if (needLocator) locatorTracer.build(mesh);

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
            } else {
                writeParticlesPPM(fn, fem, mesh, elem2dof, integ.speed(), Wpix, Hpix,
                                  rxa, rxb, rya, ryb, 0.0, 1.6 * Uinf, *ch.tracer, inDomain, bgDim);
            }
            // Stamp the filament on top of every channel.
            int lineWidth = (ch.kind == "flow") ? 2 : 2;
            int lineRGB   = (ch.kind == "vorticity") ? 0x101010 : 0xFFFFFF;
            overlayPolylineOnPPM(string(fn), rod.X, rxa, rxb, rya, ryb, lineWidth, lineRGB);
        }
    };
    auto advanceAllTracers = [&](double dtStep) {
        for (auto& ch : chans)
            if (ch.tracer)
                ch.tracer->advance(fem, mesh, elem2dof, integ.u(), integ.v(), locatorTracer, dtStep);
    };

    // ----------------------- diagnostics file -----------------------
    ofstream hist("filament_diagnostics.csv");
    hist << "t,CD,CL,xtip,ytip,maxStrain\n";
    std::vector<double> tCL, CLval, CDval, yTipHist;

    // ----------------------- FSI time loop -----------------------
    cout << "\nEvolving with IB direct-forcing FSI...\n";
    auto t0wall = chrono::high_resolution_clock::now();
    int frame = 0;
    writeAllFrames(frame);
    ++frame;
    const double qdyn = 0.5 * Uinf * Uinf * D;

    VectorXd loadFx = VectorXd::Zero(nDof);
    VectorXd loadFy = VectorXd::Zero(nDof);

    // Contact flags from the previous step's no-penetration projection.  Used
    // by the next step's drag transfer to skip rod nodes pinned on the
    // cylinder surface (otherwise the standing c_drag*(V-u) signal at those
    // nodes drifts the fluid solver after sustained grazing contact).
    std::vector<char> contactPrev(rod.N + 1, 0);

    for (int s = 1; s <= nsteps; ++s) {
        // ----------------------- partitioned FSI step -----------------------
        // (1) Sample the t^n fluid velocity at every rod node.
        auto sampleN = meshToRod(fem, mesh, elem2dof, integ.u(), integ.v(), locatorIB, rod, ibHint);
        VectorXd ds = rodArclengthWeights(rod);
        rod.Fext.setZero();
        loadFx.setZero(); loadFy.setZero();

        if (coupling == "twoway") {
            // Direct-forcing IB (Uhlmann 2005, Bhalla 2013) -- partitioned.
            //   F^IB_k = alpha * rho_f * (V_rod_k - u_h(X_k)),  alpha = ib_alpha/dt
            //   force on fluid =  +F^IB                        (pulls fluid -> rod)
            //   force on rod    =  -F^IB * ds                  (Newton's 3rd law)
            // Stable only for heavy rods (m* >= ~1) without sub-iterations.
            const double alpha = ib_alpha / dt;
            MatrixXd Fk(rod.N + 1, 2); Fk.setZero();
            // explicit Fext on the rod for two-way -- twoway is the experimental
            // option, kept simple
            rod.dragCoef.resize(0); rod.dragRef.resize(0, 0);
            for (int k = 0; k <= rod.N; ++k) {
                if (k <= 1) continue;
                if (!sampleN.alive[k]) continue;
                Fk(k, 0) = alpha * (rod.V(k, 0) - sampleN.uv(k, 0));
                Fk(k, 1) = alpha * (rod.V(k, 1) - sampleN.uv(k, 1));
                rod.Fext(2 * k    ) = -Fk(k, 0) * ds(k);
                rod.Fext(2 * k + 1) = -Fk(k, 1) * ds(k);
            }
            rodToMesh(fem, mesh, elem2dof, locatorIB, rod, Fk, ds, ibHint, loadFx, loadFy);
        } else {
            // One-way drag: fluid pushes rod, rod doesn't react on fluid.
            // Fed into the rod stepper as IMPLICIT linear drag so it can be
            // stiff without breaking Newton convergence:
            //   F_drag_k = -c_drag * (V_k - u_h(X_k)) * ds_k
            // Always stable; the rod becomes a passive flag in the wake.
            // Contact nodes from the previous step are excluded so they don't
            // re-inject persistent rod-fluid mismatch into the system.
            rod.dragCoef.setZero(rod.N + 1);
            rod.dragRef.setZero(rod.N + 1, 2);
            for (int k = 0; k <= rod.N; ++k) {
                if (k <= 1) continue;
                if (!sampleN.alive[k]) continue;
                if (k < (int)contactPrev.size() && contactPrev[k]) continue;
                rod.dragCoef(k) = c_drag * ds(k);
                rod.dragRef(k, 0) = sampleN.uv(k, 0);
                rod.dragRef(k, 1) = sampleN.uv(k, 1);
            }
            // loadFx/Fy stay zero -> the fluid is solved exactly as the bare
            // cylinder would be.
        }

        // (2) Advance the rod (implicit Newmark).
        if (!rod.step(dt)) { cout << "ERROR: rod solve failed at step " << s << "\n"; return 1; }

        // (2b) Cylinder no-penetration projection: any rod node that landed
        // inside (or within an epsilon shell of) the cylinder is pushed back
        // onto its surface, and its inward radial velocity is REFLECTED
        // (elastic bounce, restitution coefficient e).  Tangential friction
        // bleeds the grazing-contact velocity that otherwise pumps numerical
        // noise across the boundary.
        //
        // Critically, contact nodes are also marked so the *next* step's IB
        // drag transfer skips them:  a node pinned on (or inside) the
        // cylinder no-slip Dirichlet halo doesn't represent a free fluid
        // boundary anymore, so feeding c_drag * (V_node - u_h) back into the
        // rod -> fluid pipeline at those nodes produces a self-feeding
        // instability that blows the run up after sustained grazing contact.
        //
        // The contact zone is set generously to  radius + 1.5*htol  so it
        // covers the same neighbourhood that NSIntegrator's velDir() forces
        // to zero, plus a small buffer so a node oscillating right on the
        // halo boundary doesn't keep flipping in/out of contact every step.
        const double r_contact = radius + 1.5 * htol;
        const double r_safe    = radius * 1.005;
        const double r2_contact = r_contact * r_contact;
        const double r2_safe    = r_safe * r_safe;
        const double restitution = 0.3;
        const double tang_friction = 0.5;
        std::vector<char> inContact(rod.N + 1, 0);
        for (int k = 2; k <= rod.N; ++k) {
            double dx = rod.X(k, 0) - cx_;
            double dy = rod.X(k, 1) - cy_;
            double r2k = dx * dx + dy * dy;
            if (r2k >= r2_contact) continue;
            inContact[k] = 1;
            double rk = std::sqrt(std::max(r2k, 1e-300));
            double nx = (rk > 0) ? (dx / rk) : 1.0;
            double ny = (rk > 0) ? (dy / rk) : 0.0;
            // Geometric no-penetration only fires inside r_safe.
            if (r2k < r2_safe) {
                rod.X(k, 0) = cx_ + r_safe * nx;
                rod.X(k, 1) = cy_ + r_safe * ny;
                double vn = rod.V(k, 0) * nx + rod.V(k, 1) * ny;
                if (vn < 0.0) {
                    double s = (1.0 + restitution) * vn;
                    rod.V(k, 0) -= s * nx;
                    rod.V(k, 1) -= s * ny;
                }
                double vx = rod.V(k, 0), vy = rod.V(k, 1);
                double vt_x = vx - (vx * nx + vy * ny) * nx;
                double vt_y = vy - (vx * nx + vy * ny) * ny;
                rod.V(k, 0) -= tang_friction * vt_x;
                rod.V(k, 1) -= tang_friction * vt_y;
                double an = rod.A(k, 0) * nx + rod.A(k, 1) * ny;
                if (an < 0.0) {
                    rod.A(k, 0) -= an * nx;
                    rod.A(k, 1) -= an * ny;
                }
            }
            // Even if r > r_safe (no projection needed), still flag it as
            // "in contact" so the next step's drag transfer skips this node.
        }
        // Stash for next step's drag setup.
        contactPrev = inContact;
        // (3) Advance the fluid.
        if (!integ.stepWithBodyForce(s * dt, loadFx, loadFy)) {
            cout << "ERROR: NS solve failed at step " << s << "\n"; return 1;
        }
        if (needLocator)
            advanceAllTracers(dt);

        // (5) Diagnostics.
        double Fx, Fy; integ.cylinderForce(isCyl, Fx, Fy);
        double CD = Fx / qdyn, CL = Fy / qdyn;
        double t = s * dt;
        double xtip = rod.X(rod.N, 0);
        double ytip = rod.X(rod.N, 1);
        double strain = rod.maxEdgeStrain();
        hist << t << "," << CD << "," << CL << "," << xtip << "," << ytip << "," << strain << "\n";
        tCL.push_back(t); CLval.push_back(CL); CDval.push_back(CD); yTipHist.push_back(ytip - yRoot);

        if (s % save_every == 0 || s == nsteps) {
            writeAllFrames(frame);
            ++frame;
            double el = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
            cout << "  step " << setw(6) << s << "/" << nsteps << "  t=" << fixed << setprecision(2) << t
                 << "  CD=" << setprecision(3) << CD << "  CL=" << setw(7) << CL
                 << "  ytip=" << setw(7) << (ytip - yRoot) << "  strain=" << setprecision(2) << scientific << strain
                 << fixed << "  frame=" << setw(4) << (frame-1) << "  (" << setprecision(1) << el << "s)\n";
        }
    }
    hist.close();

    // ----------------------- Strouhal / tip stats -----------------------
    auto crossPeriod = [&](const std::vector<double>& y) {
        std::vector<double> cs;
        for (size_t i = 1; i < y.size(); ++i)
            if (tCL[i] > 0.5 * t_end && y[i-1] < 0 && y[i] >= 0) {
                double frac = -y[i-1] / (y[i] - y[i-1]);
                cs.push_back(tCL[i-1] + frac * (tCL[i] - tCL[i-1]));
            }
        if (cs.size() < 2) return 0.0;
        return (cs.back() - cs.front()) / (cs.size() - 1);
    };
    double T_CL  = crossPeriod(CLval);
    double T_tip = crossPeriod(yTipHist);
    double sumCD = 0; double maxYt = -1e30, minYt = 1e30; int cnt = 0;
    for (size_t i = 0; i < CLval.size(); ++i) if (tCL[i] > 0.5 * t_end) {
        sumCD += CDval[i];
        maxYt = std::max(maxYt, yTipHist[i]);
        minYt = std::min(minYt, yTipHist[i]);
        ++cnt;
    }
    double meanCD = cnt ? sumCD / cnt : 0.0;
    double tipAmp = 0.5 * (maxYt - minYt);

    double wall = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
    cout << "\nDone. " << frame << " frames in";
    for (size_t i = 0; i < chans.size(); ++i)
        cout << (i ? ", " : " ") << "'" << chans[i].dir << "/'";
    cout << ".  wall=" << fixed << setprecision(1) << wall << "s\n";
    cout << "  mean CD = " << setprecision(3) << meanCD << "\n";
    if (T_CL > 0)
        cout << "  Strouhal (lift)         St_CL  = " << setprecision(4) << D / (Uinf * T_CL)
             << "  (period " << setprecision(3) << T_CL << ")\n";
    if (T_tip > 0)
        cout << "  Strouhal (tip flap)     St_tip = " << setprecision(4) << D / (Uinf * T_tip)
             << "  (period " << setprecision(3) << T_tip << ")\n";
    cout << "  tip transverse amplitude  A/D = " << setprecision(3) << tipAmp << "\n";
    cout << "  diagnostics -> filament_diagnostics.csv\n";

    cout << "\nAssemble the movies:\n";
    for (const auto& ch : chans)
        cout << "  ffmpeg -y -framerate 25 -i " << ch.dir << "/frame_%05d.ppm"
             << " -c:v libx264 -pix_fmt yuv420p -crf 18 " << ch.mp4 << "\n";
    return 0;
}
