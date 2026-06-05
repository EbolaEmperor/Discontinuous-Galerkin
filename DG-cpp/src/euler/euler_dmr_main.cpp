// ===========================================================================
// Double Mach Reflection (Woodward & Colella, 1984) with the DG + 2nd-order IMEX
// Euler solver.  A Mach-10 shock in air (gamma=1.4), inclined 60 deg to the
// x-axis, strikes a reflecting wall at x=1/6 and produces the classic
// self-similar pattern: incident / reflected / Mach-stem shocks, two triple
// points, and a slip line whose Kelvin-Helmholtz roll-up forms the near-wall
// jet "vortex street".  High-order DG + sub-cell artificial viscosity (treated
// IMPLICITLY in the IMEX scheme) + a positivity-preserving limiter capture the
// strong shocks while keeping the slip-line instability sharp.
//
// Output: a density movie (PPM frames -> ffmpeg) over the lower domain, plus a
// final full-resolution density image, a numerical-Schlieren image, and a
// ZOOMED density+Schlieren inset on the triple-point / slip-line region that
// reveals the vortex roll-ups.  All parameters live in a JSON config.
//
//   Domain [0,4]x[0,1];  T=0.2.
//   Pre-shock  (rho,u,v,p) = (1.4, 0, 0, 1).
//   Post-shock (rho,u,v,p) = (8, 8.25 cos30, -8.25 sin30, 116.5) = (8, 7.1447, -4.125, 116.5).
//   Initial shock line x = 1/6 + y/sqrt(3); top BC shock trace x_s(t) = 1/6 + (1+20t)/sqrt(3).
// ===========================================================================
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "EulerDG.h"
#include "FEM.h"
#include "Mesh.h"
#include "Json.h"

using namespace Eigen;
using namespace std;
using namespace euler;
namespace fs = std::filesystem;

// edge tags
enum { INTERIOR = 0, LEFT = 1, RIGHT = 2, BOTTOM = 3, TOP = 4 };

int main(int argc, char** argv) {
    std::cout << std::unitbuf;   // flush each line so the run is monitorable via a redirected log
    // ----------------------- parameters (JSON-overridable) -----------------------
    int    ord       = 2;          // dP_k
    int    ny        = 100;        // cells in y (nx = round(domain_xb*ny)); h = 1/ny
    double domain_xb = 4.0;        // domain right edge (extend >4 to run past t=0.2 in-frame)
    double t_end     = 0.2;
    double cfl       = 0.15;
    double lambda_safe = 23.0;     // safe upper bound on |u|+c for the fixed dt
    int    n_frames  = 200;
    bool   full_movie = true;      // render the full-domain density frame sequence
    bool   zoom_movie = true;      // render the zoom density + zoom schlieren frame sequences
    bool   save_state = false;     // dump the final state to dmr_state.bin (re-render zooms later)
    string state_only = "";        // if set: load this state, render stills only, exit (no time stepping)
    // artificial viscosity / limiter (defaults tuned for a stable Mach-10 DMR:
    // pressure sensor keeps the slip line sharp; the per-side eps-weighted SIPG
    // form is genuinely SPSD so a standard penalty sigma_ip~20 suffices)
    double av_c      = 2.0;
    double av_kappa  = 1.0;
    double sigma_ip  = 20.0;
    int    av_refresh = 5;
    bool   use_av    = true;
    bool   use_pos   = true;
    // rendering
    int    Wpix      = 1600;       // movie width
    double rxa = 0.0, rxb = 3.2, rya = 0.0, ryb = 1.0;     // movie window
    double rho_lo = 1.4, rho_hi = 21.0;                    // density colour range
    // zoom window (triple point / slip-line vortex street at t=0.2)
    double zxa = 2.0, zxb = 2.9, zya = 0.0, zyb = 0.55;
    string framesDir = "dmr_frames";
    string cmapName  = "inferno";

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("dmr_config.json")) cfgPath = "dmr_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        ord = cfg.getInt("ord", ord);
        ny = cfg.getInt("ny", ny);
        domain_xb = cfg.getNumber("domain_xb", domain_xb);
        t_end = cfg.getNumber("t_end", t_end);
        full_movie = cfg.getBool("full_movie", full_movie);
        zoom_movie = cfg.getBool("zoom_movie", zoom_movie);
        save_state = cfg.getBool("save_state", save_state);
        state_only = cfg.getString("state_only", state_only);
        cfl = cfg.getNumber("cfl", cfl);
        lambda_safe = cfg.getNumber("lambda_safe", lambda_safe);
        n_frames = cfg.getInt("n_frames", n_frames);
        av_c = cfg.getNumber("av_c", av_c);
        av_kappa = cfg.getNumber("av_kappa", av_kappa);
        sigma_ip = cfg.getNumber("sigma_ip", sigma_ip);
        av_refresh = cfg.getInt("av_refresh", av_refresh);
        use_av = cfg.getBool("use_av", use_av);
        use_pos = cfg.getBool("use_positivity", use_pos);
        Wpix = cfg.getInt("Wpix", Wpix);
        rxa = cfg.getNumber("render_xa", rxa); rxb = cfg.getNumber("render_xb", rxb);
        rya = cfg.getNumber("render_ya", rya); ryb = cfg.getNumber("render_yb", ryb);
        rho_lo = cfg.getNumber("rho_lo", rho_lo); rho_hi = cfg.getNumber("rho_hi", rho_hi);
        zxa = cfg.getNumber("zoom_xa", zxa); zxb = cfg.getNumber("zoom_xb", zxb);
        zya = cfg.getNumber("zoom_ya", zya); zyb = cfg.getNumber("zoom_yb", zyb);
        framesDir = cfg.getString("frames_dir", framesDir);
        cmapName = cfg.getString("cmap", cmapName);
    }
    int nx = (int)std::lround(domain_xb * ny);
    double h = 1.0 / ny;
    Colormap cmap = (cmapName == "jet") ? CM_JET : (cmapName == "viridis") ? CM_VIRIDIS
                  : (cmapName == "gray") ? CM_GRAY : CM_INFERNO;

    // ----------------------- states (verified Woodward-Colella) -----------------------
    const double SQ3 = std::sqrt(3.0);
    auto preP  = []{ return Vector4d(1.4, 0.0, 0.0, 1.0); };
    auto postP = []{ return Vector4d(8.0, 8.25 * std::cos(M_PI / 6.0), -8.25 * std::sin(M_PI / 6.0), 116.5); };
    Vector4d POST = primToCons(postP()(0), postP()(1), postP()(2), postP()(3));
    Vector4d PRE  = primToCons(preP()(0),  preP()(1),  preP()(2),  preP()(3));
    auto initPrim = [&](double x, double y) {
        return (x < 1.0/6.0 + y / SQ3) ? postP() : preP();   // post behind the inclined shock
    };

    cout << "DG + 2nd-order IMEX Euler -- Double Mach Reflection (Mach 10)\n";
    cout << "  config: " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  domain [0," << domain_xb << "]x[0,1]  dP" << ord << "  nx=" << nx << " ny=" << ny << " h=" << h
         << "  (" << 2*nx*ny << " tris)\n";
    cout << "  AV=" << (use_av?"on":"off") << " (c=" << av_c << ", kappa=" << av_kappa << ", sigma=" << sigma_ip
         << ", refresh=" << av_refresh << ")  positivity=" << (use_pos?"on":"off") << "  t_end=" << t_end << "\n";

    // ----------------------- mesh + DG space -----------------------
    Mesh mesh; makeRectMesh(mesh, 0, domain_xb, 0, 1, nx, ny);
    FEM fem(ord, mesh);
    MatrixXi e2d; int nDof; fem.getDOF(mesh, e2d, nDof);
    MatrixXi edge, e2s; mesh.getEdge2Side(edge, e2s);
    int NE = edge.rows();
    VectorXi tag = VectorXi::Zero(NE);
    for (int e = 0; e < NE; ++e) {
        if (e2s(e, 0) != -1 && e2s(e, 1) != -1) continue;       // interior
        Vector2d m = 0.5 * (mesh.node.row(edge(e,0)) + mesh.node.row(edge(e,1))).transpose();
        if (std::abs(m.x() - 0.0) < 1e-9) tag(e) = LEFT;
        else if (std::abs(m.x() - domain_xb) < 1e-9) tag(e) = RIGHT;
        else if (std::abs(m.y() - 0.0) < 1e-9) tag(e) = BOTTOM;
        else if (std::abs(m.y() - 1.0) < 1e-9) tag(e) = TOP;
    }
    cout << "  nDof=" << nDof << " (x4 fields = " << 4*nDof << ")\n";

    EulerConfig cfg; cfg.use_av = use_av; cfg.use_positivity = use_pos;
    cfg.av_c = av_c; cfg.av_kappa = av_kappa; cfg.sigma_ip = sigma_ip; cfg.av_refresh = av_refresh;
    EulerDG dg(fem, mesh, e2d, edge, e2s, tag, cfg);
    dg.setState(projectInitial(fem, mesh, e2d, initPrim));

    // ----------------------- boundary ghost states -----------------------
    ExteriorStateFn bc = [&](double x, double y, double t, const Vector4d& Uin, double nx_, double ny_, int tg) -> Vector4d {
        switch (tg) {
            case LEFT:  return POST;                                    // post-shock inflow
            case RIGHT: return Uin;                                     // supersonic outflow (copy)
            case BOTTOM:
                if (x < 1.0/6.0) return POST;                          // post-shock inflow strip
                else {                                                  // reflecting/slip wall: mirror normal momentum
                    double mn = Uin(1)*nx_ + Uin(2)*ny_;
                    return Vector4d(Uin(0), Uin(1) - 2*mn*nx_, Uin(2) - 2*mn*ny_, Uin(3));
                }
            case TOP: {
                double xs = 1.0/6.0 + (1.0 + 20.0 * t) / SQ3;          // moving incident-shock trace
                return (x < xs) ? POST : PRE;
            }
            default: return Uin;
        }
    };

    // ----------------------- fixed dt -----------------------
    double dt = cfl * h / ((2.0 * ord + 1.0) * lambda_safe);
    int nsteps = std::max(1, (int)std::lround(t_end / dt));
    dt = t_end / nsteps;
    int save_every = std::max(1, nsteps / std::max(1, n_frames));
    cout << "  dt=" << scientific << setprecision(3) << dt << "  steps=" << nsteps
         << "  save_every=" << save_every << "  (lambda_safe=" << fixed << setprecision(1) << lambda_safe << ")\n";

    // ----------------------- rendering helpers -----------------------
    if (Wpix % 2) ++Wpix;
    auto writeDensity = [&](const string& path, double xa, double xb, double ya, double yb, int W) {
        int Hh = std::max(2, (int)std::lround(W * (yb - ya) / (xb - xa))); if (Hh%2) ++Hh;
        writeScalarPPM(path, fem, mesh, e2d, dg.densityField(), W, Hh, xa, xb, ya, yb, rho_lo, rho_hi, cmap, {});
    };
    auto writeSchlieren = [&](const string& path, double xa, double xb, double ya, double yb, int W) {
        int Hh = std::max(2, (int)std::lround(W * (yb - ya) / (xb - xa))); if (Hh%2) ++Hh;
        writeSchlierenPPM(path, fem, mesh, e2d, dg.densityField(), W, Hh, xa, xb, ya, yb, 8.0);
    };
    auto renderStills = [&]() {        // final-instant stills: full + zoom, density + schlieren
        writeDensity ("dmr_density.ppm",        rxa, rxb, rya, ryb, Wpix);
        writeDensity ("dmr_density_zoom.ppm",   zxa, zxb, zya, zyb, Wpix);
        writeSchlieren("dmr_schlieren.ppm",      rxa, rxb, rya, ryb, Wpix);
        writeSchlieren("dmr_schlieren_zoom.ppm", zxa, zxb, zya, zyb, Wpix);
    };

    // ----------------------- state-only mode: re-render stills from a saved state -----------------------
    if (!state_only.empty()) {
        std::ifstream in(state_only, std::ios::binary);
        if (!in) { cout << "ERROR: cannot open " << state_only << "\n"; return 1; }
        int r = 0, c = 0; in.read((char*)&r, 4); in.read((char*)&c, 4);
        MatrixXd U(r, c); in.read((char*)U.data(), (std::streamsize)r * c * sizeof(double));
        dg.setState(U);
        cout << "  [state_only] loaded " << r << "x" << c << " from " << state_only
             << "; rendering stills with zoom [" << zxa << "," << zxb << "]x[" << zya << "," << zyb << "]\n";
        renderStills();
        cout << "  stills: dmr_density.ppm dmr_density_zoom.ppm dmr_schlieren.ppm dmr_schlieren_zoom.ppm\n";
        return 0;
    }

    // ----------------------- frame-sequence output setup -----------------------
    string zoomDir = framesDir + "_zoom", schZoomDir = framesDir + "_schzoom";
    std::vector<string> dirs;
    if (full_movie) dirs.push_back(framesDir);
    if (zoom_movie) { dirs.push_back(zoomDir); dirs.push_back(schZoomDir); }
    for (const string& d : dirs) {
        fs::create_directories(d);
        for (const auto& e : fs::directory_iterator(d))
            if (e.path().extension() == ".ppm") fs::remove(e.path());
    }
    auto writeFrameSet = [&](int idx) {
        char fn[600];
        if (full_movie) { snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", framesDir.c_str(), idx); writeDensity(fn, rxa, rxb, rya, ryb, Wpix); }
        if (zoom_movie) {
            snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", zoomDir.c_str(), idx);    writeDensity(fn, zxa, zxb, zya, zyb, Wpix);
            snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", schZoomDir.c_str(), idx); writeSchlieren(fn, zxa, zxb, zya, zyb, Wpix);
        }
    };

    // ----------------------- time loop -----------------------
    cout << "\nEvolving (ARS(2,2,2) IMEX: explicit HLLC flux + implicit artificial viscosity)...\n";
    auto t0 = chrono::high_resolution_clock::now();
    int frame = 0;
    writeFrameSet(frame++);
    for (int s = 1; s <= nsteps; ++s) {
        if (!dg.step(dt, s * dt, bc)) { cout << "ERROR: step " << s << " failed\n"; return 1; }
        if (s % save_every == 0 || s == nsteps) {
            writeFrameSet(frame++);
            double lam = dg.maxWaveSpeedGlobal();
            Vector4d tot = dg.conservedTotals();
            double el = chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count();
            VectorXd rho = dg.densityField();
            cout << "  step " << setw(6) << s << "/" << nsteps << "  t=" << fixed << setprecision(4) << s*dt
                 << "  rho in [" << setprecision(3) << rho.minCoeff() << "," << rho.maxCoeff() << "]"
                 << "  CFL=" << setprecision(2) << dt*lam*(2*ord+1)/h
                 << "  mass=" << setprecision(4) << tot(0)
                 << "  (" << setprecision(1) << el << "s)\n";
        }
    }

    // ----------------------- final stills + optional state dump -----------------------
    renderStills();
    if (save_state) {
        std::ofstream out("dmr_state.bin", std::ios::binary);
        const MatrixXd& U = dg.state(); int r = (int)U.rows(), c = (int)U.cols();
        out.write((char*)&r, 4); out.write((char*)&c, 4);
        out.write((const char*)U.data(), (std::streamsize)r * c * sizeof(double));
        cout << "  final state saved to dmr_state.bin (" << r << "x" << c << ")\n";
    }
    double wall = chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count();
    cout << "\nDone. " << frame << " frames.  wall=" << fixed << setprecision(1) << wall << "s\n";
    cout << "  final stills: dmr_density.ppm, dmr_density_zoom.ppm, dmr_schlieren.ppm, dmr_schlieren_zoom.ppm\n";
    if (full_movie || zoom_movie) cout << "\nAssemble the movie(s):\n";
    if (full_movie) cout << "  ffmpeg -y -framerate 25 -i " << framesDir  << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr.mp4\n";
    if (zoom_movie) cout << "  ffmpeg -y -framerate 25 -i " << zoomDir    << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_zoom.mp4\n"
                         << "  ffmpeg -y -framerate 25 -i " << schZoomDir << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_schlieren_zoom.mp4\n";
    return 0;
}
