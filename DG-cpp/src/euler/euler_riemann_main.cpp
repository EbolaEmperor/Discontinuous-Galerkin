// ===========================================================================
// Two-dimensional Riemann problems (Schulz-Rinne / Lax-Liu configurations) with
// the DG + 2nd-order IMEX Euler solver on an adaptive conforming mesh.
//
// Four constant states fill the four quadrants of the unit square; each of the
// four interfaces is a single 1-D elementary wave (shock, rarefaction, or
// contact).  Their nonlinear interaction at the centre spawns a zoo of genuinely
// 2-D structures -- Mach stems, curved shocks, spiral contact roll-up and, in the
// four-shock configs, a Kelvin-Helmholtz "mushroom jet".  The IC is piecewise
// constant and the BCs are transmissive, so it is a clean, cheap stress test of
// contact resolution (HLLC), positivity (Zhang-Shu) and symmetry.
//
// Selectable configurations (config = ...; primitive states (rho,u,v,p)):
//   3  (four shocks; t_end 0.3): the iconic central mushroom jet with KH roll-up.
//   6  (four contacts/shears; t_end 0.3): a pinwheel of four rolling shear layers.
//   12 (two shocks + two contacts; t_end 0.25): two logarithmic-spiral contacts.
//
//   Domain [0,1]x[0,1]; gamma=1.4.  The four-state intersection is at `cross` (default
//   0.8, the Schulz-Rinne OFF-CENTER convention) so the structure grows into the large
//   lower-left sub-square and a fixed [0,1]^2 stays clean to a LONGER time (config 3 to
//   t=0.8, as in HOCUS-BVD / PyClaw).  BC is a WEAK far-field ghost state (always the
//   constant quadrant state) through the HLLC flux, handled per-characteristic -- this
//   removes both the supersonic-inflow-wall artefact and the subsonic sheared-outflow
//   instability (the right-wall blob) that a plain zero-gradient BC grows at long time.
//   NOTE: the criss-cross base mesh is not symmetric under the diagonal
//   reflections of the IC, so the roll-ups are a qualitative showpiece, not a
//   symmetry benchmark.
// ===========================================================================
#include <cmath>
#include <iostream>
#include <string>

#include "EulerDG.h"
#include "Json.h"
#include "Mesh.h"
#include "euler_amr_scene.h"

using namespace Eigen;
using namespace std;
using namespace euler;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    std::cout << std::unitbuf;

    int config = 3;
    int base_ny = 40;
    // Quadrant-intersection point.  cross=0.5 is the centered-cross benchmark (run only
    // to the standard short time, e.g. config-3 t=0.3, before the fan reaches a wall).
    // cross=0.8 is the Schulz-Rinne OFF-CENTER convention: the structure develops into
    // the large lower-left sub-square, so a fixed [0,1]^2 stays clean to a LONGER time
    // (HOCUS-BVD / PyClaw run config 3 to t=0.8 this way) -- the default here.
    double cross = 0.8;

    AMRScene S;
    S.name = "riemann"; S.framesDir = "riemann_frames"; S.title = "2D Riemann problem";
    S.ord = 2; S.max_gen = 4; S.cfl = 0.2; S.lambda_safe = 5.0; S.n_frames = 200;
    S.remesh_every = 5; S.buffer_layers = 3; S.init_passes = 10; S.th_ref = 0.25; S.th_crs = 0.06;
    // The 2D Riemann configs are CONTACT/SHEAR dominated (e.g. config 6 is four slip
    // lines, config 3 has a low-density shear-jet).  Zhang-Shu positivity preserves
    // the cell mean but NOT positive high-order overshoots, and there is no slope
    // limiter -- so the under-resolved Kelvin-Helmholtz convergence grows unbounded
    // (rho -> 1e5-1e6) unless artificial viscosity damps the contacts.  Recipe that
    // keeps all configs stable AND physical: DENSITY sensor (fires on contacts) +
    // a LOW Persson threshold av_s0=-3 (AV engages early, before the overshoot runs
    // away) + av_c=2 + cfl=0.2.  (Rusanov flux alone does NOT fix it; the early AV
    // threshold is the decisive lever.)
    S.av_c = 2.0; S.av_kappa = 1.0; S.sigma_ip = 20.0; S.av_refresh = 5; S.av_indicator = 0;
    S.av_s0 = -3.0; S.av_s0_set = true;
    S.use_av = true; S.use_pos = true;
    S.Wpix = 1400; S.cmapName = "inferno";
    S.rxa = 0.0; S.rxb = 1.0; S.rya = 0.0; S.ryb = 1.0;

    // read just `config` + `base_ny` first so config-specific defaults can be set,
    // then let the file override anything else.
    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("riemann_config.json")) cfgPath = "riemann_config.json";
    cfgjson::Json cfg; bool haveCfg = false;
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        try { cfg = cfgjson::parse(text); haveCfg = true; }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        config = cfg.getInt("config", config);
    }

    // primitive (rho,u,v,p) of the four quadrants about (0.5,0.5):
    //   NE = x>=0.5,y>=0.5 ; NW = x<0.5,y>=0.5 ; SW = x<0.5,y<0.5 ; SE = x>=0.5,y<0.5
    Vector4d NE, NW, SW, SE;
    if (config == 6) {
        NE = Vector4d(1.0,  0.75, -0.5, 1.0);
        NW = Vector4d(2.0,  0.75,  0.5, 1.0);
        SW = Vector4d(1.0, -0.75,  0.5, 1.0);
        SE = Vector4d(3.0, -0.75, -0.5, 1.0);
        S.t_end = 0.3; S.rho_lo = 0.5; S.rho_hi = 3.0;
    } else if (config == 12) {
        NE = Vector4d(0.5313, 0.0,    0.0,    0.4);
        NW = Vector4d(1.0,    0.7276, 0.0,    1.0);
        SW = Vector4d(0.8,    0.0,    0.0,    1.0);
        SE = Vector4d(1.0,    0.0,    0.7276, 1.0);
        S.t_end = 0.25; S.rho_lo = 0.4; S.rho_hi = 1.8;
    } else { // config 3 (four shocks) -- the famous mushroom jet
        config = 3;
        NE = Vector4d(1.5,    0.0,   0.0,   1.5);
        NW = Vector4d(0.5323, 1.206, 0.0,   0.3);
        SW = Vector4d(0.138,  1.206, 1.206, 0.029);
        SE = Vector4d(0.5323, 0.0,   1.206, 0.3);
        S.t_end = 0.3; S.rho_lo = 0.1; S.rho_hi = 1.8;
    }

    if (haveCfg) {
        S.ord = cfg.getInt("ord", S.ord);
        base_ny = cfg.getInt("base_ny", base_ny);
        cross = cfg.getNumber("cross", cross);
        S.cfl_rho_floor = cfg.getNumber("cfl_rho_floor", S.cfl_rho_floor);
        S.max_gen = cfg.getInt("max_gen", S.max_gen);
        S.t_end = cfg.getNumber("t_end", S.t_end);
        S.cfl = cfg.getNumber("cfl", S.cfl);
        S.lambda_safe = cfg.getNumber("lambda_safe", S.lambda_safe);
        S.n_frames = cfg.getInt("n_frames", S.n_frames);
        S.remesh_every = cfg.getInt("remesh_every", S.remesh_every);
        S.buffer_layers = cfg.getInt("buffer_layers", S.buffer_layers);
        S.th_ref = cfg.getNumber("th_ref", S.th_ref);
        S.th_crs = cfg.getNumber("th_crs", S.th_crs);
        S.init_passes = cfg.getInt("init_passes", S.init_passes);
        S.av_c = cfg.getNumber("av_c", S.av_c);
        S.sigma_ip = cfg.getNumber("sigma_ip", S.sigma_ip);
        S.av_refresh = cfg.getInt("av_refresh", S.av_refresh);
        S.av_indicator = cfg.getInt("av_indicator", S.av_indicator);
        S.use_av = cfg.getBool("use_av", S.use_av);
        S.use_pos = cfg.getBool("use_positivity", S.use_pos);
        S.use_hllc = cfg.getBool("use_hllc", S.use_hllc);
        if (cfg.hasNumber("av_s0")) { S.av_s0 = cfg.getNumber("av_s0", S.av_s0); S.av_s0_set = true; }
        S.use_tvb = cfg.getBool("use_tvb", S.use_tvb);
        S.tvb_M = cfg.getNumber("tvb_M", S.tvb_M);
        S.Wpix = cfg.getInt("Wpix", S.Wpix);
        S.rho_lo = cfg.getNumber("rho_lo", S.rho_lo);
        S.rho_hi = cfg.getNumber("rho_hi", S.rho_hi);
        S.cmapName = cfg.getString("cmap", S.cmapName);
        S.full_movie = cfg.getBool("full_movie", S.full_movie);
        S.mesh_movie = cfg.getBool("mesh_movie", S.mesh_movie);
        S.framesDir = cfg.getString("frames_dir", S.framesDir);
    }
    S.name = "riemann_cfg" + std::to_string(config);
    S.framesDir = "riemann_cfg" + std::to_string(config) + "_frames";

    cout << "2D Riemann problem, configuration " << config << "  (t_end=" << S.t_end << ")\n";

    auto quad = [=](double x, double y) -> Vector4d {
        return (x >= cross) ? (y >= cross ? NE : SE) : (y >= cross ? NW : SW);
    };
    S.primIC = quad;
    S.tagEdge = [=](double mx, double my) -> int {
        if (std::abs(mx) < 1e-9 || std::abs(mx - 1.0) < 1e-9 ||
            std::abs(my) < 1e-9 || std::abs(my - 1.0) < 1e-9) return 1;
        return 0;
    };
    // WEAK far-field ghost-state BC (Trixi.jl's BoundaryConditionDirichlet(initial_condition)
    // idiom): the ghost is ALWAYS the constant quadrant far-field state, fed through the
    // interior HLLC flux.  The Riemann solver then handles every wall PER CHARACTERISTIC --
    // supersonic outflow ignores the ghost, subsonic blends in the far-field, supersonic
    // inflow imposes it fully.  This both fixes the supersonic-INFLOW walls (left/bottom)
    // and pins the inward acoustic on the subsonic sheared OUTFLOW walls, removing the
    // spurious right-wall blob that a plain zero-gradient ("return Uin") outflow grows
    // under the SE tangential shear at longer times.  Combined with the off-center cross
    // (cross=0.8) this is the published protocol for running config 3 long on a fixed box.
    S.bc = [=](double x, double y, double, const Vector4d&, double, double, int) -> Vector4d {
        Vector4d pr = quad(x, y);
        return primToCons(pr(0), pr(1), pr(2), pr(3));
    };

    Mesh baseMesh; makeRectMesh(baseMesh, 0, 1, 0, 1, base_ny, base_ny);
    return runAMRScene(baseMesh, S);
}
