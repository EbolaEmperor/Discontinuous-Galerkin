// ===========================================================================
// Shock diffraction over a 90-degree convex corner (Zhang & Shu 2010 benchmark)
// with the DG + 2nd-order IMEX Euler solver on an adaptive conforming mesh.
//
// A strong (Mach 5.09) shock travels down a narrow inlet channel and reaches a
// backward-facing step: at the convex corner the flow expands around the 90-deg
// turn, the shock diffracts and curves, a large primary vortex rolls up off the
// corner, a shear layer goes Kelvin-Helmholtz, and -- crucially -- a NEAR-VACUUM
// region forms right at the corner.  That near-vacuum is THE canonical test of a
// positivity-preserving scheme: this solver's Zhang-Shu limiter is exactly what
// keeps rho>0 and p>0 there.  AMR tracks the moving diffracting shock cheaply.
//
//   Domain  Omega = ([0,1] x [6,11]) U ([1,13] x [0,11])   (an L / backward step;
//   the full box [0,13]x[0,11] minus the lower-left hole [0,1]x[0,6]).
//   Corner at (1,6).  Pre-shock (rho,u,v,p)=(1.4,0,0,1); a right-moving Mach-5.09
//   shock starts at x=0.5 in the channel.  BCs: post-shock inflow on the channel
//   left (x=0, 6<=y<=11); reflecting walls on the top (y=11) and the two corner
//   walls (y=6 for 0<=x<=1, x=1 for 0<=y<=6); transmissive outflow on the right
//   (x=13) and bottom (y=0).
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

enum { INTERIOR = 0, INFLOW = 1, OUTFLOW = 2, WALL = 3 };

int main(int argc, char** argv) {
    std::cout << std::unitbuf;

    // ---- geometry (Zhang-Shu corner diffraction) ----
    double Lx = 13.0, Ly = 11.0;          // full bounding box
    double chan_w = 1.0, chan_y0 = 6.0;   // channel occupies [0,chan_w] x [chan_y0, Ly]; hole = [0,chan_w]x[0,chan_y0]
    double Ms = 5.09;                      // incident-shock Mach number
    double rho_pre = 1.4, p_pre = 1.0;     // pre-shock (ambient at rest)
    double x_shock0 = 0.5;                 // initial shock position in the channel
    double base_cell = 0.25;               // base square-cell size (Lx/cell, Ly/cell integers; corner grid-aligned)

    AMRScene S;
    S.name = "corner_diffraction"; S.framesDir = "out/corner_frames"; S.title = "Shock diffraction over a 90-degree corner";
    S.ord = 2; S.max_gen = 5; S.t_end = 2.0; S.cfl = 0.3; S.lambda_safe = 8.0; S.n_frames = 240;
    S.remesh_every = 6; S.buffer_layers = 3; S.init_passes = 12; S.th_ref = 0.30; S.th_crs = 0.08;
    S.av_c = 1.2; S.av_kappa = 1.0; S.sigma_ip = 20.0; S.av_refresh = 5; S.av_indicator = 1;
    S.use_av = true; S.use_pos = true;
    S.Wpix = 1600; S.rho_lo = 0.05; S.rho_hi = 8.0; S.cmapName = "inferno";

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("corner_config.json")) cfgPath = "corner_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        S.ord = cfg.getInt("ord", S.ord);
        base_cell = cfg.getNumber("base_cell", base_cell);
        S.max_gen = cfg.getInt("max_gen", S.max_gen);
        Lx = cfg.getNumber("Lx", Lx); Ly = cfg.getNumber("Ly", Ly);
        chan_w = cfg.getNumber("chan_w", chan_w); chan_y0 = cfg.getNumber("chan_y0", chan_y0);
        Ms = cfg.getNumber("Ms", Ms);
        rho_pre = cfg.getNumber("rho_pre", rho_pre);
        p_pre = cfg.getNumber("p_pre", p_pre);
        x_shock0 = cfg.getNumber("x_shock0", x_shock0);
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
        S.Wpix = cfg.getInt("Wpix", S.Wpix);
        S.rho_lo = cfg.getNumber("rho_lo", S.rho_lo);
        S.rho_hi = cfg.getNumber("rho_hi", S.rho_hi);
        S.cmapName = cfg.getString("cmap", S.cmapName);
        S.full_movie = cfg.getBool("full_movie", S.full_movie);
        S.mesh_movie = cfg.getBool("mesh_movie", S.mesh_movie);
        S.framesDir = cfg.getString("frames_dir", S.framesDir);
    }
    S.rxa = 0.0; S.rxb = Lx; S.rya = 0.0; S.ryb = Ly;

    // ---- post-shock state behind a Mach-Ms shock moving +x into (rho_pre,p_pre) at rest ----
    const double g = GAMMA;
    double Ms2 = Ms * Ms, c1 = std::sqrt(g * p_pre / rho_pre);
    double rho2 = rho_pre * ((g + 1) * Ms2) / ((g - 1) * Ms2 + 2.0);
    double p2   = p_pre * (2.0 * g * Ms2 - (g - 1)) / (g + 1);
    double u2   = (2.0 / (g + 1)) * (Ms - 1.0 / Ms) * c1;
    Vector4d POST = primToCons(rho2, u2, 0.0, p2);

    cout << "Corner diffraction: Mach " << Ms << " shock; post-shock (rho,u,p)=("
         << rho2 << "," << u2 << "," << p2 << ")  c2=" << std::sqrt(g * p2 / rho2)
         << (u2 > std::sqrt(g * p2 / rho2) ? " (supersonic inflow)" : " (subsonic inflow)") << "\n";

    S.primIC = [=](double x, double y) -> Vector4d {
        if (x < x_shock0) return Vector4d(rho2, u2, 0.0, p2);          // shocked gas (only exists in the channel)
        return Vector4d(rho_pre, 0.0, 0.0, p_pre);                     // pre-shock gas at rest
    };
    // edge tagger over the L-shaped boundary
    S.tagEdge = [=](double mx, double my) -> int {
        const double eps = 1e-7;
        if (std::abs(mx - 0.0) < eps) return INFLOW;                                 // channel left wall
        if (std::abs(mx - Lx) < eps)  return OUTFLOW;                                // right
        if (std::abs(my - 0.0) < eps) return OUTFLOW;                                // bottom
        if (std::abs(my - Ly) < eps)  return WALL;                                   // top
        if (std::abs(my - chan_y0) < eps && mx < chan_w) return WALL;               // channel floor (y=6, x<1)
        if (std::abs(mx - chan_w) < eps && my < chan_y0) return WALL;               // step face (x=1, y<6)
        return INTERIOR;
    };
    S.bc = [=](double, double, double, const Vector4d& Uin, double nx_, double ny_, int tg) -> Vector4d {
        switch (tg) {
            case INFLOW:  return POST;
            case OUTFLOW: return Uin;
            case WALL: { double mn = Uin(1) * nx_ + Uin(2) * ny_;
                         return Vector4d(Uin(0), Uin(1) - 2 * mn * nx_, Uin(2) - 2 * mn * ny_, Uin(3)); }
            default: return Uin;
        }
    };
    // render mask: blank out the lower-left hole (no elements there anyway)
    S.inDomain = [=](double x, double y) -> bool { return !(x < chan_w && y < chan_y0); };

    // Square base cells of size 1/k.  For the standard integer geometry (Lx,Ly,
    // chan_w,chan_y0 = 13,11,1,6) this guarantees the corner (chan_w, chan_y0) and
    // the whole hole boundary fall exactly on grid lines, so deleting the hole's
    // cells leaves a conforming criss-cross base (no hanging nodes).
    int k = std::max(1, (int)std::lround(1.0 / base_cell));
    int nx = (int)std::lround(Lx * k), ny = (int)std::lround(Ly * k);
    auto keep = [=](double x, double y) -> bool { return !(x < chan_w && y < chan_y0); };
    Mesh baseMesh; makeMaskedRectMesh(baseMesh, 0, Lx, 0, Ly, nx, ny, keep);
    cout << "  base L-mesh: " << baseMesh.elem.rows() << " tris (cell " << 1.0 / k << ", "
         << nx << "x" << ny << " grid minus hole)\n";
    return runAMRScene(baseMesh, S);
}
