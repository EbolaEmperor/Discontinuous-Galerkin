// ===========================================================================
// Shock-bubble interaction (Haas & Sturtevant 1987) with the DG + 2nd-order IMEX
// Euler solver on an adaptive conforming mesh (h-AMR, newest-vertex bisection).
//
// A planar shock crosses a cylindrical pocket of a different-density gas.  The
// density jump across the bubble interface is misaligned with the pressure jump
// across the shock, so baroclinic vorticity (grad rho x grad p) is deposited on
// the interface -- the canonical "shock-accelerated inhomogeneity":
//   * a LIGHT bubble (rho_b < rho_amb, helium-like) acts as a diverging lens and
//     rolls up into a fast counter-rotating vortex PAIR / smoke ring;
//   * a HEAVY bubble (rho_b > rho_amb, R22/SF6-like) acts as a converging lens,
//     focuses the shock, and forms an axial air jet that pierces it.
//
// Single-gamma surrogate: one ideal gas (gamma=1.4); the bubble is a pure DENSITY
// jump at ambient pressure.  This reproduces the baroclinic vortex roll-up but is
// QUALITATIVE -- not the true multi-component He/R22 experiment.  HLLC keeps the
// interface crisp; the PRESSURE artificial-viscosity sensor (av_indicator=1)
// fires on the shock but not the contact, so the bubble interface stays sharp.
//
//   Domain [0, xb] x [0, 1];  ambient air (rho=1,u=0,p=1).  Planar shock Ms at
//   x = x_s0 moving +x (post-shock state from Rankine-Hugoniot).  Bubble radius
//   R centred at (cx, cy), density rho_b = ratio * rho_amb, at rest, ambient p.
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

enum { INTERIOR = 0, LEFT = 1, RIGHT = 2, BOTTOM = 3, TOP = 4 };

int main(int argc, char** argv) {
    std::cout << std::unitbuf;

    // ---- physics / geometry ----
    double domain_xb = 2.5;
    double Ms        = 1.22;    // incident-shock Mach number
    double rho_amb   = 1.0, p_amb = 1.0;
    double x_shock0  = 0.15;    // initial shock position (placed left for downstream room -> long clean run)
    double bub_cx    = 0.4, bub_cy = 0.5, bub_R = 0.2;   // centered classic layout: cx=0.9, x_shock0=0.4
    // rho_b / rho_amb.  HEAVY (R22-like ~3.0) is the default: it has a LOW sound
    // speed, no near-vacuum, so it is robust to long times, and shows the converging
    // lens -> shock focusing -> piercing air jet.  The LIGHT He-like bubble (~0.138,
    // the diverging lens -> counter-rotating vortex ring) is far stiffer (sound speed
    // ~3.2, thin near-vacuum filaments form late); it runs but needs a shorter t_end
    // / smaller cfl and is best treated as a harder variant.
    double ratio     = 3.0;

    AMRScene S;
    S.name = "shock_bubble"; S.framesDir = "out/shock_bubble_frames"; S.title = "Shock-bubble interaction";
    S.ord = 2; S.max_gen = 4; S.t_end = 4.8; S.cfl = 0.4; S.lambda_safe = 3.0; S.n_frames = 600;
    S.remesh_every = 6; S.buffer_layers = 3; S.init_passes = 12; S.th_ref = 0.25; S.th_crs = 0.06;
    // Use the DENSITY sensor (av_indicator=0), NOT the pressure sensor: the bubble
    // interface is a strong CONTACT (a pure density jump at constant p) and the
    // pressure sensor leaves it undamped.  Keep the AV GENTLE (default Persson
    // threshold, modest av_c) -- an over-aggressive low s0 over-fires near the
    // contact and the eps-weighted SIPG loses coercivity (anti-diffusion -> blow-up).
    // For the heavy bubble this is rock-solid; for the light bubble lower cfl/t_end.
    S.av_c = 1.2; S.av_kappa = 1.0; S.sigma_ip = 20.0; S.av_refresh = 5; S.av_indicator = 0;
    S.use_av = true; S.use_pos = true;
    // Cap the CFL wave-speed estimate's density at 0.05 so a spurious near-vacuum node
    // (which the late KH roll-ups spin up) can't blow up |u|+c and collapse dt -- this is
    // what lets the bubble run to long t_end.  The heavy bubble's real density never
    // approaches 0.05, so genuine cells are unaffected.  (Physics floors stay at 1e-12.)
    S.cfl_rho_floor = 0.05;
    S.Wpix = 1600; S.rho_lo = 0.1; S.rho_hi = 5.0; S.cmapName = "inferno";
    int base_ny = 40;

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("shock_bubble_config.json")) cfgPath = "shock_bubble_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        S.ord = cfg.getInt("ord", S.ord);
        base_ny = cfg.getInt("base_ny", base_ny);
        S.max_gen = cfg.getInt("max_gen", S.max_gen);
        domain_xb = cfg.getNumber("domain_xb", domain_xb);
        Ms = cfg.getNumber("Ms", Ms);
        rho_amb = cfg.getNumber("rho_amb", rho_amb);
        p_amb = cfg.getNumber("p_amb", p_amb);
        x_shock0 = cfg.getNumber("x_shock0", x_shock0);
        bub_cx = cfg.getNumber("bubble_cx", bub_cx);
        bub_cy = cfg.getNumber("bubble_cy", bub_cy);
        bub_R = cfg.getNumber("bubble_R", bub_R);
        ratio = cfg.getNumber("bubble_ratio", ratio);
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
        S.rho_floor = cfg.getNumber("rho_floor", S.rho_floor);
        S.p_floor = cfg.getNumber("p_floor", S.p_floor);
        S.cfl_rho_floor = cfg.getNumber("cfl_rho_floor", S.cfl_rho_floor);
        S.Wpix = cfg.getInt("Wpix", S.Wpix);
        S.rho_lo = cfg.getNumber("rho_lo", S.rho_lo);
        S.rho_hi = cfg.getNumber("rho_hi", S.rho_hi);
        S.cmapName = cfg.getString("cmap", S.cmapName);
        S.full_movie = cfg.getBool("full_movie", S.full_movie);
        S.mesh_movie = cfg.getBool("mesh_movie", S.mesh_movie);
        S.framesDir = cfg.getString("frames_dir", S.framesDir);
    }
    S.rxa = 0.0; S.rxb = domain_xb; S.rya = 0.0; S.ryb = 1.0;

    // ---- post-shock state behind a Mach-Ms shock moving +x into ambient air ----
    const double g = GAMMA;
    double Ms2 = Ms * Ms, c1 = std::sqrt(g * p_amb / rho_amb);
    double rho2 = rho_amb * ((g + 1) * Ms2) / ((g - 1) * Ms2 + 2.0);
    double p2   = p_amb * (2.0 * g * Ms2 - (g - 1)) / (g + 1);
    double u2   = (2.0 / (g + 1)) * (Ms - 1.0 / Ms) * c1;
    Vector4d POST = primToCons(rho2, u2, 0.0, p2);

    cout << "Shock-bubble: Ms=" << Ms << "  rho_b/rho_amb=" << ratio
         << (ratio < 1 ? " (LIGHT -> vortex pair / smoke ring)" : " (HEAVY -> focusing + air jet)") << "\n";
    cout << "  post-shock (rho,u,p)=(" << rho2 << "," << u2 << "," << p2 << ")\n";

    double R2 = bub_R * bub_R;
    S.primIC = [=](double x, double y) -> Vector4d {
        if (x < x_shock0) return Vector4d(rho2, u2, 0.0, p2);                       // shocked air
        double dx = x - bub_cx, dy = y - bub_cy;
        if (dx * dx + dy * dy < R2) return Vector4d(ratio * rho_amb, 0.0, 0.0, p_amb);  // bubble (density jump)
        return Vector4d(rho_amb, 0.0, 0.0, p_amb);                                  // ambient air
    };
    S.tagEdge = [=](double mx, double my) -> int {
        if (std::abs(mx - 0.0) < 1e-9) return LEFT;
        if (std::abs(mx - domain_xb) < 1e-9) return RIGHT;
        if (std::abs(my - 0.0) < 1e-9) return BOTTOM;
        if (std::abs(my - 1.0) < 1e-9) return TOP;
        return INTERIOR;
    };
    S.bc = [=](double, double, double, const Vector4d& Uin, double nx_, double ny_, int tg) -> Vector4d {
        switch (tg) {
            case LEFT:  return POST;                                   // post-shock inflow
            case RIGHT: {
                // CHARACTERISTIC (Riemann-invariant) non-reflecting outflow.  A plain
                // zero-gradient ghost reflects at a subsonic boundary, so once the shock
                // reaches the wall (t~1.46) and the post-shock gas (u~0.4 < c~1.26) flows
                // out, the boundary pressure drifts and spawns a rarefaction wedge (the
                // black region).  A FIXED back-pressure fails too: the far-field pressure
                // is ambient (1.0) before the shock but post-shock (1.57) after, so any
                // fixed value injects a spurious wave in the other regime.  Instead take
                // the OUTGOING invariants (entropy s, tangential v, R+=u_n+2c/(g-1)) from
                // the interior and the INCOMING R-=u_n-2c/(g-1) from the ambient far-field,
                // then reconstruct the ghost.  This self-adjusts: ambient->ambient ghost
                // (no spurious compression), post-shock->post-shock (clean exit), and R-
                // is pinned so the pressure can't drift into a rarefaction.
                const double g = GAMMA, gm1 = g - 1.0;
                double rho = Uin(0), u = Uin(1) / rho, v = Uin(2) / rho;   // n=(1,0): u normal, v tangential
                double p = std::max(pressure(Uin), 1e-12);
                double c = std::sqrt(g * p / rho);
                if (u >= c) return Uin;                                // supersonic outflow -> all chars out
                double c_inf = std::sqrt(g * p_amb / rho_amb);         // ambient far-field sound speed
                double Rp = u + 2.0 * c / gm1;                         // outgoing (from interior)
                double Rm = 0.0 - 2.0 * c_inf / gm1;                   // incoming (from ambient far-field, u_inf=0)
                double un = 0.5 * (Rp + Rm);
                double cb = 0.25 * gm1 * (Rp - Rm);
                double s  = p / std::pow(rho, g);                      // entropy from interior (outgoing)
                double rhob = std::pow(std::max(cb * cb, 1e-12) / (g * s), 1.0 / gm1);
                double pb = rhob * cb * cb / g;
                return primToCons(rhob, un, v, pb);                    // tangential v kept from interior
            }
            case BOTTOM:
            case TOP: { double mn = Uin(1) * nx_ + Uin(2) * ny_;       // slip wall
                        return Vector4d(Uin(0), Uin(1) - 2 * mn * nx_, Uin(2) - 2 * mn * ny_, Uin(3)); }
            default: return Uin;
        }
    };

    int base_nx = (int)std::lround(domain_xb * base_ny);
    Mesh baseMesh; makeRectMesh(baseMesh, 0, domain_xb, 0, 1, base_nx, base_ny);
    return runAMRScene(baseMesh, S);
}
