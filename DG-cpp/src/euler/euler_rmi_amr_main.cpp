// ===========================================================================
// Richtmyer-Meshkov instability (RMI) with the DG + 2nd-order IMEX Euler solver
// on an adaptively refined conforming mesh (h-AMR by newest-vertex bisection).
//
// A planar shock crosses a sinusoidally perturbed contact between a light and a
// heavy gas.  The misaligned density and pressure gradients deposit baroclinic
// vorticity (grad rho x grad p) on the interface, so the perturbation inverts
// phase and grows -- first linearly, then into the classic row of mushroom
// "spikes" (heavy into light) and round "bubbles" (light into heavy), with
// secondary Kelvin-Helmholtz roll-ups on the spike stems.  Unlike Rayleigh-
// Taylor it is IMPULSIVELY driven, so it needs NO gravity / body-force term.
//
// Single-gamma surrogate: both gases are the same ideal gas (gamma=1.4); the
// "two gases" are a pure DENSITY (entropy) jump at equal pressure -- the standard
// single-fluid Euler RMI model.  HLLC keeps the contact crisp; the PRESSURE-based
// artificial-viscosity sensor (av_indicator=1) stays off the contact; the
// density-gradient AMR indicator tracks both the shock and the growing interface.
//
// Optional RE-SHOCK variant (reshock=true): the right wall reflects, so the
// transmitted shock bounces back and re-shocks the interface -> violent turbulent
// mixing.  Otherwise the right wall is a transmissive outflow.
//
//   Domain [0, xb] x [0, 1];  light gas rho_L, heavy gas rho_H, both p0, at rest.
//   Interface  x_i(y) = x0 + a0 cos(2 pi y / lambda)  (meets the slip walls at
//   antinodes).  Incident shock at x = x_s0 (in the light gas), Mach Ms in +x.
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

    // ---- physics / geometry parameters ----
    double domain_xb = 4.0;
    double Ms        = 1.5;     // incident-shock Mach number (in the light gas)
    double rho_light = 1.0;
    double rho_heavy = 3.0;     // Atwood A = (rH-rL)/(rH+rL) = 0.5
    double p0        = 1.0;     // ambient pressure (both gases, pre-shock)
    double x_shock0  = 0.3;     // initial incident-shock position (placed left for room to advect)
    double iface_x0  = 0.7;     // mean interface position
    double amp       = 0.1;     // perturbation amplitude
    double wavelength = 1.0;    // perturbation wavelength (= domain height -> 1 period)
    bool   reshock   = false;   // true -> reflecting right wall (transmitted shock re-shocks)

    // ---- scene knobs (defaults tuned for RMI) ----
    AMRScene S;
    S.name = "rmi"; S.framesDir = "out/rmi_frames"; S.title = "Richtmyer-Meshkov instability";
    S.ord = 2; S.max_gen = 6; S.t_end = 3.5; S.cfl = 0.4; S.lambda_safe = 5.0; S.n_frames = 350;
    S.remesh_every = 6; S.buffer_layers = 3; S.init_passes = 12; S.th_ref = 0.30; S.th_crs = 0.08;
    // The interface is a CONTACT.  As AMR sharpens it (max_gen 5) the pressure sensor
    // (which doesn't fire on contacts) leaves it undamped and it oscillates to NaN --
    // so use the DENSITY sensor with a low Persson threshold (av_s0=-3), as for the
    // 2D Riemann configs.  This diffuses the interface a little but the baroclinic
    // mushroom/spike still develops.  (rho stays moderate here, so unlike the very
    // low-density light bubble the aggressive threshold does not break SIPG coercivity.)
    S.av_c = 2.0; S.av_kappa = 1.0; S.sigma_ip = 20.0; S.av_refresh = 5; S.av_indicator = 0;
    S.av_s0 = -3.0; S.av_s0_set = true;
    S.use_av = true; S.use_pos = true;
    S.Wpix = 1600; S.rho_lo = 0.5; S.rho_hi = 6.0; S.cmapName = "viridis";
    S.rxa = 0.0; S.rxb = domain_xb; S.rya = 0.0; S.ryb = 1.0;
    int base_ny = 25;

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("rmi_config.json")) cfgPath = "rmi_config.json";
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
        rho_light = cfg.getNumber("rho_light", rho_light);
        rho_heavy = cfg.getNumber("rho_heavy", rho_heavy);
        p0 = cfg.getNumber("p0", p0);
        x_shock0 = cfg.getNumber("x_shock0", x_shock0);
        iface_x0 = cfg.getNumber("iface_x0", iface_x0);
        amp = cfg.getNumber("amp", amp);
        wavelength = cfg.getNumber("wavelength", wavelength);
        reshock = cfg.getBool("reshock", reshock);
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
        S.av_kappa = cfg.getNumber("av_kappa", S.av_kappa);
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
    S.rxa = 0.0; S.rxb = domain_xb; S.rya = 0.0; S.ryb = 1.0;

    // ---- post-shock state behind a Mach-Ms shock moving +x into the light gas ----
    const double g = GAMMA;
    double Ms2 = Ms * Ms;
    double c1  = std::sqrt(g * p0 / rho_light);
    double rho2 = rho_light * ((g + 1) * Ms2) / ((g - 1) * Ms2 + 2.0);
    double p2   = p0 * (2.0 * g * Ms2 - (g - 1)) / (g + 1);
    double u2   = (2.0 / (g + 1)) * (Ms - 1.0 / Ms) * c1;     // post-shock particle velocity (+x)
    Vector4d POST = primToCons(rho2, u2, 0.0, p2);
    Vector4d PRE_L = primToCons(rho_light, 0.0, 0.0, p0);     // (unused; documents states)
    (void)PRE_L;

    cout << "Richtmyer-Meshkov: Ms=" << Ms << "  light rho=" << rho_light << " heavy rho=" << rho_heavy
         << " (A=" << (rho_heavy - rho_light) / (rho_heavy + rho_light) << ")\n";
    cout << "  post-shock (rho,u,p)=(" << rho2 << "," << u2 << "," << p2 << ")  "
         << (u2 + std::sqrt(g * p2 / rho2) > 0 && u2 > std::sqrt(g * p2 / rho2) ? "supersonic" : "subsonic")
         << " inflow;  re-shock=" << (reshock ? "ON" : "off") << "\n";

    // ---- scene hooks ----
    const double TWO_PI = 2.0 * M_PI;
    S.primIC = [=](double x, double y) -> Vector4d {
        if (x < x_shock0) return Vector4d(rho2, u2, 0.0, p2);            // shocked light gas (primitive)
        double xi = iface_x0 + amp * std::cos(TWO_PI * y / wavelength);  // perturbed interface
        if (x < xi) return Vector4d(rho_light, 0.0, 0.0, p0);           // unshocked light gas
        return Vector4d(rho_heavy, 0.0, 0.0, p0);                       // heavy gas
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
            case LEFT:  return POST;                                     // post-shock inflow
            case RIGHT:
                if (!reshock) return Uin;                               // transmissive outflow
                /* fallthrough: reflecting wall for the re-shock variant */
            case BOTTOM:
            case TOP:
                return slipWallExterior(Uin, nx_, ny_);
            default: return Uin;
        }
    };

    int base_nx = (int)std::lround(domain_xb * base_ny);
    Mesh baseMesh; makeRectMesh(baseMesh, 0, domain_xb, 0, 1, base_nx, base_ny);
    return runAMRScene(baseMesh, S);
}
