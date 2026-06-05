// ===========================================================================
// Double Mach Reflection (Mach 10) with the DG + 2nd-order IMEX Euler solver,
// driven on an ADAPTIVELY REFINED conforming mesh (h-AMR by newest-vertex
// bisection).  Same physics, fluxes, artificial viscosity, positivity limiter
// and ARS(2,2,2) time stepping as euler_dmr -- only the mesh is now adaptive:
//
//   * a scaled density-gradient + inter-cell-jump indicator refines the
//     incident / reflected / Mach-stem shocks AND the slip line (the density
//     jump there drives the Kelvin-Helmholtz vortex street) to a target level,
//     and coarsens the smooth post-shock / undisturbed regions, with hysteresis
//     (refine at th_ref > coarsen at th_crs) and an N-layer buffer so the moving
//     shocks cannot outrun the fine zone between remeshes;
//   * the mesh stays CONFORMING (NVB) so the entire proven DG machinery is
//     reused unchanged; the DG state is transferred between meshes exactly on
//     refine (polynomial restriction) and conservatively on coarsen (L2);
//   * CFL CONTROL: after every remesh the global time step is recomputed from
//     the current finest cell, dt = cfl * h_min / ((2k+1) * lambda_safe), with
//     h_min = min over leaves of 2*area/longest_edge (altitude -- robust to the
//     bisection aspect ratio), so refinement automatically shrinks dt.
//
// The win is FEW ELEMENTS: the features are 1-D curves, so an adaptive mesh
// resolves them at the same effective resolution as a fine uniform grid with a
// small fraction of the cells (and hence per-step cost), matching the uniform
// ny=150 result on a fraction of the work.
//
// Output: a density movie + a mesh-overlay movie (the adaptive grid drawn over
// the density) + final density / Schlieren stills.  All parameters live in JSON.
// ===========================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "AdaptiveMesh.h"
#include "EulerDG.h"
#include "FEM.h"
#include "Json.h"
#include "Mesh.h"

using namespace Eigen;
using namespace std;
using namespace euler;
namespace fs = std::filesystem;

enum { INTERIOR = 0, LEFT = 1, RIGHT = 2, BOTTOM = 3, TOP = 4 };

// ---- per-element geometry helpers ----------------------------------------
static double longestEdge(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0)), p1 = m.node.row(m.elem(t, 1)), p2 = m.node.row(m.elem(t, 2));
    return std::max({(p1 - p0).norm(), (p2 - p1).norm(), (p0 - p2).norm()});
}
static double triArea(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0)), p1 = m.node.row(m.elem(t, 1)), p2 = m.node.row(m.elem(t, 2));
    return 0.5 * std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) - (p2.x() - p0.x()) * (p1.y() - p0.y()));
}
// CFL length scale = altitude to the longest edge (robust to thin children).
static double hCFL(const Mesh& m, int t) { return 2.0 * triArea(m, t) / std::max(longestEdge(m, t), 1e-300); }

// ---- boundary tags by geometry (rebuilt every remesh) --------------------
static VectorXi computeTags(const Mesh& mesh, const MatrixXi& edge, const MatrixXi& e2s, double xb) {
    int NE = static_cast<int>(edge.rows());
    VectorXi tag = VectorXi::Zero(NE);
    for (int e = 0; e < NE; ++e) {
        if (e2s(e, 0) != -1 && e2s(e, 1) != -1) continue;
        Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1))).transpose();
        if (std::abs(m.x() - 0.0) < 1e-9) tag(e) = LEFT;
        else if (std::abs(m.x() - xb) < 1e-9) tag(e) = RIGHT;
        else if (std::abs(m.y() - 0.0) < 1e-9) tag(e) = BOTTOM;
        else if (std::abs(m.y() - 1.0) < 1e-9) tag(e) = TOP;
    }
    return tag;
}

// ---- strong geometric conformity test: every 1-sided edge is on the box ----
static bool checkGeomConforming(const Mesh& mesh, const MatrixXi& edge, const MatrixXi& e2s,
                                double xb, int& nHanging) {
    nHanging = 0;
    for (int e = 0; e < edge.rows(); ++e) {
        bool oneSided = (e2s(e, 0) == -1 || e2s(e, 1) == -1);
        if (!oneSided) continue;
        Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1))).transpose();
        bool onBox = std::abs(m.x()) < 1e-9 || std::abs(m.x() - xb) < 1e-9 ||
                     std::abs(m.y()) < 1e-9 || std::abs(m.y() - 1.0) < 1e-9;
        if (!onBox) ++nHanging;   // a 1-sided edge in the interior == a hanging node
    }
    return nHanging == 0;
}

// ---- minimal PPM read/write + mesh-edge overlay --------------------------
static bool readPPM(const string& path, int& W, int& H, vector<unsigned char>& img) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    string magic; in >> magic; if (magic != "P6") return false;
    int maxv; in >> W >> H >> maxv; in.get();
    img.resize((size_t)W * H * 3);
    in.read(reinterpret_cast<char*>(img.data()), (std::streamsize)img.size());
    return true;
}
static void writePPM(const string& path, int W, int H, const vector<unsigned char>& img) {
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}
static void drawLine(vector<unsigned char>& img, int W, int H, int x0, int y0, int x1, int y1,
                     unsigned char r, unsigned char g, unsigned char b, double alpha) {
    int dx = std::abs(x1 - x0), dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1, err = dx + dy;
    while (true) {
        if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) {
            size_t k = ((size_t)y0 * W + x0) * 3;
            img[k]   = (unsigned char)(alpha * r + (1 - alpha) * img[k]);
            img[k+1] = (unsigned char)(alpha * g + (1 - alpha) * img[k+1]);
            img[k+2] = (unsigned char)(alpha * b + (1 - alpha) * img[k+2]);
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}
// Draw the triangulation edges onto an existing density PPM (one frame).
static void overlayMesh(const string& path, const Mesh& mesh, double xa, double xb, double ya, double yb) {
    int W, H; vector<unsigned char> img;
    if (!readPPM(path, W, H, img)) return;
    auto toPix = [&](const Vector2d& p, int& px, int& py) {
        px = (int)std::lround((p.x() - xa) / (xb - xa) * W);
        py = (int)std::lround((yb - p.y()) / (yb - ya) * H);
    };
    for (int t = 0; t < mesh.elem.rows(); ++t) {
        Vector2d c = (mesh.node.row(mesh.elem(t,0)) + mesh.node.row(mesh.elem(t,1)) + mesh.node.row(mesh.elem(t,2))).transpose() / 3.0;
        if (c.x() < xa - 0.02 || c.x() > xb + 0.02 || c.y() < ya - 0.02 || c.y() > yb + 0.02) continue;
        for (int k = 0; k < 3; ++k) {
            Vector2d A = mesh.node.row(mesh.elem(t, k)), B = mesh.node.row(mesh.elem(t, (k + 1) % 3));
            int ax, ay, bx, by; toPix(A, ax, ay); toPix(B, bx, by);
            drawLine(img, W, H, ax, ay, bx, by, 80, 200, 255, 0.32);   // faint cyan grid
        }
    }
    writePPM(path, W, H, img);
}

int main(int argc, char** argv) {
    std::cout << std::unitbuf;
    // --------------------------- parameters (JSON) ---------------------------
    int    ord        = 2;
    int    base_ny    = 50;        // base-mesh cells in y  (h_base = 1/base_ny)
    int    max_gen    = 4;         // max NVB bisections from base (h ~ h_base * 2^(-gen/2))
    double domain_xb  = 4.5;
    double t_end      = 0.30;
    double cfl        = 0.40;
    double lambda_safe = 19.0;
    int    n_frames   = 300;
    // AMR controls
    int    remesh_every = 6;       // R: steps between remeshes
    int    buffer_layers = 3;      // N: refine-buffer dilation around flagged cells
    double th_ref     = 0.30;      // refine if indicator > th_ref
    double th_crs     = 0.08;      // coarsen if indicator < th_crs  (hysteresis: th_crs<th_ref)
    int    init_passes = 12;       // initial-adaptation sweeps (refine the IC shock to max_gen)
    // artificial viscosity / limiter (match euler_dmr defaults)
    double av_c       = 1.2;
    double av_kappa   = 1.0;
    double sigma_ip   = 20.0;
    int    av_refresh = 5;
    bool   use_av     = true;
    bool   use_pos    = true;
    // rendering
    int    Wpix       = 1920;
    double rxa = 0.0, rxb = 4.4, rya = 0.0, ryb = 1.0;
    double rho_lo = 1.4, rho_hi = 21.0;
    double zxa = 3.32, zxb = 4.12, zya = 0.0, zyb = 0.42;
    bool   full_movie = true, mesh_movie = true, save_state = false;
    bool   check_conservation = false;   // opt-in: verify transfer conservation each remesh (2 extra integrals)
    string state_only = "";
    string framesDir = "dmr_amr_frames", cmapName = "inferno";

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("dmr_amr_config.json")) cfgPath = "dmr_amr_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        ord = cfg.getInt("ord", ord);
        base_ny = cfg.getInt("base_ny", base_ny);
        max_gen = cfg.getInt("max_gen", max_gen);
        domain_xb = cfg.getNumber("domain_xb", domain_xb);
        t_end = cfg.getNumber("t_end", t_end);
        cfl = cfg.getNumber("cfl", cfl);
        lambda_safe = cfg.getNumber("lambda_safe", lambda_safe);
        n_frames = cfg.getInt("n_frames", n_frames);
        remesh_every = cfg.getInt("remesh_every", remesh_every);
        buffer_layers = cfg.getInt("buffer_layers", buffer_layers);
        th_ref = cfg.getNumber("th_ref", th_ref);
        th_crs = cfg.getNumber("th_crs", th_crs);
        init_passes = cfg.getInt("init_passes", init_passes);
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
        full_movie = cfg.getBool("full_movie", full_movie);
        mesh_movie = cfg.getBool("mesh_movie", mesh_movie);
        save_state = cfg.getBool("save_state", save_state);
        check_conservation = cfg.getBool("check_conservation", check_conservation);
        state_only = cfg.getString("state_only", state_only);
        framesDir = cfg.getString("frames_dir", framesDir);
        cmapName = cfg.getString("cmap", cmapName);
    }
    int locDof = (ord + 1) * (ord + 2) / 2;
    int base_nx = (int)std::lround(domain_xb * base_ny);
    Colormap cmap = (cmapName == "jet") ? CM_JET : (cmapName == "viridis") ? CM_VIRIDIS
                  : (cmapName == "gray") ? CM_GRAY : CM_INFERNO;

    // --------------------------- states (Woodward-Colella) ---------------------------
    const double SQ3 = std::sqrt(3.0);
    auto preP  = [] { return Vector4d(1.4, 0.0, 0.0, 1.0); };
    auto postP = [] { return Vector4d(8.0, 8.25 * std::cos(M_PI / 6.0), -8.25 * std::sin(M_PI / 6.0), 116.5); };
    Vector4d POST = primToCons(postP()(0), postP()(1), postP()(2), postP()(3));
    Vector4d PRE  = primToCons(preP()(0),  preP()(1),  preP()(2),  preP()(3));
    auto initPrim = [&](double x, double y) { return (x < 1.0/6.0 + y / SQ3) ? postP() : preP(); };

    cout << "DG + 2nd-order IMEX Euler -- Double Mach Reflection (Mach 10) with h-AMR (NVB)\n";
    cout << "  config: " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  domain [0," << domain_xb << "]x[0,1]  dP" << ord << "  base nx=" << base_nx << " ny=" << base_ny
         << "  max_gen=" << max_gen << "  (finest h ~ " << (1.0/base_ny)/std::pow(2.0, max_gen/2.0) << ")\n";
    cout << "  AMR: th_ref=" << th_ref << " th_crs=" << th_crs << " buffer=" << buffer_layers
         << " remesh_every=" << remesh_every << "\n";
    cout << "  AV=" << (use_av?"on":"off") << " (c=" << av_c << ", sigma=" << sigma_ip << ", refresh=" << av_refresh
         << ")  positivity=" << (use_pos?"on":"off") << "  t_end=" << t_end << "\n";

    EulerConfig dgcfg; dgcfg.use_av = use_av; dgcfg.use_positivity = use_pos;
    dgcfg.av_c = av_c; dgcfg.av_kappa = av_kappa; dgcfg.sigma_ip = sigma_ip; dgcfg.av_refresh = av_refresh;

    // --------------------------- base mesh + forest ---------------------------
    Mesh baseMesh; makeRectMesh(baseMesh, 0, domain_xb, 0, 1, base_nx, base_ny);
    AdaptiveForest forest(baseMesh, ord);

    // driver-owned, persistent objects (EulerDG references these; must outlive it)
    Mesh mesh;
    std::unique_ptr<FEM> fem;
    MatrixXi e2d, edge, e2s; int nDof = 0; VectorXi tag;
    std::unique_ptr<EulerDG> dg;

    // ---- lightweight profiling (AMR_PROFILE=1): where does the wall time go? ----
    const bool PROFILE = std::getenv("AMR_PROFILE") != nullptr;
    double Tstep = 0, Tremesh = 0, Trender = 0;                 // top-level buckets
    double Tbuild = 0, Tfem = 0, Tdof = 0, Tedge = 0, Tdgctor = 0;  // rebuildSolver sub-parts
    double Tflags = 0, Tadapt = 0, Tgs = 0, Tcons = 0;         // other remesh sub-parts
    auto clk = [] { return std::chrono::high_resolution_clock::now(); };
    auto secs = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };

    // Rebuild fem / dofs / edges / tags / solver for the CURRENT forest mesh.
    auto rebuildSolver = [&]() {
        auto a = clk(); forest.buildMesh(mesh);                 if (PROFILE) Tbuild   += secs(a, clk());
        a = clk(); fem = std::make_unique<FEM>(ord, mesh, /*withHessian=*/false);  if (PROFILE) Tfem += secs(a, clk());
        a = clk(); fem->getDOF(mesh, e2d, nDof);                if (PROFILE) Tdof     += secs(a, clk());
        a = clk(); mesh.getEdge2Side(edge, e2s);
                   tag = computeTags(mesh, edge, e2s, domain_xb); if (PROFILE) Tedge  += secs(a, clk());
        a = clk(); dg = std::make_unique<EulerDG>(*fem, mesh, e2d, edge, e2s, tag, dgcfg);
                                                                if (PROFILE) Tdgctor += secs(a, clk());
    };

    // ----- indicator: scaled density gradient + inter-cell jump -> per-cell flag -----
    auto computeFlags = [&](bool allowCoarsen) {
        int NT = static_cast<int>(mesh.elem.rows());
        const MatrixXd& U = dg->state();
        VectorXd rbar(NT), eta(NT);
        // |grad rho| at the centroid only (the inter-cell jump term below catches
        // shocks sitting on a cell face). The reference d/dlam of the basis at the
        // centroid is order-only -> evaluate ONCE; per element grad = Dlam_t * (dphiC
        // * rho_block). Avoids a per-element-per-point polynomial-gradient build
        // (std::pow + allocs) -- the single biggest remesh cost in the baseline.
        MatrixXd dphiC = fem->computeBasisDlam_all(Vector3d(1.0/3, 1.0/3, 1.0/3));   // 3 x locDof
        VectorXd rblk(locDof);
        for (int t = 0; t < NT; ++t) {
            double rb = 0;
            for (int i = 0; i < locDof; ++i) { rblk(i) = U(e2d(t, i), 0); rb += rblk(i); }
            rb /= locDof; rbar(t) = rb;
            Vector2d g = fem->Dlam[t] * (dphiC * rblk);                   // physical grad rho at centroid
            eta(t) = longestEdge(mesh, t) * g.norm() / std::max(rb, 1e-9);  // ~ relative dr/r across the cell
        }
        // inter-cell density-jump term (catches shocks sitting on a cell face / wall)
        for (int e = 0; e < e2s.rows(); ++e) {
            int a = e2s(e, 0), b = e2s(e, 1);
            if (a < 0 || b < 0) continue;
            double j = std::abs(rbar(a) - rbar(b)) / std::max(std::min(rbar(a), rbar(b)), 1e-9);
            eta(a) = std::max(eta(a), j); eta(b) = std::max(eta(b), j);
        }
        // base flags with hysteresis
        std::vector<int> flag(NT, 0);
        for (int t = 0; t < NT; ++t) {
            if (eta(t) > th_ref && forest.gen(t) < max_gen) flag[t] = 1;
            else if (allowCoarsen && eta(t) < th_crs && forest.gen(t) > 0) flag[t] = -1;
        }
        // buffer: dilate the refine set by N face-neighbour layers (overrides coarsen/hold)
        std::vector<std::vector<int>> nbr(NT);
        for (int e = 0; e < e2s.rows(); ++e) {
            int a = e2s(e, 0), b = e2s(e, 1);
            if (a >= 0 && b >= 0) { nbr[a].push_back(b); nbr[b].push_back(a); }
        }
        std::vector<int> front;
        for (int t = 0; t < NT; ++t) if (flag[t] == 1) front.push_back(t);
        for (int layer = 0; layer < buffer_layers; ++layer) {
            std::vector<int> next;
            for (int t : front) for (int n : nbr[t])
                if (flag[n] != 1 && forest.gen(n) < max_gen) { flag[n] = 1; next.push_back(n); }
            front.swap(next);
            if (front.empty()) break;
        }
        // a cell adjacent to a refine cell must not coarsen (keep the grading smooth)
        for (int t = 0; t < NT; ++t)
            if (flag[t] == -1)
                for (int n : nbr[t]) if (flag[n] == 1) { flag[t] = 0; break; }
        return flag;
    };

    auto hMinCFL = [&]() {
        double h = 1e300; for (int t = 0; t < mesh.elem.rows(); ++t) h = std::min(h, hCFL(mesh, t));
        return h;
    };
    // CFL control: size dt from the ACTUAL global max wave speed (|u|+c) each remesh,
    // with a safety margin, floored by lambda_safe. The fixed-lambda_safe estimate
    // underestimates the real peak |u|+c (~25 here) and lets the effective CFL drift
    // up; using the measured speed pins the effective CFL near `cfl` -- robust when
    // refinement shrinks h_min (the user's "control CFL on refinement" requirement).
    auto dtFromCFL = [&]() {
        double lam = std::max(lambda_safe, 1.15 * dg->maxWaveSpeedGlobal());
        return cfl * hMinCFL() / ((2.0 * ord + 1.0) * lam);
    };
    auto maxAspect = [&]() {
        double a = 1.0; for (int t = 0; t < mesh.elem.rows(); ++t)
            a = std::max(a, longestEdge(mesh, t) / std::max(hCFL(mesh, t), 1e-300));
        return a;
    };

    // BC ghost states (mesh-independent; defined once)
    ExteriorStateFn bc = [&](double x, double y, double t, const Vector4d& Uin, double nx_, double ny_, int tg) -> Vector4d {
        switch (tg) {
            case LEFT:  return POST;
            case RIGHT: return Uin;
            case BOTTOM:
                if (x < 1.0/6.0) return POST;
                else { double mn = Uin(1)*nx_ + Uin(2)*ny_; return Vector4d(Uin(0), Uin(1)-2*mn*nx_, Uin(2)-2*mn*ny_, Uin(3)); }
            case TOP: { double xs = 1.0/6.0 + (1.0 + 20.0 * t) / SQ3; return (x < xs) ? POST : PRE; }
            default: return Uin;
        }
    };

    // ----------------------- rendering helpers -----------------------
    if (Wpix % 2) ++Wpix;
    auto writeDensity = [&](const string& path, double xa, double xb, double ya, double yb, int W) {
        int Hh = std::max(2, (int)std::lround(W * (yb - ya) / (xb - xa))); if (Hh%2) ++Hh;
        writeScalarPPM(path, *fem, mesh, e2d, dg->densityField(), W, Hh, xa, xb, ya, yb, rho_lo, rho_hi, cmap, {});
    };
    auto writeSchlieren = [&](const string& path, double xa, double xb, double ya, double yb, int W) {
        int Hh = std::max(2, (int)std::lround(W * (yb - ya) / (xb - xa))); if (Hh%2) ++Hh;
        writeSchlierenPPM(path, *fem, mesh, e2d, dg->densityField(), W, Hh, xa, xb, ya, yb, 8.0);
    };
    auto renderStills = [&]() {
        writeDensity ("dmr_amr_density.ppm",      rxa, rxb, rya, ryb, Wpix);
        writeDensity ("dmr_amr_density_zoom.ppm", zxa, zxb, zya, zyb, Wpix);
        writeSchlieren("dmr_amr_schlieren.ppm",    rxa, rxb, rya, ryb, Wpix);
        writeSchlieren("dmr_amr_schlieren_zoom.ppm", zxa, zxb, zya, zyb, Wpix);
        string meshStill = "dmr_amr_mesh.ppm";
        writeDensity(meshStill, rxa, rxb, rya, ryb, Wpix);
        overlayMesh(meshStill, mesh, rxa, rxb, rya, ryb);
    };

    // ----------------------- state-only mode -----------------------
    if (!state_only.empty()) {
        // rebuild the SAVED mesh + state (saved as forest is not persisted; we just
        // reload a uniform render of the saved nodal state on a matching mesh).
        cout << "  [state_only] not supported for AMR (mesh is adaptive); rerun the solver.\n";
        return 1;
    }

    // =======================================================================
    // Transfer self-test (AMR_SELFTEST=1): a degree-1 conservative field must be
    // preserved EXACTLY through refine-all -> coarsen-all -- both the integral
    // (conservation) AND the L2 energy (exact restriction + exact projection).
    // =======================================================================
    if (std::getenv("AMR_SELFTEST")) {
        // FULL degree-2 field (constant velocity keeps primToCons within dP2): rho and
        // p are quadratic in x,y so all 6 dP2 basis modes (1,x,y,x^2,xy,y^2) participate
        // -> the test actually exercises the quadratic rows of the transfer matrices.
        auto smooth = [](double x, double y) {
            return Vector4d(2.0 + 0.3*x + 0.2*y + 0.10*x*x - 0.05*x*y + 0.02*y*y, 0.5, -0.3,
                            5.0 + 0.4*x + 0.1*y + 0.08*x*x + 0.03*x*y - 0.04*y*y);
        };
        rebuildSolver();
        dg->setState(projectInitial(*fem, mesh, e2d, smooth));
        auto integ = [&](int comp, int pw) {
            const MatrixXd& U = dg->state(); MatrixXd qL; VectorXd w; fem->quad2d(qL, w);
            MatrixXd phi = fem->computeBasisValue_all(qL); int nq = (int)w.size(); double s = 0;
            for (int tt = 0; tt < mesh.elem.rows(); ++tt) {
                double ar = triArea(mesh, tt);
                for (int q = 0; q < nq; ++q) {
                    double val = 0; for (int i = 0; i < locDof; ++i) val += phi(q, i) * U(e2d(tt, i), comp);
                    s += w(q) * ar * std::pow(val, pw);
                }
            }
            return s;
        };
        auto flagAll = [&](int v) { return std::vector<int>(mesh.elem.rows(), v); };
        double m0 = integ(0, 1), q0 = integ(0, 2);
        for (int p = 0; p < max_gen + 2; ++p) {
            forest.syncFromState(dg->state(), e2d, locDof); auto fl = flagAll(1);
            dg.reset(); auto rc = forest.adapt(fl, max_gen); rebuildSolver();
            dg->setState(forest.gatherState(e2d, locDof, nDof));
            if (rc.first == 0) break;
        }
        double m1 = integ(0, 1), q1 = integ(0, 2); int nfine = (int)mesh.elem.rows();
        for (int p = 0; p < max_gen + 2; ++p) {
            forest.syncFromState(dg->state(), e2d, locDof); auto fl = flagAll(-1);
            dg.reset(); auto rc = forest.adapt(fl, max_gen); rebuildSolver();
            dg->setState(forest.gatherState(e2d, locDof, nDof));
            if (rc.second == 0) break;
        }
        double m2 = integ(0, 1), q2 = integ(0, 2);
        cout << "\n[AMR_SELFTEST] degree-2 field through refine-all -> coarsen-all\n";
        cout << "  refine: " << nfine << " tris   coarsen back to " << mesh.elem.rows() << " tris\n";
        cout << scientific << setprecision(3);
        cout << "  int(rho) drift  refine=" << std::abs(m1-m0)/std::abs(m0)
             << "  cycle=" << std::abs(m2-m0)/std::abs(m0) << "   (conservation)\n";
        cout << "  int(rho^2) drift refine=" << std::abs(q1-q0)/std::abs(q0)
             << "  cycle=" << std::abs(q2-q0)/std::abs(q0) << "   (exactness)\n";
        return 0;
    }

    // =======================================================================
    // Initial condition + initial adaptation (refine the IC shock to max_gen).
    // The analytic IC is RE-PROJECTED onto each refined mesh (so no resolution
    // is lost to the coarse base), then the indicator drives the next refine.
    // =======================================================================
    rebuildSolver();
    dg->setState(projectInitial(*fem, mesh, e2d, initPrim));
    cout << "\nInitial adaptation:\n";
    for (int pass = 0; pass < init_passes; ++pass) {
        std::vector<int> flag = computeFlags(/*allowCoarsen=*/false);
        int nref = 0; for (int f : flag) if (f == 1) ++nref;
        if (nref == 0) break;
        forest.syncFromState(dg->state(), e2d, locDof);     // (unused for IC; we re-project)
        dg.reset();                                          // free refs before mutating mesh
        auto [r, c] = forest.adapt(flag, max_gen);
        rebuildSolver();
        dg->setState(projectInitial(*fem, mesh, e2d, initPrim));   // sharp IC on the new mesh
        string cm; forest.checkConforming(cm);
        int nh; checkGeomConforming(mesh, edge, e2s, domain_xb, nh);
        cout << "  pass " << setw(2) << pass << ": +" << r << " refine  -> " << mesh.elem.rows()
             << " tris, nDof=" << nDof << "  hanging=" << nh << "\n";
        if (nh) { cout << "  ERROR: non-conforming mesh (hanging nodes)!\n"; return 2; }
    }

    // ----------------------- frame dirs -----------------------
    string meshDir = framesDir + "_mesh";
    std::vector<string> dirs; if (full_movie) dirs.push_back(framesDir); if (mesh_movie) dirs.push_back(meshDir);
    for (const string& d : dirs) {
        fs::create_directories(d);
        for (const auto& e : fs::directory_iterator(d)) if (e.path().extension() == ".ppm") fs::remove(e.path());
    }
    auto writeFrameSet = [&](int idx) {
        char fn[600];
        if (full_movie) { snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", framesDir.c_str(), idx); writeDensity(fn, rxa, rxb, rya, ryb, Wpix); }
        if (mesh_movie) {
            snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", meshDir.c_str(), idx);
            writeDensity(fn, rxa, rxb, rya, ryb, Wpix); overlayMesh(fn, mesh, rxa, rxb, rya, ryb);
        }
    };

    // ----------------------- time loop (variable dt; remesh every R steps) -----------------------
    double dt = dtFromCFL();
    cout << "\nEvolving (ARS(2,2,2) IMEX, HLLC + implicit AV) on the adaptive mesh...\n";
    cout << "  start: " << mesh.elem.rows() << " tris, nDof=" << nDof << ", h_min=" << scientific << setprecision(3)
         << hMinCFL() << ", dt=" << dt << fixed << "\n";
    auto t0 = chrono::high_resolution_clock::now();
    double t = 0.0; int step = 0, frame = 0;
    double frame_dt = t_end / std::max(1, n_frames), next_frame = 0.0;
    double maxRemeshDrift = 0.0;     // max relative change of the conserved totals across a remesh
    int lastRef = 0, lastCrs = 0;
    writeFrameSet(frame++); next_frame += frame_dt;

    while (t < t_end - 1e-12) {
        double dt_step = std::min(dt, t_end - t);
        auto aStep = clk();
        if (!dg->step(dt_step, t + dt_step, bc)) { cout << "ERROR: step " << step << " failed\n"; return 1; }
        if (PROFILE) Tstep += secs(aStep, clk());
        t += dt_step; ++step;

        // ---- adaptive remesh ----
        if (step % remesh_every == 0 && t < t_end - 1e-12) {
            auto aRem = clk();
            // conservation diagnostic is OPT-IN (check_conservation): it costs two
            // full-domain integral sweeps per remesh; the transfer is already proven
            // conservative (drift ~1e-13), so production skips it.
            Vector4d totBefore;
            if (check_conservation) { auto a = clk(); totBefore = dg->conservedTotals(); if (PROFILE) Tcons += secs(a, clk()); }
            forest.syncFromState(dg->state(), e2d, locDof);
            auto aF = clk(); std::vector<int> flag = computeFlags(/*allowCoarsen=*/true); if (PROFILE) Tflags += secs(aF, clk());
            dg.reset();
            auto aA = clk(); auto rc = forest.adapt(flag, max_gen); if (PROFILE) Tadapt += secs(aA, clk());
            lastRef = rc.first; lastCrs = rc.second;
            rebuildSolver();
            auto aG = clk(); dg->setState(forest.gatherState(e2d, locDof, nDof)); if (PROFILE) Tgs += secs(aG, clk());
            if (check_conservation) {
                auto a = clk(); Vector4d totAfter = dg->conservedTotals(); if (PROFILE) Tcons += secs(a, clk());
                maxRemeshDrift = std::max(maxRemeshDrift,
                                          (totAfter - totBefore).norm() / std::max(totBefore.norm(), 1e-30));
            }
            dt = dtFromCFL();                                  // CFL control: recompute dt from h_min + actual |u|+c
            if (PROFILE) Tremesh += secs(aRem, clk());
        }

        // ---- frame output (time-based) ----
        if (t >= next_frame - 1e-12 || t >= t_end - 1e-12) {
            auto aR = clk();
            writeFrameSet(frame++); next_frame += frame_dt;
            if (PROFILE) Trender += secs(aR, clk());
            double lam = dg->maxWaveSpeedGlobal();
            Vector4d tot = dg->conservedTotals();
            VectorXd rho = dg->densityField();
            double el = chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count();
            int gmax = 0; for (int tt = 0; tt < mesh.elem.rows(); ++tt) gmax = std::max(gmax, forest.gen(tt));
            cout << "  t=" << fixed << setprecision(4) << t << "  step=" << setw(6) << step
                 << "  tris=" << setw(7) << mesh.elem.rows()
                 << "  rho[" << setprecision(2) << rho.minCoeff() << "," << rho.maxCoeff() << "]"
                 << "  CFL=" << dt*lam*(2*ord+1)/hMinCFL()
                 << "  gmax=" << gmax << "  remesh(+" << lastRef << ",-" << lastCrs << ")"
                 << "  mass=" << setprecision(3) << tot(0)
                 << "  (" << setprecision(1) << el << "s)\n";
        }
    }

    // ----------------------- final stills + optional state -----------------------
    renderStills();
    if (save_state) {
        std::ofstream out("dmr_amr_state.bin", std::ios::binary);
        const MatrixXd& U = dg->state(); int r = (int)U.rows(), c = (int)U.cols();
        out.write((char*)&r, 4); out.write((char*)&c, 4);
        out.write((const char*)U.data(), (std::streamsize)r * c * sizeof(double));
        // also dump the mesh so the stills can be regenerated
        std::ofstream mo("dmr_amr_mesh.bin", std::ios::binary);
        int nn = (int)mesh.node.rows(), ne = (int)mesh.elem.rows();
        mo.write((char*)&nn, 4); mo.write((char*)&ne, 4);
        mo.write((const char*)mesh.node.data(), (std::streamsize)nn*2*sizeof(double));
        mo.write((const char*)mesh.elem.data(), (std::streamsize)ne*3*sizeof(int));
    }
    double wall = chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count();
    cout << "\nDone. " << frame << " frames, " << step << " steps.  final tris=" << mesh.elem.rows()
         << "  wall=" << fixed << setprecision(1) << wall << "s\n";
    if (check_conservation)
        cout << "  max conserved-totals drift across any remesh transfer = "
             << scientific << setprecision(2) << maxRemeshDrift << fixed
             << "  (machine-zero => the AMR transfer is conservative)\n";
    if (PROFILE) {
        double acc = Tstep + Tremesh + Trender;
        cout << fixed << setprecision(1)
             << "  [profile] step=" << Tstep << "s (" << 100*Tstep/acc << "%)  remesh=" << Tremesh
             << "s (" << 100*Tremesh/acc << "%)  render=" << Trender << "s (" << 100*Trender/acc << "%)\n"
             << "  [profile] remesh: rebuild[buildMesh=" << Tbuild << " FEM=" << Tfem << " getDOF=" << Tdof
             << " edge+tags=" << Tedge << " EulerDG=" << Tdgctor << "]  flags=" << Tflags
             << " adapt=" << Tadapt << " gather+set=" << Tgs << " consDiag=" << Tcons << " (s)\n";
    }
    cout << "  stills: dmr_amr_density.ppm, dmr_amr_density_zoom.ppm, dmr_amr_schlieren.ppm,\n"
         << "          dmr_amr_schlieren_zoom.ppm, dmr_amr_mesh.ppm\n";
    if (full_movie) cout << "  ffmpeg -y -framerate 25 -i " << framesDir << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr.mp4\n";
    if (mesh_movie) cout << "  ffmpeg -y -framerate 25 -i " << meshDir   << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr_mesh.mp4\n";
    return 0;
}
