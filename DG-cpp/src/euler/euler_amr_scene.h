#ifndef EULER_AMR_SCENE_H
#define EULER_AMR_SCENE_H

// ===========================================================================
// Reusable h-AMR scene driver for the DG + 2nd-order IMEX compressible-Euler
// solver.  This is the generic engine behind euler_dmr_amr, factored out so the
// shock-driven demos (Richtmyer-Meshkov, shock-bubble, corner diffraction, ...)
// share ONE proven adaptive driver and differ only in:
//
//   * the base mesh (a rectangle from makeRectMesh, or an L-shape built from it),
//   * the initial PRIMITIVE field  primIC(x,y) -> (rho,u,v,p),
//   * the boundary ghost states    bc(...)  keyed on an integer edge tag,
//   * the edge tagger              tagEdge(mx,my) -> tag (>0 boundary, 0 interior),
//   * an optional render mask      inDomain(x,y) (for non-rectangular domains).
//
// Everything else -- newest-vertex-bisection refine/coarsen, conservative state
// transfer, the density-gradient + inter-cell-jump indicator, CFL-controlled dt,
// HLLC + Persson-Peraire artificial viscosity + Zhang-Shu positivity, frame /
// mesh-overlay movies and final stills -- is identical to euler_dmr_amr and is
// reused byte-for-byte.  (LTS and the DMR-specific self-tests are intentionally
// omitted; LTS does not pay off on a growing instability region.)
// ===========================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "AdaptiveMesh.h"
#include "EulerDG.h"
#include "FEM.h"
#include "Mesh.h"

namespace euler {

using namespace Eigen;

// ---- per-element geometry helpers ----------------------------------------
inline double scene_longestEdge(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0)), p1 = m.node.row(m.elem(t, 1)), p2 = m.node.row(m.elem(t, 2));
    return std::max({(p1 - p0).norm(), (p2 - p1).norm(), (p0 - p2).norm()});
}
inline double scene_triArea(const Mesh& m, int t) {
    Vector2d p0 = m.node.row(m.elem(t, 0)), p1 = m.node.row(m.elem(t, 1)), p2 = m.node.row(m.elem(t, 2));
    return 0.5 * std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) - (p2.x() - p0.x()) * (p1.y() - p0.y()));
}
// CFL length scale = altitude to the longest edge (robust to thin NVB children).
inline double scene_hCFL(const Mesh& m, int t) {
    return 2.0 * scene_triArea(m, t) / std::max(scene_longestEdge(m, t), 1e-300);
}

// ---- minimal PPM read/write + mesh-edge overlay --------------------------
inline bool scene_readPPM(const std::string& path, int& W, int& H, std::vector<unsigned char>& img) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string magic; in >> magic; if (magic != "P6") return false;
    int maxv; in >> W >> H >> maxv; in.get();
    img.resize((size_t)W * H * 3);
    in.read(reinterpret_cast<char*>(img.data()), (std::streamsize)img.size());
    return true;
}
inline void scene_writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img) {
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}
inline void scene_drawLine(std::vector<unsigned char>& img, int W, int H, int x0, int y0, int x1, int y1,
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
inline void scene_overlayMesh(const std::string& path, const Mesh& mesh,
                              double xa, double xb, double ya, double yb) {
    int W, H; std::vector<unsigned char> img;
    if (!scene_readPPM(path, W, H, img)) return;
    auto toPix = [&](const Vector2d& p, int& px, int& py) {
        px = (int)std::lround((p.x() - xa) / (xb - xa) * W);
        py = (int)std::lround((yb - p.y()) / (yb - ya) * H);
    };
    for (int t = 0; t < mesh.elem.rows(); ++t) {
        Vector2d c = (mesh.node.row(mesh.elem(t,0)) + mesh.node.row(mesh.elem(t,1)) + mesh.node.row(mesh.elem(t,2))).transpose() / 3.0;
        if (c.x() < xa - 0.02*(xb-xa) || c.x() > xb + 0.02*(xb-xa) ||
            c.y() < ya - 0.02*(yb-ya) || c.y() > yb + 0.02*(yb-ya)) continue;
        for (int k = 0; k < 3; ++k) {
            Vector2d A = mesh.node.row(mesh.elem(t, k)), B = mesh.node.row(mesh.elem(t, (k + 1) % 3));
            int ax, ay, bx, by; toPix(A, ax, ay); toPix(B, bx, by);
            scene_drawLine(img, W, H, ax, ay, bx, by, 80, 200, 255, 0.32);   // faint cyan grid
        }
    }
    scene_writePPM(path, W, H, img);
}

// Build an L-shaped (or any axis-aligned multi-rectangle) base mesh by deleting,
// from the full criss-cross makeRectMesh, every CELL (= pair of triangles) whose
// centroid falls outside `keep`, then compacting the referenced vertices.  Both
// triangles of a deleted cell go together, so the survivors stay a valid
// criss-cross mesh (every diagonal is still the longest edge of both its tris)
// -> the AdaptiveForest NVB labelling and conformity closure are unaffected.
inline void makeMaskedRectMesh(Mesh& mesh, double x0, double x1, double y0, double y1,
                               int nx, int ny, const std::function<bool(double, double)>& keep) {
    Mesh full; makeRectMesh(full, x0, x1, y0, y1, nx, ny);
    std::vector<int> keepTri; keepTri.reserve(full.elem.rows());
    for (int t = 0; t < full.elem.rows(); ++t) {
        Vector2d c = (full.node.row(full.elem(t,0)) + full.node.row(full.elem(t,1)) + full.node.row(full.elem(t,2))).transpose() / 3.0;
        if (keep(c.x(), c.y())) keepTri.push_back(t);
    }
    std::vector<int> remap(full.node.rows(), -1);
    int nn = 0;
    for (int t : keepTri) for (int k = 0; k < 3; ++k) { int v = full.elem(t, k); if (remap[v] < 0) remap[v] = nn++; }
    mesh.node.resize(nn, 2);
    for (int v = 0; v < full.node.rows(); ++v) if (remap[v] >= 0) mesh.node.row(remap[v]) = full.node.row(v);
    mesh.elem.resize((int)keepTri.size(), 3);
    for (int i = 0; i < (int)keepTri.size(); ++i)
        for (int k = 0; k < 3; ++k) mesh.elem(i, k) = remap[full.elem(keepTri[i], k)];
}

// ---------------------------------------------------------------------------
// Scene description: numeric knobs (typically parsed from JSON in the driver
// main) + the four scene-defining hooks.
// ---------------------------------------------------------------------------
struct AMRScene {
    // identity / output
    std::string name      = "scene";        // file prefix for stills
    std::string framesDir = "out/scene_frames"; // density frame directory
    std::string cmapName  = "inferno";
    std::string title     = "Euler AMR scene";

    // discretization / time
    int    ord         = 2;
    int    max_gen     = 5;
    double t_end       = 1.0;
    double cfl         = 0.5;
    double lambda_safe = 6.0;
    int    n_frames    = 240;
    // dt-robustness: density floor used ONLY in the CFL wave-speed estimate (not the
    // physics).  0 = off (use dg->maxWaveSpeedGlobal()).  When >0, a near-vacuum nodal
    // value can no longer spike |u|+c and collapse dt -> the long shock-bubble run
    // survives.  Genuine cells (rho >= floor) are estimated exactly, so no CFL risk.
    double cfl_rho_floor = 0.0;

    // AMR controls
    int    remesh_every  = 6;
    int    buffer_layers = 3;
    int    init_passes   = 12;
    double th_ref = 0.30, th_crs = 0.08;

    // shock capturing
    bool   use_hllc = true;                      // HLLC (sharp contacts) vs Rusanov (more dissipative)
    double av_c = 1.2, av_kappa = 1.0, sigma_ip = 20.0;
    int    av_refresh = 5, av_indicator = 1;     // 1 = pressure sensor (keeps contacts sharp), 0 = density
    double av_s0 = -3.0; bool av_s0_set = false; // Persson sensor threshold override (lower => AV fires sooner)
    bool   use_av = true, use_pos = true;
    bool   use_tvb = false; double tvb_M = 0.0;  // minmod/TVB moment limiter (kills high-order overshoots)
    double rho_floor = 1e-12, p_floor = 1e-12;   // positivity floors

    // rendering
    int    Wpix = 1600;
    double rxa = 0, rxb = 1, rya = 0, ryb = 1;   // main render window
    double rho_lo = 0, rho_hi = 1;               // density colour range
    bool   zoom = false; double zxa = 0, zxb = 1, zya = 0, zyb = 1;
    bool   full_movie = true, mesh_movie = true, schlieren_still = true;

    // hooks (REQUIRED)
    std::function<Vector4d(double, double)> primIC;     // primitive IC (rho,u,v,p)
    ExteriorStateFn                         bc;         // ghost states keyed on tag
    std::function<int(double, double)>      tagEdge;    // (mx,my)->tag (>0 boundary, 0 interior/unknown)
    std::function<bool(double, double)>     inDomain;   // render mask (optional; empty = full window)
};

// ---------------------------------------------------------------------------
// Run a scene to t_end on an adaptive mesh built from `baseMesh`.  Returns 0 on
// success, nonzero on a failed step / non-conforming mesh.
// ---------------------------------------------------------------------------
inline int runAMRScene(const Mesh& baseMesh, const AMRScene& S) {
    namespace fs = std::filesystem;
    using std::cout; using std::string;
    const int ord = S.ord, locDof = (ord + 1) * (ord + 2) / 2, max_gen = S.max_gen;
    Colormap cmap = (S.cmapName == "jet") ? CM_JET : (S.cmapName == "viridis") ? CM_VIRIDIS
                  : (S.cmapName == "gray") ? CM_GRAY : (S.cmapName == "coolwarm") ? CM_COOLWARM : CM_INFERNO;

    cout << S.title << " -- DG dP" << ord << " + ARS(2,2,2) IMEX, HLLC + implicit AV, h-AMR (NVB)\n";
    cout << "  base mesh: " << baseMesh.elem.rows() << " tris, max_gen=" << max_gen
         << "  AV=" << (S.use_av ? "on" : "off") << " (c=" << S.av_c << ", sigma=" << S.sigma_ip
         << ", sensor=" << (S.av_indicator ? "p" : "rho") << ")  positivity=" << (S.use_pos ? "on" : "off")
         << "  t_end=" << S.t_end << "\n";
    cout << "  AMR: th_ref=" << S.th_ref << " th_crs=" << S.th_crs << " buffer=" << S.buffer_layers
         << " remesh_every=" << S.remesh_every << " init_passes=" << S.init_passes << "\n";

    EulerConfig dgcfg;
    dgcfg.use_hllc = S.use_hllc;
    dgcfg.use_av = S.use_av; dgcfg.use_positivity = S.use_pos;
    dgcfg.av_c = S.av_c; dgcfg.av_kappa = S.av_kappa; dgcfg.sigma_ip = S.sigma_ip;
    dgcfg.av_refresh = S.av_refresh; dgcfg.av_indicator = S.av_indicator;
    dgcfg.av_s0 = S.av_s0; dgcfg.av_s0_set = S.av_s0_set;
    dgcfg.use_tvb = S.use_tvb; dgcfg.tvb_M = S.tvb_M;
    dgcfg.rho_floor = S.rho_floor; dgcfg.p_floor = S.p_floor;

    AdaptiveForest forest(baseMesh, ord);

    // driver-owned persistent objects (EulerDG references these; must outlive it)
    Mesh mesh;
    std::unique_ptr<FEM> fem;
    MatrixXi e2d, edge, e2s; int nDof = 0; VectorXi tag;
    std::unique_ptr<EulerDG> dg;

    auto rebuildSolver = [&]() {
        forest.buildMesh(mesh);
        fem = std::make_unique<FEM>(ord, mesh, /*withHessian=*/false);
        fem->getDOF(mesh, e2d, nDof);
        mesh.getEdge2Side(edge, e2s);
        int NE = static_cast<int>(edge.rows());
        tag = VectorXi::Zero(NE);
        for (int e = 0; e < NE; ++e) {
            if (e2s(e, 0) != -1 && e2s(e, 1) != -1) continue;   // interior edge
            Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1))).transpose();
            tag(e) = S.tagEdge(m.x(), m.y());
        }
        dg = std::make_unique<EulerDG>(*fem, mesh, e2d, edge, e2s, tag, dgcfg);
    };

    // strong geometric conformity test: every 1-sided edge must be a real boundary
    // (tagEdge > 0); a 1-sided interior edge (tag 0) is a hanging node.
    auto checkConformGeom = [&](int& nHanging) {
        nHanging = 0;
        for (int e = 0; e < edge.rows(); ++e) {
            if (e2s(e, 0) != -1 && e2s(e, 1) != -1) continue;
            Vector2d m = 0.5 * (mesh.node.row(edge(e, 0)) + mesh.node.row(edge(e, 1))).transpose();
            if (S.tagEdge(m.x(), m.y()) == 0) ++nHanging;
        }
    };

    // indicator: scaled density gradient (centroid) + inter-cell density jump
    auto computeFlags = [&](bool allowCoarsen) {
        int NT = static_cast<int>(mesh.elem.rows());
        const MatrixXd& U = dg->state();
        VectorXd rbar(NT), eta(NT);
        MatrixXd dphiC = fem->computeBasisDlam_all(Vector3d(1.0/3, 1.0/3, 1.0/3));   // 3 x locDof
        VectorXd rblk(locDof);
        for (int t = 0; t < NT; ++t) {
            double rb = 0;
            for (int i = 0; i < locDof; ++i) { rblk(i) = U(e2d(t, i), 0); rb += rblk(i); }
            rb /= locDof; rbar(t) = rb;
            Vector2d g = fem->Dlam[t] * (dphiC * rblk);
            eta(t) = scene_longestEdge(mesh, t) * g.norm() / std::max(rb, 1e-9);
        }
        for (int e = 0; e < e2s.rows(); ++e) {
            int a = e2s(e, 0), b = e2s(e, 1);
            if (a < 0 || b < 0) continue;
            double j = std::abs(rbar(a) - rbar(b)) / std::max(std::min(rbar(a), rbar(b)), 1e-9);
            eta(a) = std::max(eta(a), j); eta(b) = std::max(eta(b), j);
        }
        std::vector<int> flag(NT, 0);
        for (int t = 0; t < NT; ++t) {
            if (eta(t) > S.th_ref && forest.gen(t) < max_gen) flag[t] = 1;
            else if (allowCoarsen && eta(t) < S.th_crs && forest.gen(t) > 0) flag[t] = -1;
        }
        std::vector<std::vector<int>> nbr(NT);
        for (int e = 0; e < e2s.rows(); ++e) {
            int a = e2s(e, 0), b = e2s(e, 1);
            if (a >= 0 && b >= 0) { nbr[a].push_back(b); nbr[b].push_back(a); }
        }
        std::vector<int> front;
        for (int t = 0; t < NT; ++t) if (flag[t] == 1) front.push_back(t);
        for (int layer = 0; layer < S.buffer_layers; ++layer) {
            std::vector<int> next;
            for (int t : front) for (int n : nbr[t])
                if (flag[n] != 1 && forest.gen(n) < max_gen) { flag[n] = 1; next.push_back(n); }
            front.swap(next);
            if (front.empty()) break;
        }
        for (int t = 0; t < NT; ++t)
            if (flag[t] == -1)
                for (int n : nbr[t]) if (flag[n] == 1) { flag[t] = 0; break; }
        return flag;
    };

    auto hMinCFL = [&]() {
        double h = 1e300; for (int t = 0; t < mesh.elem.rows(); ++t) h = std::min(h, scene_hCFL(mesh, t));
        return h;
    };
    // robust global max wave speed |u|+c.  With cfl_rho_floor>0 the per-dof density is
    // floored ONLY here (for the dt estimate), so a spurious near-vacuum node can't
    // blow up the sound speed and collapse dt; genuine cells are unchanged.
    auto maxWaveRobust = [&]() {
        if (S.cfl_rho_floor <= 0.0) return dg->maxWaveSpeedGlobal();
        const MatrixXd& U = dg->state();
        const double rf = S.cfl_rho_floor, g = GAMMA;
        double lam = 0.0;
        for (int i = 0; i < U.rows(); ++i) {
            double rho = std::max(U(i, 0), rf);
            double mom2 = U(i, 1) * U(i, 1) + U(i, 2) * U(i, 2);
            double p = std::max((g - 1.0) * (U(i, 3) - 0.5 * mom2 / rho), 0.0);
            double sp = std::sqrt(mom2) / rho + std::sqrt(g * p / rho);
            lam = std::max(lam, sp);
        }
        return lam;
    };
    auto dtFromCFL = [&]() {
        double lam = std::max(S.lambda_safe, 1.15 * maxWaveRobust());
        return S.cfl * hMinCFL() / ((2.0 * ord + 1.0) * lam);
    };

    // rendering helpers
    int Wpix = S.Wpix; if (Wpix % 2) ++Wpix;
    auto writeDensity = [&](const string& path, double xa, double xb, double ya, double yb, int W) {
        int Hh = std::max(2, (int)std::lround(W * (yb - ya) / (xb - xa))); if (Hh % 2) ++Hh;
        writeScalarPPM(path, *fem, mesh, e2d, dg->densityField(), W, Hh, xa, xb, ya, yb, S.rho_lo, S.rho_hi, cmap, S.inDomain);
    };
    auto writeSchlieren = [&](const string& path, double xa, double xb, double ya, double yb, int W) {
        int Hh = std::max(2, (int)std::lround(W * (yb - ya) / (xb - xa))); if (Hh % 2) ++Hh;
        writeSchlierenPPM(path, *fem, mesh, e2d, dg->densityField(), W, Hh, xa, xb, ya, yb, 8.0);
    };

    // ----- initial condition + initial adaptation (re-project the sharp IC) -----
    rebuildSolver();
    dg->setState(projectInitial(*fem, mesh, e2d, S.primIC));
    cout << "\nInitial adaptation:\n";
    for (int pass = 0; pass < S.init_passes; ++pass) {
        std::vector<int> flag = computeFlags(/*allowCoarsen=*/false);
        int nref = 0; for (int f : flag) if (f == 1) ++nref;
        if (nref == 0) break;
        forest.syncFromState(dg->state(), e2d, locDof);
        dg.reset();
        auto [r, c] = forest.adapt(flag, max_gen); (void)c;
        rebuildSolver();
        dg->setState(projectInitial(*fem, mesh, e2d, S.primIC));
        int nh; checkConformGeom(nh);
        cout << "  pass " << std::setw(2) << pass << ": +" << r << " refine -> " << mesh.elem.rows()
             << " tris, nDof=" << nDof << "  hanging=" << nh << "\n";
        if (nh) { cout << "  ERROR: non-conforming mesh (hanging nodes)!\n"; return 2; }
    }

    // ----- frame directories -----
    string meshDir = S.framesDir + "_mesh";
    std::vector<string> dirs;
    if (S.full_movie) dirs.push_back(S.framesDir);
    if (S.mesh_movie) dirs.push_back(meshDir);
    for (const string& d : dirs) {
        fs::create_directories(d);
        for (const auto& e : fs::directory_iterator(d)) if (e.path().extension() == ".ppm") fs::remove(e.path());
    }
    auto writeFrameSet = [&](int idx) {
        char fn[600];
        if (S.full_movie) { snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", S.framesDir.c_str(), idx); writeDensity(fn, S.rxa, S.rxb, S.rya, S.ryb, Wpix); }
        if (S.mesh_movie) {
            snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", meshDir.c_str(), idx);
            writeDensity(fn, S.rxa, S.rxb, S.rya, S.ryb, Wpix); scene_overlayMesh(fn, mesh, S.rxa, S.rxb, S.rya, S.ryb);
        }
    };

    // ----- time loop (variable dt; remesh every R steps) -----
    double dt = dtFromCFL();
    cout << "\nEvolving on the adaptive mesh...\n";
    cout << "  start: " << mesh.elem.rows() << " tris, nDof=" << nDof << ", h_min="
         << std::scientific << std::setprecision(3) << hMinCFL() << ", dt=" << dt << std::fixed << "\n";
    auto t0 = std::chrono::high_resolution_clock::now();
    double t = 0.0; int step = 0, frame = 0;
    double frame_dt = S.t_end / std::max(1, S.n_frames), next_frame = 0.0;
    int lastRef = 0, lastCrs = 0; double effCfl = 0.0;
    writeFrameSet(frame++); next_frame += frame_dt;

    const int max_steps = 400000;   // backstop against a runaway run
    while (t < S.t_end - 1e-12) {
        // dt-collapse / NaN guard: a near-vacuum or degenerate cell can spike the wave
        // speed so dtFromCFL() returns a vanishing (or non-finite) dt; without this the
        // loop crawls forever (t advances by ~0) instead of finishing.  Stop cleanly and
        // keep the frames written so far.
        if (!std::isfinite(dt) || dt < (S.t_end > 0 ? S.t_end : 1.0) * 1e-9) {
            cout << "  [stop] dt collapsed to " << std::scientific << dt << std::fixed
                 << " at t=" << t << " (step " << step << ") -- near-vacuum stiffness; "
                 << "stopping with " << frame << " frames written.\n";
            break;
        }
        if (step >= max_steps) { cout << "  [stop] hit max_steps=" << max_steps << " at t=" << t << "\n"; break; }
        double dt_step = std::min(dt, S.t_end - t);
        if (!dg->step(dt_step, t + dt_step, S.bc)) { cout << "ERROR: step " << step << " failed\n"; return 1; }
        effCfl = dt_step * dg->maxWaveSpeedGlobal() * (2 * ord + 1) / hMinCFL();
        t += dt_step; ++step;

        if (step % S.remesh_every == 0 && t < S.t_end - 1e-12) {
            forest.syncFromState(dg->state(), e2d, locDof);
            std::vector<int> flag = computeFlags(/*allowCoarsen=*/true);
            dg.reset();
            auto rc = forest.adapt(flag, max_gen);
            lastRef = rc.first; lastCrs = rc.second;
            rebuildSolver();
            dg->setState(forest.gatherState(e2d, locDof, nDof));
            dt = dtFromCFL();
        }

        if (t >= next_frame - 1e-12 || t >= S.t_end - 1e-12) {
            writeFrameSet(frame++); next_frame += frame_dt;
            Vector4d tot = dg->conservedTotals();
            VectorXd rho = dg->densityField();
            int gmax = 0; for (int tt = 0; tt < mesh.elem.rows(); ++tt) gmax = std::max(gmax, forest.gen(tt));
            double el = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
            cout << "  t=" << std::fixed << std::setprecision(4) << t << "  step=" << std::setw(6) << step
                 << "  tris=" << std::setw(7) << mesh.elem.rows()
                 << "  rho[" << std::setprecision(2) << rho.minCoeff() << "," << rho.maxCoeff() << "]"
                 << "  CFL=" << std::setprecision(2) << effCfl
                 << "  gmax=" << gmax << "  remesh(+" << lastRef << ",-" << lastCrs << ")"
                 << "  mass=" << std::setprecision(3) << tot(0)
                 << "  (" << std::setprecision(1) << el << "s)\n";
        }
    }

    // ----- final stills -----
    writeDensity(S.name + "_density.ppm", S.rxa, S.rxb, S.rya, S.ryb, Wpix);
    if (S.schlieren_still) writeSchlieren(S.name + "_schlieren.ppm", S.rxa, S.rxb, S.rya, S.ryb, Wpix);
    {
        string meshStill = S.name + "_mesh.ppm";
        writeDensity(meshStill, S.rxa, S.rxb, S.rya, S.ryb, Wpix);
        scene_overlayMesh(meshStill, mesh, S.rxa, S.rxb, S.rya, S.ryb);
    }
    if (S.zoom) {
        writeDensity(S.name + "_density_zoom.ppm", S.zxa, S.zxb, S.zya, S.zyb, Wpix);
        if (S.schlieren_still) writeSchlieren(S.name + "_schlieren_zoom.ppm", S.zxa, S.zxb, S.zya, S.zyb, Wpix);
    }
    double wall = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    cout << "\nDone. " << frame << " frames, " << step << " steps.  final tris=" << mesh.elem.rows()
         << "  wall=" << std::fixed << std::setprecision(1) << wall << "s\n";
    cout << "  stills: " << S.name << "_density.ppm";
    if (S.schlieren_still) cout << ", " << S.name << "_schlieren.ppm";
    cout << ", " << S.name << "_mesh.ppm\n";
    if (S.full_movie) cout << "  ffmpeg -y -framerate 25 -i " << S.framesDir << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 " << S.name << ".mp4\n";
    if (S.mesh_movie) cout << "  ffmpeg -y -framerate 25 -i " << meshDir << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 " << S.name << "_mesh.mp4\n";
    return 0;
}

} // namespace euler

#endif
