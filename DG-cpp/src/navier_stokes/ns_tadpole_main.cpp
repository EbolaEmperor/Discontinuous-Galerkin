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
#include "Tadpole.h"
#include "Json.h"

using namespace Eigen;
using namespace std;
using namespace ns;
namespace fs = std::filesystem;

// ===========================================================================
// Cylinder + a passive-drift "tadpole" (rigid head + rigid tail).  Starts
// out in the far field; once it enters the wake of the cylinder it is
// drawn toward the low-speed refuge zone (Karman-gait style).  Bounces
// elastically off the cylinder wall and the four outer walls.
//
// Three rigid-body DOFs: head centre (x, y) and orientation theta.  The
// flow is the same DG NS solver as the bare-cylinder demo; the tadpole
// does NOT push back on the fluid (one-way coupling) -- we treat it as a
// small floating object whose feedback on the wake is negligible.
// ===========================================================================

namespace {

// Stamp a filled disc (anti-aliased) into an existing PPM file.  Used to
// draw the tadpole's head on top of every rendered frame.
void overlayDiscOnPPM(const string& path, double xc, double yc, double rad,
                      double xmin, double xmax, double ymin, double ymax,
                      int rgb)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return;
    string magic; in >> magic;
    if (magic != "P6") return;
    auto skip = [&](){
        int c = in.peek();
        while (c != EOF) {
            if (std::isspace(c)) { in.get(); c = in.peek(); }
            else if (c == '#') { string d; std::getline(in, d); c = in.peek(); }
            else break;
        }
    };
    skip();
    int W = 0, H = 0, maxv = 0;
    in >> W; skip(); in >> H; skip(); in >> maxv;
    in.get();
    std::vector<unsigned char> img((size_t)W * H * 3);
    in.read(reinterpret_cast<char*>(img.data()), (std::streamsize)img.size());
    in.close();
    unsigned char R = (unsigned char)((rgb >> 16) & 0xFF);
    unsigned char G = (unsigned char)((rgb >> 8) & 0xFF);
    unsigned char B = (unsigned char)( rgb        & 0xFF);
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    double col_c = (xc - xmin) / dx - 0.5;
    double row_c = (ymax - yc) / dy - 0.5;
    double rad_pix_x = rad / dx, rad_pix_y = rad / dy;
    double rad_pix = 0.5 * (rad_pix_x + rad_pix_y);
    int cmin = std::max(0, (int)std::floor(col_c - rad_pix - 1));
    int cmax = std::min(W - 1, (int)std::ceil(col_c + rad_pix + 1));
    int rmin = std::max(0, (int)std::floor(row_c - rad_pix - 1));
    int rmax = std::min(H - 1, (int)std::ceil(row_c + rad_pix + 1));
    for (int row = rmin; row <= rmax; ++row)
        for (int col = cmin; col <= cmax; ++col) {
            double dxp = col - col_c, dyp = row - row_c;
            double dist = std::hypot(dxp, dyp);
            double a = std::min(1.0, std::max(0.0, rad_pix + 0.5 - dist));
            if (a <= 0) continue;
            size_t idx = ((size_t)row * W + col) * 3;
            img[idx]     = (unsigned char)((1 - a) * img[idx]     + a * R);
            img[idx + 1] = (unsigned char)((1 - a) * img[idx + 1] + a * G);
            img[idx + 2] = (unsigned char)((1 - a) * img[idx + 2] + a * B);
        }
    std::ofstream out(path, std::ios::binary);
    if (!out) return;
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

} // anon

int main(int argc, char** argv) {
    // ----------------------- flow params (copied from ns_main defaults) ------
    int    ord    = 3;
    double radius = 0.5;
    double cx_ = 0.0, cy_ = 0.0;
    double xa = -5.0, xb = 18.0;
    double ya = -6.0, yb = 6.0;
    double h_      = 0.14;
    double farRatio = 9.0;
    double grade  = 0.08;
    double Re     = 100.0;
    double Uinf   = 1.0;
    double cfl    = 0.40;
    double dt     = 0.0;
    double t_end  = 250.0;
    int    save_every = 0;
    int    n_frames = 1000;
    double sigmaFac = 8.0;
    double gradDiv = 0.0;
    double ppeDivDamping = 30.0;
    string pressureModeName = "direct_ppe";
    double rampTime = 0.5;
    double perturb = 0.30;
    // rendering
    int    Wpix = 1000;
    double rxa = -3.0, rxb = 16.0;
    double rya = -5.0, ryb = 5.0;
    double vortClip = 3.0;
    std::vector<string> outputs;
    string framesDirVort  = "ns_tadpole_frames";
    string framesDirFlow  = "ns_tadpole_flow_frames";
    int    nParticles     = 1500;
    unsigned int particleSeed = 12345u;
    int    trailLen       = 28;
    int    trailStride    = 2;
    double bgDim          = 0.12;

    // ----------------------- tadpole params ----------------------------------
    double tad_x0 = 10.0, tad_y0 = 1.0;     // far-field start
    double tad_theta0 = 0.0;                // initially heading +x (away from cylinder)
    double tad_rHead = 0.15;
    double tad_Ltail = 0.6;
    double tad_rho   = 1.0;
    double tad_cDrag = 6.0;
    double tad_kRefuge = 6.0;
    double tad_swimForce = 0.0;
    double tad_dampLin = 0.4;
    double tad_dampAng = 1.0;
    double tad_maxSpeed = 1.5;
    int    tad_nSampleHead = 8;
    int    tad_nSampleTail = 12;
    // collision
    double restitution = 0.6;
    double tang_friction = 0.4;

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("ns_tadpole_config.json")) cfgPath = "ns_tadpole_config.json";
    if (!cfgPath.empty()) {
        string text;
        if (!cfgjson::readFile(cfgPath, text)) { cout << "ERROR: cannot open " << cfgPath << "\n"; return 1; }
        cfgjson::Json cfg;
        try { cfg = cfgjson::parse(text); }
        catch (const std::exception& e) { cout << "ERROR parsing " << cfgPath << ": " << e.what() << "\n"; return 1; }
        ord = cfg.getInt("ord", ord);
        radius = cfg.getNumber("radius", radius);
        cx_ = cfg.getNumber("cx", cx_); cy_ = cfg.getNumber("cy", cy_);
        xa = cfg.getNumber("xa", xa); xb = cfg.getNumber("xb", xb);
        ya = cfg.getNumber("ya", ya); yb = cfg.getNumber("yb", yb);
        h_ = cfg.getNumber("h", h_);
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
        // tadpole
        tad_x0   = cfg.getNumber("tad_x0", tad_x0);
        tad_y0   = cfg.getNumber("tad_y0", tad_y0);
        tad_theta0 = cfg.getNumber("tad_theta0", tad_theta0);
        tad_rHead = cfg.getNumber("tad_rHead", tad_rHead);
        tad_Ltail = cfg.getNumber("tad_Ltail", tad_Ltail);
        tad_rho   = cfg.getNumber("tad_rho", tad_rho);
        tad_cDrag = cfg.getNumber("tad_cDrag", tad_cDrag);
        tad_kRefuge = cfg.getNumber("tad_kRefuge", tad_kRefuge);
        tad_swimForce = cfg.getNumber("tad_swimForce", tad_swimForce);
        tad_dampLin = cfg.getNumber("tad_dampLin", tad_dampLin);
        tad_dampAng = cfg.getNumber("tad_dampAng", tad_dampAng);
        tad_maxSpeed = cfg.getNumber("tad_maxSpeed", tad_maxSpeed);
        tad_nSampleHead = cfg.getInt("tad_nSampleHead", tad_nSampleHead);
        tad_nSampleTail = cfg.getInt("tad_nSampleTail", tad_nSampleTail);
        restitution = cfg.getNumber("restitution", restitution);
        tang_friction = cfg.getNumber("tang_friction", tang_friction);
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
    if (dt <= 0.0) dt = cfl * h_ / ((2.0 * ord + 1.0) * Uinf);
    int nsteps = std::max(1, (int)std::lround(t_end / dt));
    if (save_every <= 0) save_every = std::max(1, nsteps / std::max(1, n_frames));
    double sigma = sigmaFac * (ord + 1) * (ord + 1);
    int pressureMode = (pressureModeName == "projection") ? NSPRESSURE_PROJECTION : NSPRESSURE_DIRECT_PPE;

    cout << "DG NS + tadpole rigid-body drift demo\n";
    cout << "  config:  " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  domain:  [" << xa << "," << xb << "] x [" << ya << "," << yb << "]"
         << "  cylinder r=" << radius << " at (" << cx_ << "," << cy_ << ")\n";
    cout << "  Re=" << Re << "  Uinf=" << Uinf << "  D=" << D << "  nu=" << nu
         << "  dP" << ord << "  h=" << h_ << "\n";
    cout << "  dt=" << dt << "  t_end=" << t_end << "  steps=" << nsteps
         << "  save_every=" << save_every << "\n";
    cout << "  tadpole: head r=" << tad_rHead << " tail L=" << tad_Ltail
         << "  start (x,y,theta)=(" << tad_x0 << "," << tad_y0 << "," << tad_theta0 << ")"
         << "  c_drag=" << tad_cDrag << " k_refuge=" << tad_kRefuge << "\n";

    // ----------------------- mesh + DG space ---------------------------------
    CylinderGeom geom{xa, xb, ya, yb, cx_, cy_, radius};
    Mesh mesh;
    generateCylinderMesh(mesh, geom, h_, farRatio, grade);
    FEM fem(ord, mesh);
    MatrixXi elem2dof; int nDof; fem.getDOF(mesh, elem2dof, nDof);
    MatrixXi edge, edge2side; mesh.getEdge2Side(edge, edge2side);
    VectorXi tag = classifyEdges(mesh, edge, edge2side, geom);
    cout << "  nDof=" << nDof << "\n";

    // ----------------------- BCs (same as the bare cylinder demo) -----------
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
    double rcyl = radius, ccx = cx_, ccy = cy_, htol = 0.3 * h_, Uin = Uinf, tr = rampTime;
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

    // ----------------------- IC ----------------------------------------------
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
    {
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
    }
    integ.setInitial(u0, v0);

    // ----------------------- tadpole + IB locator ----------------------------
    Tadpole tad;
    tad.x = tad_x0; tad.y = tad_y0; tad.theta = tad_theta0;
    tad.vx = 0.0; tad.vy = 0.0; tad.omega = 0.0;
    tad.rHead = tad_rHead; tad.Ltail = tad_Ltail; tad.rho = tad_rho;
    tad.cDrag = tad_cDrag; tad.kRefuge = tad_kRefuge;
    tad.swimForce = tad_swimForce;
    tad.dampLin = tad_dampLin; tad.dampAng = tad_dampAng;
    tad.maxSpeed = tad_maxSpeed;
    tad.nSampleHead = tad_nSampleHead; tad.nSampleTail = tad_nSampleTail;
    tad.initInertia();

    MeshLocator locatorIB; locatorIB.build(mesh);
    std::vector<int> ibHint;

    // ----------------------- output channels ---------------------------------
    int Hpix = std::max(2, (int)std::lround(Wpix * (ryb - rya) / (rxb - rxa)));
    if (Hpix % 2) ++Hpix; if (Wpix % 2) ++Wpix;
    auto inDomain = [=](double x, double y) {
        return std::hypot(x - cx_, y - cy_) > radius * 1.001 &&
               x > xa - 1e-9 && x < xb + 1e-9 && y > ya - 1e-9 && y < yb + 1e-9;
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
        if (kind == "vorticity") { ch.dir = framesDirVort; ch.mp4 = "tadpole_vortex.mp4"; }
        else if (kind == "speed"){ ch.dir = framesDirVort; ch.mp4 = "tadpole_speed.mp4"; }
        else /* flow */          {
            ch.dir = framesDirFlow + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx)));
            ch.mp4 = "tadpole_flow"   + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx))) + ".mp4";
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
        // Build a 2-point polyline for the tail (head centre -> tail tip).
        Vector2d hc = tad.headCentre();
        Vector2d tt = tad.tailTip();
        MatrixXd tailLine(2, 2);
        tailLine(0, 0) = hc.x(); tailLine(0, 1) = hc.y();
        tailLine(1, 0) = tt.x(); tailLine(1, 1) = tt.y();
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
            // Tail (white polyline) + head (white disc).
            int lineRGB = (ch.kind == "vorticity") ? 0x101010 : 0xFFFFFF;
            int discRGB = (ch.kind == "vorticity") ? 0x101010 : 0xFFFFFF;
            overlayPolylineOnPPM(string(fn), tailLine, rxa, rxb, rya, ryb, 2, lineRGB);
            overlayDiscOnPPM(string(fn), tad.x, tad.y, tad.rHead, rxa, rxb, rya, ryb, discRGB);
        }
    };
    auto advanceAllTracers = [&](double dtStep) {
        for (auto& ch : chans)
            if (ch.tracer)
                ch.tracer->advance(fem, mesh, elem2dof, integ.u(), integ.v(), locatorTracer, dtStep);
    };

    // ----------------------- collision helpers -------------------------------
    // Reflect a body point's velocity against an outward-pointing normal and
    // adjust the head-centre velocity accordingly.  The point may be the
    // head centre (offset (0,0)) or the tail tip (offset r2body).
    auto resolveCollision = [&](double rxC, double ryC, double nx, double ny,
                                double& cx_pt, double& cy_pt) {
        // Push the body point out of the obstacle.
        double cdx = cx_pt - rxC;
        double cdy = cy_pt - ryC;
        cx_pt = rxC + nx * std::abs(cdx * nx + cdy * ny + 0.0);
        cy_pt = ryC + ny * std::abs(cdx * nx + cdy * ny + 0.0);
        // Reflect velocity.
        double vn = tad.vx * nx + tad.vy * ny;
        if (vn < 0.0) {
            tad.vx -= (1.0 + restitution) * vn * nx;
            tad.vy -= (1.0 + restitution) * vn * ny;
            // Bleed tangential component.
            double vt_x = tad.vx - (tad.vx * nx + tad.vy * ny) * nx;
            double vt_y = tad.vy - (tad.vx * nx + tad.vy * ny) * ny;
            tad.vx -= tang_friction * vt_x;
            tad.vy -= tang_friction * vt_y;
        }
        tad.omega *= 0.7;
    };
    auto enforceCollisions = [&]() {
        // 1) Cylinder: head + tail tip kept outside r = radius + r_head.
        double rmin_h = radius + tad.rHead + 1e-3;
        double dx = tad.x - cx_, dy = tad.y - cy_;
        double dist = std::hypot(dx, dy);
        if (dist < rmin_h) {
            double nx = (dist > 1e-12) ? dx / dist : 1.0;
            double ny = (dist > 1e-12) ? dy / dist : 0.0;
            tad.x = cx_ + rmin_h * nx;
            tad.y = cy_ + rmin_h * ny;
            double vn = tad.vx * nx + tad.vy * ny;
            if (vn < 0.0) {
                tad.vx -= (1.0 + restitution) * vn * nx;
                tad.vy -= (1.0 + restitution) * vn * ny;
                double vt_x = tad.vx - (tad.vx * nx + tad.vy * ny) * nx;
                double vt_y = tad.vy - (tad.vx * nx + tad.vy * ny) * ny;
                tad.vx -= tang_friction * vt_x;
                tad.vy -= tang_friction * vt_y;
            }
            tad.omega *= 0.7;
        }
        // 2) Outer rectangle walls.  We treat each wall like a half-space.
        auto wallReflect = [&](double pos, double low, double hi, double nx, double ny) {
            if (pos < low) {
                double penetration = low - pos;
                tad.x += penetration * nx;
                tad.y += penetration * ny;
                double vn = tad.vx * nx + tad.vy * ny;
                if (vn < 0.0) {
                    tad.vx -= (1.0 + restitution) * vn * nx;
                    tad.vy -= (1.0 + restitution) * vn * ny;
                }
            } else if (pos > hi) {
                double penetration = pos - hi;
                tad.x -= penetration * nx;
                tad.y -= penetration * ny;
                double vn = tad.vx * nx + tad.vy * ny;
                if (vn > 0.0) {
                    tad.vx -= (1.0 + restitution) * vn * nx;
                    tad.vy -= (1.0 + restitution) * vn * ny;
                }
            }
        };
        // Head margin from each wall is rHead.
        wallReflect(tad.x, xa + tad.rHead, xb - tad.rHead, 1.0, 0.0);
        wallReflect(tad.y, ya + tad.rHead, yb - tad.rHead, 0.0, 1.0);

        // 3) Tail tip vs cylinder + walls.  We resolve by nudging the head
        // centre opposite to the tail offset (a coarse but stable approach).
        Vector2d tt = tad.tailTip();
        double tdx = tt.x() - cx_, tdy = tt.y() - cy_;
        double tdist = std::hypot(tdx, tdy);
        if (tdist < radius + 1e-3) {
            double nx = (tdist > 1e-12) ? tdx / tdist : 1.0;
            double ny = (tdist > 1e-12) ? tdy / tdist : 0.0;
            // Push the head opposite the tail's penetration direction.
            double pen = (radius + 1e-3) - tdist;
            tad.x += pen * nx;
            tad.y += pen * ny;
            tad.omega *= 0.5;
        }
        // Tail vs outer walls.
        if (tt.x() < xa) tad.x += (xa - tt.x());
        if (tt.x() > xb) tad.x -= (tt.x() - xb);
        if (tt.y() < ya) tad.y += (ya - tt.y());
        if (tt.y() > yb) tad.y -= (tt.y() - yb);
        (void)resolveCollision;
    };

    // ----------------------- diagnostics file --------------------------------
    ofstream hist("tadpole_diagnostics.csv");
    hist << "t,xtad,ytad,theta,vx,vy,omega,CD,CL\n";

    // ----------------------- time loop ---------------------------------------
    cout << "\nEvolving tadpole drift...\n";
    auto t0wall = chrono::high_resolution_clock::now();
    int frame = 0;
    writeAllFrames(frame); ++frame;
    const double qdyn = 0.5 * Uinf * Uinf * D;

    for (int s = 1; s <= nsteps; ++s) {
        // Advance the fluid (no body force; one-way coupling).
        if (!integ.step(s * dt)) {
            cout << "ERROR: NS solve failed at step " << s << "\n"; return 1;
        }
        // Advance the tadpole on the new fluid field.
        if (!tad.step(fem, mesh, elem2dof, integ.u(), integ.v(), locatorIB, ibHint, dt)) {
            cout << "WARN: tadpole entirely outside the fluid at step " << s << "\n";
        }
        // Resolve collisions with the cylinder + outer walls.
        enforceCollisions();
        // Tracers advance on the same fluid field.
        if (needLocator) advanceAllTracers(dt);

        // Diagnostics.
        double Fx, Fy; integ.cylinderForce(isCyl, Fx, Fy);
        double CD = Fx / qdyn, CL = Fy / qdyn;
        double t = s * dt;
        hist << t << "," << tad.x << "," << tad.y << "," << tad.theta << ","
             << tad.vx << "," << tad.vy << "," << tad.omega << ","
             << CD << "," << CL << "\n";

        if (s % save_every == 0 || s == nsteps) {
            writeAllFrames(frame); ++frame;
            double el = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
            cout << "  step " << setw(6) << s << "/" << nsteps
                 << "  t=" << fixed << setprecision(2) << t
                 << "  xy=(" << setprecision(2) << tad.x << "," << tad.y << ")"
                 << "  v=(" << setprecision(2) << tad.vx << "," << tad.vy << ")"
                 << "  CD=" << setprecision(3) << CD
                 << "  frame=" << setw(4) << (frame-1)
                 << "  (" << setprecision(1) << el << "s)\n";
        }
    }
    hist.close();

    double wall = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
    cout << "\nDone. " << frame << " frames in";
    for (size_t i = 0; i < chans.size(); ++i)
        cout << (i ? ", " : " ") << "'" << chans[i].dir << "/'";
    cout << ".  wall=" << fixed << setprecision(1) << wall << "s\n";
    cout << "  diagnostics -> tadpole_diagnostics.csv\n";
    cout << "\nAssemble the movies:\n";
    for (const auto& ch : chans)
        cout << "  ffmpeg -y -framerate 25 -i " << ch.dir << "/frame_%05d.ppm"
             << " -c:v libx264 -pix_fmt yuv420p -crf 18 " << ch.mp4 << "\n";
    return 0;
}
