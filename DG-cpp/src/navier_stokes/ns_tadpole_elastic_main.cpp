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
#include "IBCoupler.h"
#include "ElasticTadpole.h"
#include "Json.h"

using namespace Eigen;
using namespace std;
using namespace ns;
namespace fs = std::filesystem;

// ===========================================================================
// Cylinder + a passive-drift tadpole with an immersed elastic tail.
//
// The head remains the small rigid body from the tadpole drift demo.  The tail
// is a Cosserat/Kirchhoff rod clamped to the moving head and coupled to the
// flow by FE immersed-boundary transfer: interpolate u_h to rod nodes, impose
// a relaxed no-slip penalty, and scatter the equal reaction into the DG RHS.
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
    string framesDirVort  = "ns_tadpole_elastic_frames";
    string framesDirFlow  = "ns_tadpole_elastic_flow_frames";
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
    double tad_rhoHead = 1.0;
    double tad_rhoTail = 1.0;
    int    tad_Ntail   = 24;
    double tad_KB      = 0.01;
    double tad_KS      = 1e4;
    double tad_tailDamp = 0.20;
    double tad_cDragHead = 6.0;
    double tad_cDragTail = 18.0;
    double tad_kRefuge = 0.0;
    double tad_swimForce = 0.0;
    double tad_dampLin = 0.4;
    double tad_dampAng = 10.0;
    double tad_maxSpeed = 1.5;
    double tad_maxOmega = 0.5;
    double tad_maxOmegaAccel = 1.0;
    int    tad_nSampleHead = 16;
    double tad_ibAlpha = 0.15;        // relaxed direct-forcing strength, F=(alpha/dt)(V-u)
    double tad_ibForceCap = 80.0;     // per-unit-length cap on IB force
    int    tad_ibSmooth = 1;          // 1-2-1 smoothing passes along the tail
    // collision
    double restitution = 0.6;
    double tang_friction = 0.4;

    string cfgPath;
    if (argc > 1) cfgPath = argv[1];
    else if (fs::exists("ns_tadpole_elastic_config.json")) cfgPath = "ns_tadpole_elastic_config.json";
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
        tad_rhoHead = cfg.getNumber("tad_rhoHead", tad_rhoHead);
        tad_rhoTail = cfg.getNumber("tad_rhoTail", tad_rhoTail);
        tad_Ntail   = cfg.getInt("tad_Ntail", tad_Ntail);
        tad_KB      = cfg.getNumber("tad_KB", tad_KB);
        tad_KS      = cfg.getNumber("tad_KS", tad_KS);
        tad_tailDamp = cfg.getNumber("tad_tailDamp", tad_tailDamp);
        tad_cDragHead = cfg.getNumber("tad_cDragHead", tad_cDragHead);
        tad_cDragTail = cfg.getNumber("tad_cDragTail", tad_cDragTail);
        tad_kRefuge = cfg.getNumber("tad_kRefuge", tad_kRefuge);
        tad_swimForce = cfg.getNumber("tad_swimForce", tad_swimForce);
        tad_dampLin = cfg.getNumber("tad_dampLin", tad_dampLin);
        tad_dampAng = cfg.getNumber("tad_dampAng", tad_dampAng);
        tad_maxSpeed = cfg.getNumber("tad_maxSpeed", tad_maxSpeed);
        tad_maxOmega = cfg.getNumber("tad_maxOmega", tad_maxOmega);
        tad_maxOmegaAccel = cfg.getNumber("tad_maxOmegaAccel", tad_maxOmegaAccel);
        tad_nSampleHead = cfg.getInt("tad_nSampleHead", tad_nSampleHead);
        tad_ibAlpha = cfg.getNumber("tad_ibAlpha", tad_ibAlpha);
        tad_ibForceCap = cfg.getNumber("tad_ibForceCap", tad_ibForceCap);
        tad_ibSmooth = cfg.getInt("tad_ibSmooth", tad_ibSmooth);
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

    cout << "DG NS + elastic-tail tadpole immersed-boundary demo\n";
    cout << "  config:  " << (cfgPath.empty() ? string("(built-in defaults)") : cfgPath) << "\n";
    cout << "  domain:  [" << xa << "," << xb << "] x [" << ya << "," << yb << "]"
         << "  cylinder r=" << radius << " at (" << cx_ << "," << cy_ << ")\n";
    cout << "  Re=" << Re << "  Uinf=" << Uinf << "  D=" << D << "  nu=" << nu
         << "  dP" << ord << "  h=" << h_ << "\n";
    cout << "  dt=" << dt << "  t_end=" << t_end << "  steps=" << nsteps
         << "  save_every=" << save_every << "\n";
    cout << "  tadpole: head r=" << tad_rHead << " tail L=" << tad_Ltail
         << "  start (x,y,theta)=(" << tad_x0 << "," << tad_y0 << "," << tad_theta0 << ")"
         << "  cDragH=" << tad_cDragHead << " cDragT=" << tad_cDragTail
         << " swim=" << tad_swimForce << " K_B=" << tad_KB
         << " ibAlpha=" << tad_ibAlpha << "\n";

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
    ElasticTadpole tad;
    tad.x = tad_x0; tad.y = tad_y0; tad.theta = tad_theta0;
    tad.vx = 0.0; tad.vy = 0.0; tad.omega = 0.0;
    tad.rHead = tad_rHead; tad.Ltail = tad_Ltail;
    tad.rhoHead = tad_rhoHead; tad.rhoTail = tad_rhoTail;
    tad.Ntail = tad_Ntail; tad.KB = tad_KB; tad.KS = tad_KS;
    tad.tailDamp = tad_tailDamp;
    tad.cDragHead = tad_cDragHead; tad.cDragTail = tad_cDragTail;
    tad.kRefuge = tad_kRefuge;
    tad.swimForce = tad_swimForce;
    tad.dampLin = tad_dampLin; tad.dampAng = tad_dampAng;
    tad.maxSpeed = tad_maxSpeed; tad.maxOmega = tad_maxOmega;
    tad.maxOmegaAccel = tad_maxOmegaAccel;
    tad.nSampleHead = tad_nSampleHead;
    tad.init();

    MeshLocator locatorIB; locatorIB.build(mesh);
    std::vector<int> headHint, tailHint, tailForceHint;
    std::vector<char> contactPrev(tad.tail ? tad.tail->N + 1 : 0, 0);

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
        if (kind == "vorticity") { ch.dir = framesDirVort; ch.mp4 = "tadpole_elastic_vortex.mp4"; }
        else if (kind == "speed"){ ch.dir = framesDirVort; ch.mp4 = "tadpole_elastic_speed.mp4"; }
        else /* flow */          {
            ch.dir = framesDirFlow + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx)));
            ch.mp4 = "tadpole_elastic_flow" + (flowIdx == 0 ? string("") : ("_" + std::to_string(flowIdx))) + ".mp4";
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
        // The tail polyline is the full Cosserat rod node trail.
        MatrixXd tailLine;
        if (tad.tail) tailLine = tad.tail->X;
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
            int lineRGB = (ch.kind == "vorticity") ? 0x101010 : 0xFFFFFF;
            int discRGB = (ch.kind == "vorticity") ? 0x101010 : 0xFFFFFF;
            if (tailLine.rows() >= 2)
                overlayPolylineOnPPM(string(fn), tailLine, rxa, rxb, rya, ryb, 2, lineRGB);
            overlayDiscOnPPM(string(fn), tad.x, tad.y, tad.rHead, rxa, rxb, rya, ryb, discRGB);
        }
    };
    auto advanceAllTracers = [&](double dtStep) {
        for (auto& ch : chans)
            if (ch.tracer)
                ch.tracer->advance(fem, mesh, elem2dof, integ.u(), integ.v(), locatorTracer, dtStep);
    };

    auto tailShapeMetrics = [&]() {
        struct Metrics { double rmsCurv, maxCurv, roughness; };
        Metrics m{0.0, 0.0, 0.0};
        if (!tad.tail || tad.tail->N < 2) return m;
        std::vector<double> curv;
        curv.reserve(std::max(0, tad.tail->N - 1));
        for (int i = 1; i < tad.tail->N; ++i) {
            Vector2d a(tad.tail->X(i, 0)     - tad.tail->X(i - 1, 0),
                       tad.tail->X(i, 1)     - tad.tail->X(i - 1, 1));
            Vector2d b(tad.tail->X(i + 1, 0) - tad.tail->X(i, 0),
                       tad.tail->X(i + 1, 1) - tad.tail->X(i, 1));
            double la = a.norm(), lb = b.norm();
            if (la < 1e-14 || lb < 1e-14) continue;
            double cr = a.x() * b.y() - a.y() * b.x();
            double dtp = a.x() * b.x() + a.y() * b.y();
            double th = std::atan2(cr, dtp);
            double kappa = th / (0.5 * (la + lb));
            curv.push_back(kappa);
            m.rmsCurv += kappa * kappa;
            m.maxCurv = std::max(m.maxCurv, std::abs(kappa));
        }
        if (!curv.empty()) m.rmsCurv = std::sqrt(m.rmsCurv / (double)curv.size());
        if (curv.size() >= 2) {
            for (size_t i = 1; i < curv.size(); ++i) {
                double d = curv[i] - curv[i - 1];
                m.roughness += d * d;
            }
            m.roughness = std::sqrt(m.roughness / (double)(curv.size() - 1));
        }
        return m;
    };

    VectorXd loadFx = VectorXd::Zero(nDof);
    VectorXd loadFy = VectorXd::Zero(nDof);
    auto buildElasticTailIBLoad = [&]() {
        loadFx.setZero();
        loadFy.setZero();
        if (!tad.tail || tad_ibAlpha <= 0.0) return;
        tad.syncTailClamp();
        auto sampleN = meshToRod(fem, mesh, elem2dof, integ.u(), integ.v(),
                                 locatorIB, *tad.tail, tailForceHint);
        VectorXd ds = rodArclengthWeights(*tad.tail);
        MatrixXd Fk(tad.tail->N + 1, 2);
        Fk.setZero();
        const double alpha = tad_ibAlpha / dt;
        const double cap = std::max(0.0, tad_ibForceCap);
        for (int k = 2; k <= tad.tail->N; ++k) {
            if (!sampleN.alive[k]) continue;
            if (k < (int)contactPrev.size() && contactPrev[k]) continue;
            double fx = alpha * (tad.tail->V(k, 0) - sampleN.uv(k, 0));
            double fy = alpha * (tad.tail->V(k, 1) - sampleN.uv(k, 1));
            double nrm = std::hypot(fx, fy);
            if (cap > 0.0 && nrm > cap) {
                double s = cap / nrm;
                fx *= s; fy *= s;
            }
            Fk(k, 0) = fx;
            Fk(k, 1) = fy;
        }
        int passes = std::max(0, tad_ibSmooth);
        for (int pass = 0; pass < passes; ++pass) {
            MatrixXd old = Fk;
            for (int k = 3; k < tad.tail->N; ++k) {
                if (k < (int)contactPrev.size() && contactPrev[k]) continue;
                Fk.row(k) = 0.25 * old.row(k - 1) + 0.5 * old.row(k) + 0.25 * old.row(k + 1);
            }
        }
        rodToMesh(fem, mesh, elem2dof, locatorIB, *tad.tail, Fk, ds,
                  tailForceHint, loadFx, loadFy);
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
        std::vector<char> nextContact(tad.tail ? tad.tail->N + 1 : 0, 0);
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

        // 3) Elastic tail nodes vs cylinder + outer walls.  Project each rod
        // node out of the obstacle and reflect its inward radial velocity, so
        // the tail bends and slaps off the cylinder rather than tunnelling
        // through it.  Skip the two clamped nodes (already pinned to the head).
        if (tad.tail) {
            const double r2_cyl = (radius + 1e-3) * (radius + 1e-3);
            for (int k = 2; k <= tad.tail->N; ++k) {
                // Cylinder.
                double tdx = tad.tail->X(k, 0) - cx_, tdy = tad.tail->X(k, 1) - cy_;
                double r2k = tdx * tdx + tdy * tdy;
                if (r2k < r2_cyl) {
                    if (k < (int)nextContact.size()) nextContact[k] = 1;
                    double rk = std::sqrt(std::max(r2k, 1e-300));
                    double nx = (rk > 0) ? tdx / rk : 1.0;
                    double ny = (rk > 0) ? tdy / rk : 0.0;
                    tad.tail->X(k, 0) = cx_ + (radius + 1e-3) * nx;
                    tad.tail->X(k, 1) = cy_ + (radius + 1e-3) * ny;
                    double vn = tad.tail->V(k, 0) * nx + tad.tail->V(k, 1) * ny;
                    if (vn < 0.0) {
                        double s = (1.0 + restitution) * vn;
                        tad.tail->V(k, 0) -= s * nx;
                        tad.tail->V(k, 1) -= s * ny;
                    }
                }
                // Outer rectangle walls.
                if (tad.tail->X(k, 0) < xa) {
                    if (k < (int)nextContact.size()) nextContact[k] = 1;
                    tad.tail->X(k, 0) = xa;
                    if (tad.tail->V(k, 0) < 0) tad.tail->V(k, 0) = -restitution * tad.tail->V(k, 0);
                } else if (tad.tail->X(k, 0) > xb) {
                    if (k < (int)nextContact.size()) nextContact[k] = 1;
                    tad.tail->X(k, 0) = xb;
                    if (tad.tail->V(k, 0) > 0) tad.tail->V(k, 0) = -restitution * tad.tail->V(k, 0);
                }
                if (tad.tail->X(k, 1) < ya) {
                    if (k < (int)nextContact.size()) nextContact[k] = 1;
                    tad.tail->X(k, 1) = ya;
                    if (tad.tail->V(k, 1) < 0) tad.tail->V(k, 1) = -restitution * tad.tail->V(k, 1);
                } else if (tad.tail->X(k, 1) > yb) {
                    if (k < (int)nextContact.size()) nextContact[k] = 1;
                    tad.tail->X(k, 1) = yb;
                    if (tad.tail->V(k, 1) > 0) tad.tail->V(k, 1) = -restitution * tad.tail->V(k, 1);
                }
            }
            tad.syncTailClamp();
        }
        contactPrev.swap(nextContact);
        (void)resolveCollision;
    };

    // ----------------------- diagnostics file --------------------------------
    ofstream hist("tadpole_elastic_diagnostics.csv");
    hist << "t,xtad,ytad,theta,vx,vy,omega,tailTipX,tailTipY,maxTailStrain,"
         << "tailRmsCurv,tailMaxCurv,tailRoughness,CD,CL\n";

    // ----------------------- time loop ---------------------------------------
    cout << "\nEvolving elastic-tail tadpole drift...\n";
    auto t0wall = chrono::high_resolution_clock::now();
    int frame = 0;
    writeAllFrames(frame); ++frame;
    const double qdyn = 0.5 * Uinf * Uinf * D;

    for (int s = 1; s <= nsteps; ++s) {
        // Advance the fluid with the tail's immersed-boundary reaction force.
        buildElasticTailIBLoad();
        if (!integ.stepWithBodyForce(s * dt, loadFx, loadFy)) {
            cout << "ERROR: NS solve failed at step " << s << "\n"; return 1;
        }
        // Advance the tadpole on the new fluid field.
        if (!tad.step(fem, mesh, elem2dof, integ.u(), integ.v(), locatorIB, headHint, tailHint, dt)) {
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
        Vector2d tt = tad.tailTip();
        double maxStrain = tad.tail ? tad.tail->maxEdgeStrain() : 0.0;
        auto tm = tailShapeMetrics();
        hist << t << "," << tad.x << "," << tad.y << "," << tad.theta << ","
             << tad.vx << "," << tad.vy << "," << tad.omega << ","
             << tt.x() << "," << tt.y() << "," << maxStrain << ","
             << tm.rmsCurv << "," << tm.maxCurv << "," << tm.roughness << ","
             << CD << "," << CL << "\n";

        if (s % save_every == 0 || s == nsteps) {
            writeAllFrames(frame); ++frame;
            double el = chrono::duration<double>(chrono::high_resolution_clock::now() - t0wall).count();
            cout << "  step " << setw(6) << s << "/" << nsteps
                 << "  t=" << fixed << setprecision(2) << t
                 << "  xy=(" << setprecision(2) << tad.x << "," << tad.y << ")"
                 << "  v=(" << setprecision(2) << tad.vx << "," << tad.vy << ")"
                 << "  strain=" << scientific << setprecision(2) << maxStrain << fixed
                 << "  rough=" << scientific << setprecision(2) << tm.roughness << fixed
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
    cout << "  diagnostics -> tadpole_elastic_diagnostics.csv\n";
    cout << "\nAssemble the movies:\n";
    for (const auto& ch : chans)
        cout << "  ffmpeg -y -framerate 25 -i " << ch.dir << "/frame_%05d.ppm"
             << " -c:v libx264 -pix_fmt yuv420p -crf 18 " << ch.mp4 << "\n";
    return 0;
}
