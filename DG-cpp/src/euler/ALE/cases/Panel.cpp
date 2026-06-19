#include "Panel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace euler_ale {

using namespace Eigen;

double PanelMap::mode(double x) const {
    return std::sin(M_PI * std::clamp(x, 0.0, 1.0));
}

double PanelMap::q(double t) const {
    double den = std::max(1e-14, t1 - t0);
    double a = std::clamp((t - t0) / den, 0.0, 1.0);
    return (1.0 - a) * q0 + a * q1;
}

double PanelMap::qd(double) const {
    return (q1 - q0) / std::max(1e-14, t1 - t0);
}

double PanelMap::topY(double x, double t) const {
    return 1.0 + q(t) * mode(x);
}

Vector2d PanelMap::refToPhys(const Vector2d& X, double t) const {
    return Vector2d(X.x(), X.y() * topY(X.x(), t));
}

Vector2d PanelMap::velocityAt(double x, double y, double t) const {
    double h = std::max(1e-8, topY(x, t));
    double Y = y / h;
    return Vector2d(0.0, Y * qd(t) * mode(x));
}

void PanelMap::setLinearMotion(double a, double b, double qa, double qb) {
    t0 = a;
    t1 = b;
    q0 = qa;
    q1 = qb;
}

double PanelMap::maxMeshSpeed(double t) const { return std::abs(qd(t)); }

double panelAcceleration(const ModalState& s, double fluidForce, const PanelParams& p) {
    return (fluidForce - p.stiffness * s.q - p.damping * s.qd) / p.mass;
}

ModalState advancePanelSymplectic(ModalState s, double fluidForce,
                                  const PanelParams& p, double dt) {
    double a = panelAcceleration(s, fluidForce, p);
    s.qd += dt * a;
    s.q += dt * s.qd;
    if (s.q < p.qMin) { s.q = p.qMin; s.qd = std::max(0.0, s.qd); }
    if (s.q > p.qMax) { s.q = p.qMax; s.qd = std::min(0.0, s.qd); }
    return s;
}

int panelBoundaryTag(double x, double y, double t, const PanelMap& map) {
    double tol = 8e-4;
    if (std::abs(y) < tol) return TAG_SLIP_WALL;
    if (std::abs(x) < tol) return TAG_EXACT;
    if (std::abs(x - 1.0) < tol) return TAG_OUTFLOW;
    if (y > 0.45 && std::abs(y - map.topY(x, t)) < tol) return TAG_MOVING_WALL;
    return TAG_INTERIOR;
}

double panelGeneralizedForce(const Space& sp, const MatrixXd& U, const PanelMap& map,
                             double time, double pExt, double* meanPressure) {
    MatrixXd q1d;
    VectorXd w1d;
    sp.fem->quad1d(q1d, w1d);
    int nqe = static_cast<int>(w1d.size());
    int locDof = sp.fem->locDof;
    double force = 0.0, area = 0.0, pInt = 0.0;
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) != TAG_MOVING_WALL) continue;
        int t = (sp.e2s(e, 0) >= 0) ? sp.e2s(e, 0) : sp.e2s(e, 1);
        if (t < 0) continue;
        int n1 = sp.edge(e, 0), n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, t, n1, n2);
        Matrix<double, Dynamic, 4> Ue(locDof, 4);
        for (int i = 0; i < locDof; ++i) Ue.row(i) = U.row(sp.e2d(t, i));
        for (int q = 0; q < nqe; ++q) {
            double l1 = q1d(q, 0), l2 = q1d(q, 1);
            Vector3d lam = Vector3d::Zero();
            if (eo.et == 0) { lam(0) = (eo.dir == 0) ? l1 : l2; lam(1) = (eo.dir == 0) ? l2 : l1; }
            else if (eo.et == 1) { lam(1) = (eo.dir == 0) ? l1 : l2; lam(2) = (eo.dir == 0) ? l2 : l1; }
            else { lam(2) = (eo.dir == 0) ? l1 : l2; lam(0) = (eo.dir == 0) ? l2 : l1; }
            RowVectorXd phi = sp.fem->computeBasisValue_all(lam.transpose()).row(0);
            Vector4d Uq = (phi * Ue).transpose();
            double p = euler::pressure(Uq);
            Vector2d Pp = l1 * sp.mesh.node.row(n1).transpose() +
                          l2 * sp.mesh.node.row(n2).transpose();
            double ds = w1d(q) * eo.he;
            force += (p - pExt) * eo.nout.y() * map.mode(Pp.x()) * ds;
            pInt += p * ds;
            area += ds;
        }
    }
    if (meanPressure) *meanPressure = (area > 0.0) ? pInt / area : pExt;
    return force;
}
int runPanel(bool quick) {
    namespace fs = std::filesystem;
    int ord = 1;
    int nx = quick ? 26 : 36;
    int ny = quick ? 8 : 10;
    int maxGen = quick ? 1 : 2;
    int remeshEvery = quick ? 8 : 8;
    int nFrames = quick ? 32 : 96;
    double tEnd = quick ? 0.25 : 0.55;
    double cfl = 0.18;
    double uIn = 0.32;
    double thRef = 0.006, thCrs = 0.002;

    PanelParams params;
    ModalState panel;
    panel.q = 0.035;
    panel.qd = 0.0;

    PanelMap map;
    map.setLinearMotion(0.0, 1e-12, panel.q, panel.q);
    RefMapFn refMap = [&](const Vector2d& X, double time) { return map.refToPhys(X, time); };
    MaxMeshSpeedFn maxSpeed = [&](double time) { return map.maxMeshSpeed(time); };
    Tagger tagger = [&](double x, double y, double t) { return panelBoundaryTag(x, y, t, map); };
    MeshVelocityFn meshVel = [&](double x, double y, double t) { return map.velocityAt(x, y, t); };
    ALEBCFn bc = [&](double, double, double, const Vector4d& Um,
                     double nxn, double nyn, int tag, double wn) {
        if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL)
            return movingWallGhost(Um, nxn, nyn, (tag == TAG_MOVING_WALL) ? wn : 0.0);
        if (tag == TAG_EXACT)
            return euler::primToCons(1.0, uIn, 0.0, 1.0);
        return Um;
    };

    Mesh base;
    euler::makeRectMesh(base, 0.0, 1.0, 0.0, 1.0, nx, ny);
    ALEAdaptiveForest forest(base, ord, 4);
    Space sp;
    rebuildSpace(forest, ord, refMap, 0.0, tagger, sp);
    MatrixXd U = euler::projectInitial(*sp.fem, sp.mesh, sp.e2d,
                                       [=](double, double) { return Vector4d(1.0, uIn, 0.0, 1.0); });

    std::string dir = quick ? "out/panel_quick_frames" : "out/panel_frames";
    std::string csvPath = quick ? "out/panel_quick_diagnostics.csv" : "out/panel_diagnostics.csv";
    fs::create_directories(dir);
    fs::create_directories("out");
    for (const auto& e : fs::directory_iterator(dir))
        if (e.path().extension() == ".ppm") fs::remove(e.path());
    std::ofstream diag(csvPath);
    diag << "time,q,qd,generalized_force,mean_panel_pressure,triangles,rho_min,rho_max\n";

    auto setCurrentMap = [&](double time) {
        map.setLinearMotion(time, time + 1e-12, panel.q, panel.q);
    };
    auto writeFrame = [&](int idx, double time) {
        rebuildSpace(forest, ord, refMap, time, tagger, sp);
        char fn[512];
        std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", dir.c_str(), idx);
        int W = quick ? 800 : 1300;
        int H = quick ? 420 : 680;
        euler::writeScalarPPM(fn, *sp.fem, sp.mesh, sp.e2d, U.col(0), W, H,
                              -0.03, 1.03, -0.04, 1.10, 0.88, 1.12, euler::CM_COOLWARM);
        overlayMesh(fn, sp.mesh, -0.03, 1.03, -0.04, 1.10);
    };

    std::cout << "Body-fitted ALE-AMR aeroelastic panel demo: first-mode flexible top wall\n";
    std::cout << "  base=" << nx << "x" << ny << " max_gen=" << maxGen << " dP" << ord
              << " inflow u=" << uIn << "\n";
    std::cout << "  panel: modal m=" << params.mass << " k=" << params.stiffness
              << " c=" << params.damping << " q0=" << panel.q << "\n";

    double t = 0.0, nextFrame = 0.0, frameDt = tEnd / std::max(1, nFrames);
    int step = 0, frame = 0, lastRef = 0, lastCrs = 0;
    setCurrentMap(t);
    writeFrame(frame++, t);
    nextFrame += frameDt;

    while (t < tEnd - 1e-14) {
        setCurrentMap(t);
        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        double pMean = params.pExt;
        double fFluid = panelGeneralizedForce(sp, U, map, t, params.pExt, &pMean);
        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        dt = std::min(dt, 0.0025);

        ModalState nextPanel = advancePanelSymplectic(panel, fFluid, params, dt);
        map.setLinearMotion(t, t + dt, panel.q, nextPanel.q);
        U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
        panel = nextPanel;
        t += dt;
        ++step;

        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        if (step % remeshEvery == 0 && t < tEnd - 1e-14) {
            std::vector<int> flag = computeAMRFlags(sp, U, forest, maxGen, thRef, thCrs, 2, true);
            forest.syncFromState(U, sp.e2d, sp.fem->locDof);
            auto rc = forest.adapt(flag, maxGen);
            lastRef = rc.first;
            lastCrs = rc.second;
            rebuildSpace(forest, ord, refMap, t, tagger, sp);
            U = forest.gatherState(sp.e2d, sp.fem->locDof, sp.nDof);
        }

        double rmin = U.col(0).minCoeff(), rmax = U.col(0).maxCoeff();
        diag << std::setprecision(12) << t << "," << panel.q << "," << panel.qd << ","
             << fFluid << "," << pMean << "," << sp.mesh.elem.rows() << ","
             << rmin << "," << rmax << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            int gmax = 0;
            for (int k = 0; k < sp.mesh.elem.rows(); ++k) gmax = std::max(gmax, forest.gen(k));
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step << " q=" << panel.q << " qd=" << panel.qd
                      << " Fq=" << fFluid << " ptop=" << pMean
                      << " tris=" << sp.mesh.elem.rows()
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " gmax=" << gmax << " remesh(+" << lastRef << ",-" << lastCrs << ")\n";
            nextFrame += frameDt;
        }
    }

    std::string still = quick ? "out/panel_quick.ppm" : "out/panel.ppm";
    rebuildSpace(forest, ord, refMap, t, tagger, sp);
    euler::writeScalarPPM(still, *sp.fem, sp.mesh, sp.e2d, U.col(0), quick ? 800 : 1300, quick ? 420 : 680,
                          -0.03, 1.03, -0.04, 1.10, 0.88, 1.12, euler::CM_COOLWARM);
    overlayMesh(still, sp.mesh, -0.03, 1.03, -0.04, 1.10);
    std::cout << "Done. frames=" << frame << ", still=" << still
              << ", diagnostics=" << csvPath << "\n";
    std::cout << "ffmpeg -y -framerate 25 -i " << dir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 "
              << (quick ? "out/panel_quick.mp4" : "out/panel.mp4") << "\n";
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
    }
    return euler_ale::runPanel(quick);
}
