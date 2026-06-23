#include "FSI.h"

#include "PistonModel.h"

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

int runFSI(bool quick) {
    namespace fs = std::filesystem;
    int ord = 1;
    int nx = quick ? 20 : 34;
    int ny = quick ? 6 : 10;
    int maxGen = quick ? 1 : 2;
    int remeshEvery = quick ? 8 : 6;
    int nFrames = quick ? 28 : 120;
    double tEnd = quick ? 0.28 : 0.85;
    double cfl = 0.20;
    double thRef = 0.010, thCrs = 0.0035;

    FSIParams fsi;
    fsi.mass = 1.0;
    fsi.stiffness = 52.0;
    fsi.damping = 0.85;
    fsi.pExt = 1.0;
    PistonState body;
    body.x = 0.045;
    body.v = 0.0;

    PistonMap map;
    map.setLinearMotion(0.0, 1e-12, body.x, body.x);
    RefMapFn refMap = [&](const Vector2d& X, double time) { return map.refToPhys(X, time); };
    MaxMeshSpeedFn maxSpeed = [&](double time) { return map.maxMeshSpeed(time); };
    Tagger tagger = [&](double x, double y, double t) { return pistonBoundaryTag(x, y, t, map); };
    MeshVelocityFn meshVel = [&](double x, double y, double t) { return map.velocityAt(x, y, t); };
    ALEBCFn bc = [&](double, double, double, const Vector4d& Um,
                     double nxn, double nyn, int tag, double wn) {
        if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL)
            return movingWallGhost(Um, nxn, nyn, (tag == TAG_MOVING_WALL) ? wn : 0.0);
        return Um;
    };

    Mesh base;
    euler::makeRectMesh(base, 0.0, 1.0, 0.0, 1.0, nx, ny);
    ALEAdaptiveForest forest(base, ord, 4);
    Space sp;
    rebuildSpace(forest, ord, refMap, 0.0, tagger, sp);
    MatrixXd U = euler::projectInitial(*sp.fem, sp.mesh, sp.e2d,
                                       [](double, double) { return Vector4d(1.0, 0.0, 0.0, 1.0); });

    std::string dir = quick ? "out/fsi_quick_frames" : "out/fsi_frames";
    std::string csvPath = quick ? "out/fsi_quick_diagnostics.csv" : "out/fsi_diagnostics.csv";
    fs::create_directories(dir);
    fs::create_directories("out");
    for (const auto& e : fs::directory_iterator(dir))
        if (e.path().extension() == ".ppm") fs::remove(e.path());
    std::ofstream diag(csvPath);
    diag << "time,x,v,fluid_force,mean_wall_pressure,triangles,rho_min,rho_max\n";

    auto setCurrentMap = [&](double time) {
        map.setLinearMotion(time, time + 1e-12, body.x, body.x);
    };
    auto writeFrame = [&](int idx, double time) {
        rebuildSpace(forest, ord, refMap, time, tagger, sp);
        char fn[512];
        std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", dir.c_str(), idx);
        int W = quick ? 720 : 1200;
        int H = quick ? 360 : 600;
        euler::writeScalarPPM(fn, *sp.fem, sp.mesh, sp.e2d, U.col(0), W, H,
                              -0.10, 1.04, -0.04, 1.04, 0.82, 1.18, euler::CM_COOLWARM);
        overlayMesh(fn, sp.mesh, -0.10, 1.04, -0.04, 1.04);
    };

    std::cout << "Two-way body-fitted ALE-AMR FSI demo: mass-spring piston, base "
              << nx << "x" << ny << ", max_gen=" << maxGen << ", dP" << ord << "\n";
    std::cout << "  structure: m=" << fsi.mass << " k=" << fsi.stiffness
              << " c=" << fsi.damping << " x0=" << body.x << " p_ext=" << fsi.pExt << "\n";

    double t = 0.0, nextFrame = 0.0, frameDt = tEnd / std::max(1, nFrames);
    int step = 0, frame = 0, lastRef = 0, lastCrs = 0;
    setCurrentMap(t);
    writeFrame(frame++, t);
    nextFrame += frameDt;

    while (t < tEnd - 1e-14) {
        setCurrentMap(t);
        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        double pMean = fsi.pExt;
        double fFluid = pistonPressureForce(sp, U, fsi.pExt, &pMean);
        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        dt = std::min(dt, 0.004);

        PistonState nextBody = advanceStructureSymplectic(body, fFluid, fsi, dt);
        map.setLinearMotion(t, t + dt, body.x, nextBody.x);
        U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
        body = nextBody;
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
        diag << std::setprecision(12) << t << "," << body.x << "," << body.v << ","
             << fFluid << "," << pMean << "," << sp.mesh.elem.rows() << ","
             << rmin << "," << rmax << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            int gmax = 0;
            for (int k = 0; k < sp.mesh.elem.rows(); ++k) gmax = std::max(gmax, forest.gen(k));
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step << " x=" << std::setprecision(4) << body.x
                      << " v=" << body.v << " F=" << fFluid << " pwall=" << pMean
                      << " tris=" << sp.mesh.elem.rows()
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " gmax=" << gmax << " remesh(+" << lastRef << ",-" << lastCrs << ")\n";
            nextFrame += frameDt;
        }
    }

    std::string still = quick ? "out/fsi_quick.ppm" : "out/fsi.ppm";
    rebuildSpace(forest, ord, refMap, t, tagger, sp);
    euler::writeScalarPPM(still, *sp.fem, sp.mesh, sp.e2d, U.col(0), quick ? 720 : 1200, quick ? 360 : 600,
                          -0.10, 1.04, -0.04, 1.04, 0.82, 1.18, euler::CM_COOLWARM);
    overlayMesh(still, sp.mesh, -0.10, 1.04, -0.04, 1.04);
    std::cout << "Done. frames=" << frame << ", still=" << still
              << ", diagnostics=" << csvPath << "\n";
    std::cout << "ffmpeg -y -framerate 25 -i " << dir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 "
              << (quick ? "out/fsi_quick.mp4" : "out/fsi.mp4") << "\n";
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
    }
    return euler_ale::runFSI(quick);
}
