#include "Piston.h"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

namespace euler_ale {

using namespace Eigen;

int runPiston(bool quick) {
    namespace fs = std::filesystem;
    int ord = 1;
    int nx = quick ? 18 : 36;
    int ny = quick ? 6 : 10;
    int maxGen = 1;
    int remeshEvery = quick ? 8 : 6;
    int nFrames = quick ? 18 : 90;
    double tEnd = quick ? 0.15 : 0.35;
    double cfl = 0.22;
    double thRef = 0.012, thCrs = 0.004;
    PistonMap map;
    map.amp = 0.055;
    map.period = 0.55;
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

    std::string dir = quick ? "out/piston_quick_frames" : "out/piston_frames";
    fs::create_directories(dir);
    for (const auto& e : fs::directory_iterator(dir))
        if (e.path().extension() == ".ppm") fs::remove(e.path());

    auto writeFrame = [&](int idx, double time) {
        rebuildSpace(forest, ord, refMap, time, tagger, sp);
        char fn[512];
        std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", dir.c_str(), idx);
        int W = quick ? 720 : 1200;
        int H = quick ? 360 : 600;
        euler::writeScalarPPM(fn, *sp.fem, sp.mesh, sp.e2d, U.col(0), W, H,
                              -0.08, 1.04, -0.04, 1.04, 0.92, 1.10, euler::CM_INFERNO);
        overlayMesh(fn, sp.mesh, -0.08, 1.04, -0.04, 1.04);
    };

    std::cout << "Body-fitted ALE-AMR piston demo: base " << nx << "x" << ny
              << ", max_gen=" << maxGen << ", dP" << ord << "\n";
    double t = 0.0, nextFrame = 0.0, frameDt = tEnd / std::max(1, nFrames);
    int step = 0, frame = 0, lastRef = 0, lastCrs = 0;
    writeFrame(frame++, t);
    nextFrame += frameDt;
    while (t < tEnd - 1e-14) {
        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
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

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            double rmin = U.col(0).minCoeff(), rmax = U.col(0).maxCoeff();
            int gmax = 0;
            for (int k = 0; k < sp.mesh.elem.rows(); ++k) gmax = std::max(gmax, forest.gen(k));
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step << " tris=" << sp.mesh.elem.rows()
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " gmax=" << gmax << " remesh(+" << lastRef << ",-" << lastCrs << ")\n";
            nextFrame += frameDt;
        }
    }

    std::string still = quick ? "out/piston_quick.ppm" : "out/piston.ppm";
    euler::writeScalarPPM(still, *sp.fem, sp.mesh, sp.e2d, U.col(0), quick ? 720 : 1200, quick ? 360 : 600,
                          -0.08, 1.04, -0.04, 1.04, 0.92, 1.10, euler::CM_INFERNO);
    overlayMesh(still, sp.mesh, -0.08, 1.04, -0.04, 1.04);
    std::cout << "Done. frames=" << frame << ", still=" << still << "\n";
    std::cout << "ffmpeg -y -framerate 25 -i " << dir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 "
              << (quick ? "out/piston_quick.mp4" : "out/piston.mp4") << "\n";
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
    }
    return euler_ale::runPiston(quick);
}
