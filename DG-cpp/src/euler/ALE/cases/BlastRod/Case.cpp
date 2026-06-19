#include "Case.h"

#include "MeshGen.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace euler_ale {

int runBlastRod(bool quick, bool freshStart) {
    namespace fs = std::filesystem;
    std::cout << std::unitbuf;

    int ord = 1;
    int nFrames = quick ? 36 : 900;
    double tEnd = quick ? 0.24 : 2.16;
    double cfl = quick ? 0.18 : 0.15;
    double h = quick ? 0.024 : 0.012;
    int meshIter = quick ? 120 : 260;
    int renderW = quick ? 900 : 1500;
    int renderH = quick ? 600 : 1000;
    int ssaa = quick ? 2 : 3;
    const std::string outputPrefix = quick ? "blast_rod_quick" : "blast_rod";

    BlastRodGeom geom;
    SolidMaterial material;
    material.density = 70.0;
    material.young = 2.8e2;
    material.poisson = 0.34;
    material.damping = 0.55;
    material.velocityLimit = 2.6;
    material.displacementLimit = 0.24;

    ElasticSolid2D solid;
    int solidNx = quick ? 5 : 8;
    int solidNy = quick ? 56 : 112;
    auto buildInitialSolid = [&]() {
        solid.buildRectangularBeam(geom.rodLeft(), geom.rodRight(), geom.rodBaseY,
                                   geom.rodTipY(), solidNx, solidNy, material);
    };
    buildInitialSolid();

    SolidALEMap map;
    map.setSolid(&solid);
    map.setDomain(geom.xa, geom.xb, geom.ya, geom.yb);
    map.setInfluence(0.070, 0.075);
    MatrixXd aleReferenceNodes = solid.referenceNodes();
    map.setReferenceNodes(aleReferenceNodes);
    map.setCurrent(0.0, solid.currentNodes(), solid.velocities());
    const double pExt = 1.0;

    std::cout << "Body-fitted ALE static-mesh high-pressure gas over 2D elastic solid beam\n";
    std::cout << "  channel=" << geom.xb << "x" << geom.yb
              << " rod x=" << geom.rodX << " length=" << geom.rodL
              << " width=" << geom.rodW << "\n";
    std::cout << "  DG order=P" << ord << "\n";
    std::cout << "  generating common SDF mesh h_far=" << h << "...\n";
    Mesh base = makeBlastRodMesh(geom, h, meshIter, true);
    double minAng = 0.0, meanAng = 0.0, minArea = 0.0, maxArea = 0.0;
    meshQuality(base, minAng, meanAng, minArea, maxArea);

    RefMapFn refMap = [&](const Vector2d& X, double time) { return map.refToPhys(X, time); };
    MaxMeshSpeedFn maxSpeed = [&](double time) { return map.maxMeshSpeed(time); };
    auto inflowPrim = [](double time) {
        double r = blastRodSmoothRamp(time / 0.035);
        double rho = 1.0 + 0.72 * r;
        double u = 0.90 * r;
        double p = 1.0 + 2.35 * r;
        return Vector4d(rho, u, 0.0, p);
    };
    Tagger tagger = [&](double x, double y, double t) {
        return blastRodBoundaryTag(x, y, t, map, geom);
    };
    MeshVelocityFn meshVel = [&](double x, double y, double t) {
        return map.velocityAt(x, y, t);
    };
    ALEBCFn bc = [&](double, double, double time, const Vector4d& Um,
                     double nxn, double nyn, int tag, double wn) {
        if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL)
            return movingWallGhost(Um, nxn, nyn, (tag == TAG_MOVING_WALL) ? wn : 0.0);
        if (tag == TAG_EXACT) {
            Vector4d pr = inflowPrim(time);
            return euler::primToCons(pr(0), pr(1), pr(2), pr(3));
        }
        return Um;
    };

    ALEAdaptiveForest forest(base, ord, 4);
    StaticSpace stepSpace;
    StaticSpace nextSpace;
    StaticSpace renderSpace;
    initializeStaticSpace(forest, ord, refMap, 0.0, tagger, stepSpace);
    initializeStaticSpace(forest, ord, refMap, 0.0, tagger, nextSpace);
    initializeStaticSpace(forest, ord, refMap, 0.0, tagger, renderSpace);
    Space& sp = stepSpace.space;

    auto initialPrim = [=](double x, double) {
        if (x < geom.rodLeft() - 0.055) return Vector4d(1.72, 0.10, 0.0, 3.35);
        return Vector4d(1.0, 0.0, 0.0, 1.0);
    };
    MatrixXd U;
    auto resetFluidToInitial = [&]() {
        U = MatrixXd::Zero(sp.nDof, 4);
        MatrixXd dofLam = sp.fem->lagrangeNodes();
        for (int elem = 0; elem < sp.mesh.elem.rows(); ++elem) {
            Vector2d p0 = sp.mesh.node.row(sp.mesh.elem(elem, 0));
            Vector2d p1 = sp.mesh.node.row(sp.mesh.elem(elem, 1));
            Vector2d p2 = sp.mesh.node.row(sp.mesh.elem(elem, 2));
            for (int i = 0; i < sp.fem->locDof; ++i) {
                Vector3d lam = dofLam.row(i).transpose();
                Vector2d p = lam(0) * p0 + lam(1) * p1 + lam(2) * p2;
                Vector4d pr = initialPrim(p.x(), p.y());
                U.row(sp.e2d(elem, i)) =
                    euler::primToCons(pr(0), pr(1), pr(2), pr(3)).transpose();
            }
        }
    };
    const double rhoFloor = quick ? 0.05 : 0.08;
    const double pFloor = quick ? 0.04 : 0.06;
    const double speedMax = 4.0;
    resetFluidToInitial();
    applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);

    std::string dir = "out/" + outputPrefix + "_frames";
    std::string dirNoMesh = "out/" + outputPrefix + "_nomesh_frames";
    std::string csvPath = "out/" + outputPrefix + "_diagnostics.csv";
    fs::create_directories("out");
    fs::create_directories(dir);
    fs::create_directories(dirNoMesh);
    const std::string diagHeader =
        "time,tip_displacement,tip_velocity,fluid_force_x,drag,mean_rod_pressure,"
        "fluid_triangles,solid_nodes,solid_triangles,min_h,remesh_count,rho_min,rho_max\n";
    const double initialMinH = sp.minH;

    auto setCurrentMap = [&](double time) {
        map.setCurrent(time, solid.currentNodes(), solid.velocities());
    };
    auto renderViews = [&](const std::string& noMeshPath,
                           const std::string& meshPath,
                           double time) {
        updateStaticSpace(renderSpace, refMap, time);
        Space& rsp = renderSpace.space;
        int hiW = renderW * ssaa;
        int hiH = renderH * ssaa;
        std::vector<unsigned char> hiNoMesh =
            euler::renderScalarPPMImage(*rsp.fem, rsp.mesh, rsp.e2d, U.col(0),
                                        hiW, hiH, geom.xa - 0.02, geom.xb + 0.02,
                                        geom.ya - 0.025, geom.yb + 0.025, 0.55, 3.25,
                                        euler::CM_INFERNO);
        int outW = 0, outH = 0;
        std::vector<unsigned char> noMesh =
            downsampleImage(hiNoMesh, hiW, hiH, ssaa, outW, outH);
        if (!noMesh.empty()) {
            writePPM(noMeshPath, outW, outH, noMesh);
        }

        std::vector<unsigned char> hiMesh = hiNoMesh;
        overlayMesh(hiMesh, hiW, hiH, rsp.mesh, geom.xa - 0.02, geom.xb + 0.02,
                    geom.ya - 0.025, geom.yb + 0.025);
        overlaySolidMesh(hiMesh, hiW, hiH, solid, geom.xa - 0.02, geom.xb + 0.02,
                         geom.ya - 0.025, geom.yb + 0.025, true);
        outW = 0;
        outH = 0;
        std::vector<unsigned char> meshImage =
            downsampleImage(hiMesh, hiW, hiH, ssaa, outW, outH);
        if (!meshImage.empty()) {
            writePPM(meshPath, outW, outH, meshImage);
        }
    };
    auto writeFrame = [&](int idx, double time) {
        char fn[512];
        char fnNoMesh[512];
        std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", dir.c_str(), idx);
        std::snprintf(fnNoMesh, sizeof(fnNoMesh), "%s/frame_%05d.ppm",
                      dirNoMesh.c_str(), idx);
        renderViews(fnNoMesh, fn, time);
    };

    std::cout << "  base nodes=" << base.node.rows() << " elems=" << base.elem.rows()
              << " min_angle=" << std::setprecision(3) << minAng
              << " mean_angle=" << meanAng << " static local refinement near rod\n";
    std::cout << "  solid beam FEM: nodes=" << solid.numNodes()
              << " tris=" << solid.numElements()
              << " mass=" << solid.totalMass()
              << " E=" << material.young
              << " nu=" << material.poisson
              << " damping=" << material.damping
              << " p_high=3.35\n";
    std::cout << "  render SSAA=" << ssaa << "x"
              << " final=" << renderW << "x" << renderH
              << " internal=" << renderW * ssaa << "x" << renderH * ssaa
              << "\n";
    const double remeshMinH = (quick ? 0.35 : 0.50) * initialMinH;
    int remeshCount = 0;
    std::array<int, 3> checkpointDone{{0, 0, 0}};
    double t = 0.0;
    double nextFrame = 0.0;
    double frameDt = tEnd / std::max(1, nFrames);
    int step = 0;
    int frame = 0;
    bool resumed = false;
    std::cout << "  remesh guard: initial_min_h=" << initialMinH
              << " remesh_min_h=" << remeshMinH << "\n";

    auto restoreColdStartState = [&]() {
        buildInitialSolid();
        aleReferenceNodes = solid.referenceNodes();
        map.setReferenceNodes(aleReferenceNodes);
        map.setCurrent(0.0, solid.currentNodes(), solid.velocities());
        forest = ALEAdaptiveForest(base, ord, 4);
        initializeStaticSpace(forest, ord, refMap, 0.0, tagger, stepSpace);
        initializeStaticSpace(forest, ord, refMap, 0.0, tagger, nextSpace);
        initializeStaticSpace(forest, ord, refMap, 0.0, tagger, renderSpace);
        resetFluidToInitial();
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
    };

    if (!freshStart) {
        std::optional<RunCheckpoint> cp =
            loadLatestCheckpoint("blast_rod", quick, ord, nFrames, tEnd, h,
                                 solid.numNodes());
        if (cp.has_value()) {
            solid.setState(cp->solidNodes, cp->solidVelocities);
            aleReferenceNodes = cp->aleReferenceNodes;
            map.setReferenceNodes(aleReferenceNodes);
            map.setCurrent(cp->time, solid.currentNodes(), solid.velocities());
            forest = ALEAdaptiveForest(cp->referenceMesh, ord, 4);
            initializeStaticSpace(forest, ord, refMap, cp->time, tagger, stepSpace);
            initializeStaticSpace(forest, ord, refMap, cp->time, tagger, nextSpace);
            initializeStaticSpace(forest, ord, refMap, cp->time, tagger, renderSpace);
            if (cp->U.rows() == sp.nDof && cp->U.cols() == 4) {
                U = cp->U;
                applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
                t = cp->time;
                nextFrame = cp->nextFrame;
                step = cp->step;
                frame = cp->frame;
                remeshCount = cp->remeshCount;
                checkpointDone = cp->milestoneDone;
                pruneFramesFrom(dir, frame);
                pruneFramesFrom(dirNoMesh, frame);
                trimDiagnosticsToTime(csvPath, t);
                resumed = true;
                std::cout << "  resumed from checkpoint at t=" << std::fixed
                          << std::setprecision(5) << t
                          << " step=" << step
                          << " next_frame_index=" << frame
                          << " remesh=" << remeshCount << "\n";
            } else {
                std::cerr << "Warning: checkpoint DG state has incompatible size; cold start\n";
                restoreColdStartState();
            }
        }
    } else {
        std::cout << "  fresh start requested; checkpoint auto-resume disabled\n";
    }

    std::ofstream diag;
    if (resumed) {
        diag.open(csvPath, std::ios::app);
        std::error_code ec;
        if (!fs::exists(csvPath, ec) || fs::file_size(csvPath, ec) == 0) diag << diagHeader;
    } else {
        clearFrameDirectory(dir);
        clearFrameDirectory(dirNoMesh);
        diag.open(csvPath, std::ios::trunc);
        diag << diagHeader;
        setCurrentMap(t);
        writeFrame(frame++, t);
        nextFrame += frameDt;
    }
    if (!diag) {
        std::cerr << "Error: cannot open diagnostics " << csvPath << "\n";
        return 2;
    }

    auto remeshFluidAtCurrent = [&](double time, const std::string& reason) {
        Space oldSpace = std::move(stepSpace.space);
        MatrixXd oldU = U;
        aleReferenceNodes = solid.currentNodes();
        map.setReferenceNodes(aleReferenceNodes);
        map.setCurrent(time, solid.currentNodes(), solid.velocities());

        Mesh newBase = makeCurrentSolidBlastRodMesh(geom, solid, h, meshIter, false);
        double rMinAng = 0.0, rMeanAng = 0.0, rMinArea = 0.0, rMaxArea = 0.0;
        meshQuality(newBase, rMinAng, rMeanAng, rMinArea, rMaxArea);
        forest = ALEAdaptiveForest(newBase, ord, 4);
        initializeStaticSpace(forest, ord, refMap, time, tagger, stepSpace);
        initializeStaticSpace(forest, ord, refMap, time, tagger, nextSpace);
        initializeStaticSpace(forest, ord, refMap, time, tagger, renderSpace);
        U = interpolateDGToSpace(oldSpace, oldU, stepSpace.space);
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
        ++remeshCount;
        std::cout << "  remesh#" << remeshCount << " at t=" << std::fixed
                  << std::setprecision(5) << time
                  << " reason=" << reason
                  << " elems=" << newBase.elem.rows()
                  << " min_angle=" << std::setprecision(3) << rMinAng
                  << " minH=" << sp.minH << "\n";
    };

    auto makeCheckpoint = [&]() {
        RunCheckpoint cp;
        cp.quick = quick;
        cp.ord = ord;
        cp.nFrames = nFrames;
        cp.tEnd = tEnd;
        cp.h = h;
        cp.time = t;
        cp.nextFrame = nextFrame;
        cp.step = step;
        cp.frame = frame;
        cp.remeshCount = remeshCount;
        cp.milestoneDone = checkpointDone;
        cp.referenceMesh = stepSpace.referenceMesh;
        cp.U = U;
        cp.solidNodes = solid.currentNodes();
        cp.solidVelocities = solid.velocities();
        cp.aleReferenceNodes = aleReferenceNodes;
        return cp;
    };

    auto writeMilestoneCheckpoints = [&]() {
        double frac = (tEnd > 0.0) ? t / tEnd : 1.0;
        for (int i = 0; i < static_cast<int>(kCheckpointFractions.size()); ++i) {
            if (checkpointDone[i]) continue;
            if (frac + 1e-12 < kCheckpointFractions[i]) continue;
            checkpointDone[i] = 1;
            diag.flush();
            RunCheckpoint cp = makeCheckpoint();
            cp.milestoneDone = checkpointDone;
            fs::path path = checkpointPath("blast_rod", quick, kCheckpointLabels[i]);
            if (writeCheckpointAtomic(path, cp)) {
                std::cout << "  checkpoint " << kCheckpointLabels[i] << "% written at t="
                          << std::fixed << std::setprecision(5) << t
                          << " -> " << path.string() << "\n";
            }
        }
    };

    while (t < tEnd - 1e-14) {
        setCurrentMap(t);
        updateStaticSpace(stepSpace, refMap, t);
        if (sp.minH < remeshMinH) {
            remeshFluidAtCurrent(t, "pre_step_minH");
        }

        double pMean = pExt;
        double drag = 0.0;
        solid.clearExternalForces();
        double fFluid = blastRodLoadSolidFromFluid(sp, U, solid, pExt, &pMean, &drag);
        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        dt = std::min(dt, solid.stableTimeStep(0.38));
        dt = std::min(dt, quick ? 0.0012 : 0.0010);

        MatrixXd solidX0 = solid.currentNodes();
        MatrixXd solidV0 = solid.velocities();
        solid.advanceExplicit(dt);
        MatrixXd solidX1 = solid.currentNodes();
        MatrixXd solidV1 = solid.velocities();
        map.setMotion(t, t + dt, solidX0, solidV0, solidX1, solidV1);
        U = advanceOneStatic(stepSpace, nextSpace, refMap, t, dt, meshVel, bc, U);
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
        t += dt;
        ++step;

        updateStaticSpace(stepSpace, refMap, t);
        if (sp.minH < remeshMinH) {
            remeshFluidAtCurrent(t, "post_step_minH");
        }

        double rmin = U.col(0).minCoeff();
        double rmax = U.col(0).maxCoeff();
        diag << std::setprecision(12) << t << "," << solid.tipDisplacementX()
             << "," << solid.tipVelocityX() << ","
             << fFluid << "," << drag << "," << pMean << ","
             << sp.mesh.elem.rows() << "," << solid.numNodes() << ","
             << solid.numElements() << "," << sp.minH << "," << remeshCount
             << "," << rmin << "," << rmax << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step << " tipX=" << solid.tipDisplacementX()
                      << " tipV=" << solid.tipVelocityX()
                      << " Fq=" << fFluid << " drag=" << drag
                      << " tris=" << sp.mesh.elem.rows()
                      << " minH=" << sp.minH
                      << " remesh=" << remeshCount
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " solid_fem static_mesh\n";
            nextFrame += frameDt;
        }
        writeMilestoneCheckpoints();
    }

    std::string still = "out/" + outputPrefix + ".ppm";
    std::string stillNoMesh = "out/" + outputPrefix + "_nomesh.ppm";
    std::string stillPng = "out/" + outputPrefix + ".png";
    std::string stillNoMeshPng = "out/" + outputPrefix + "_nomesh.png";
    std::string video = "out/" + outputPrefix + ".mp4";
    std::string videoNoMesh = "out/" + outputPrefix + "_nomesh.mp4";
    renderViews(stillNoMesh, still, t);
    std::cout << "Done. frames=" << frame << ", still=" << still
              << ", no-mesh still=" << stillNoMesh
              << ", diagnostics=" << csvPath << "\n";

    int videoFps = quick ? 25 : 60;
    runOutputCommand("ffmpeg -y -i " + still + " -frames:v 1 -update 1 " + stillPng);
    runOutputCommand("ffmpeg -y -i " + stillNoMesh + " -frames:v 1 -update 1 " + stillNoMeshPng);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(videoFps) + " -i " + dir +
                     "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 " + video);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(videoFps) + " -i " + dirNoMesh +
                     "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 " + videoNoMesh);
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    bool freshStart = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
        if (std::string(argv[i]) == "--fresh") freshStart = true;
    }
    return euler_ale::runBlastRod(quick, freshStart);
}
