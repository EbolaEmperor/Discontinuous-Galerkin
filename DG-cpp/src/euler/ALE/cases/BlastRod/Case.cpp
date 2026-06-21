#include "Case.h"

#include "MeshGen.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace euler_ale {
namespace {

bool checkpointMatchesGeometry(const RunCheckpoint& cp, const BlastRodGeom& geom) {
    if (cp.referenceMesh.node.rows() == 0) return false;
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    for (int i = 0; i < cp.referenceMesh.node.rows(); ++i) {
        xmin = std::min(xmin, cp.referenceMesh.node(i, 0));
        xmax = std::max(xmax, cp.referenceMesh.node(i, 0));
    }
    double tol = std::max(1e-8, 0.02 * std::max(cp.h, 1e-12));
    return std::abs(xmin - geom.xa) <= tol && std::abs(xmax - geom.xb) <= tol;
}

int fullDomainRenderWidth(int physicalWidth, const BlastRodGeom& geom) {
    double physicalLength = std::max(geom.flowXb() - geom.xa, 1e-14);
    double computationalLength = std::max(geom.xb - geom.xa, physicalLength);
    int w = static_cast<int>(std::lround(physicalWidth * computationalLength / physicalLength));
    w = std::max(physicalWidth, w);
    if (w % 2) --w;
    return std::max(2, w);
}

struct MeshStats {
    double minH = 0.0;
    double minAngle = 0.0;
    double meanAngle = 0.0;
    double minArea = 0.0;
    double maxArea = 0.0;
    int minHElem = -1;
    Vector2d minHCentroid = Vector2d::Zero();
    double minHArea = 0.0;
    double minHLongest = 0.0;
    std::array<double, 3> minHEdges{{0.0, 0.0, 0.0}};
};

struct SolidBoundaryStats {
    double minMovingEdge = 0.0;
    int side = 0;
    int a = -1;
    int b = -1;
    Vector2d midpoint = Vector2d::Zero();
};

MeshStats currentMeshStats(const Space& space) {
    MeshStats s;
    s.minH = space.minH;
    meshQuality(space.mesh, s.minAngle, s.meanAngle, s.minArea, s.maxArea);
    double best = std::numeric_limits<double>::max();
    for (int e = 0; e < space.mesh.elem.rows(); ++e) {
        double h = hCFL(space.mesh, e);
        if (!(h < best)) continue;
        Vector2d p[3] = {
            space.mesh.node.row(space.mesh.elem(e, 0)).transpose(),
            space.mesh.node.row(space.mesh.elem(e, 1)).transpose(),
            space.mesh.node.row(space.mesh.elem(e, 2)).transpose()
        };
        best = h;
        s.minHElem = e;
        s.minH = h;
        s.minHCentroid = (p[0] + p[1] + p[2]) / 3.0;
        s.minHArea = triArea(space.mesh, e);
        s.minHLongest = longestEdge(space.mesh, e);
        s.minHEdges = {{(p[1] - p[0]).norm(), (p[2] - p[1]).norm(),
                        (p[0] - p[2]).norm()}};
    }
    return s;
}

SolidBoundaryStats movingBoundaryStats(const ElasticSolid2D& solid) {
    SolidBoundaryStats stats;
    stats.minMovingEdge = std::numeric_limits<double>::max();
    const MatrixXd& x = solid.currentNodes();
    for (const auto& seg : solid.movingBoundarySegments()) {
        if (seg.a < 0 || seg.b < 0 || seg.a >= x.rows() || seg.b >= x.rows()) continue;
        Vector2d a = x.row(seg.a).transpose();
        Vector2d b = x.row(seg.b).transpose();
        double len = (b - a).norm();
        if (len >= stats.minMovingEdge) continue;
        stats.minMovingEdge = len;
        stats.side = seg.side;
        stats.a = seg.a;
        stats.b = seg.b;
        stats.midpoint = 0.5 * (a + b);
    }
    if (stats.minMovingEdge == std::numeric_limits<double>::max())
        stats.minMovingEdge = 0.0;
    return stats;
}

std::vector<BlastRodMeshOptions> remeshCandidateOptions(double h) {
    const double hNear = 0.25 * h;
    const std::array<std::pair<double, double>, 8> offsets{{
        {0.00, 0.00}, {0.37, 0.11}, {0.13, 0.41}, {0.61, 0.29},
        {0.23, 0.67}, {0.79, 0.53}, {0.47, 0.83}, {0.89, 0.17}
    }};
    std::vector<BlastRodMeshOptions> opts;
    opts.reserve(offsets.size() + 2);
    for (int i = 0; i < static_cast<int>(offsets.size()); ++i) {
        BlastRodMeshOptions o;
        o.seedOffsetX = offsets[i].first * hNear;
        o.seedOffsetY = offsets[i].second * hNear;
        o.randomSeed = 12345u + 104729u * static_cast<unsigned>(i);
        if (i >= 4) o.gradeRadius = 0.42;
        opts.push_back(o);
    }
    const std::array<double, 3> coarseFactors{{0.30, 0.35, 0.40}};
    for (int i = 0; i < static_cast<int>(coarseFactors.size()); ++i) {
        BlastRodMeshOptions o;
        o.hNearFactor = coarseFactors[i];
        o.gradeRadius = 0.40 + 0.06 * i;
        o.seedOffsetX = (0.19 + 0.21 * i) * hNear;
        o.seedOffsetY = (0.31 + 0.17 * i) * hNear;
        o.randomSeed = 83492791u + 65537u * static_cast<unsigned>(i);
        opts.push_back(o);
    }
    BlastRodMeshOptions finer;
    finer.hNearFactor = 0.22;
    finer.gradeRadius = 0.42;
    finer.seedOffsetX = 0.31 * hNear;
    finer.seedOffsetY = 0.73 * hNear;
    finer.randomSeed = 998244353u;
    opts.push_back(finer);
    return opts;
}

} // namespace

int runBlastRodRemeshBench(bool quick) {
    int ord = 1;
    int nFrames = quick ? 36 : 900;
    double tEnd = quick ? 0.24 : 2.16;
    double h = quick ? 0.024 : 0.012;
    int meshIter = quick ? 120 : 260;

    BlastRodGeom geom;
    SolidMaterial material;
    material.density = 70.0;
    material.young = 2.8e2;
    material.poisson = 0.34;
    material.damping = 0.55;

    ElasticSolid2D solid;
    int solidNx = quick ? 8 : 12;
    int solidNy = quick ? 84 : 140;
    solid.buildRoundedRootBeam(geom.rodLeft(), geom.rodRight(), geom.rodBaseY,
                               geom.rodTipY(), geom.rodRootRadius(),
                               solidNx, solidNy, material);

    std::optional<RunCheckpoint> cp =
        loadLatestCheckpoint("blast_rod", quick, ord, nFrames, tEnd, h, solid.numNodes());
    if (cp.has_value() && !checkpointMatchesGeometry(*cp, geom)) {
        std::cerr << "Checkpoint geometry does not match current sponge domain; skipping\n";
        cp.reset();
    }
    if (!cp.has_value()) {
        std::cerr << "No compatible checkpoint found for remesh bench\n";
        return 2;
    }
    if (cp->solidReferenceMesh.node.rows() > 0) {
        solid.resetReferenceMesh(cp->solidReferenceMesh, material);
    }
    solid.setState(cp->solidNodes, cp->solidVelocities);
    SolidMeshQuality solidQ = solid.currentMeshQuality();
    SolidBoundaryStats bStats = movingBoundaryStats(solid);
    std::cout << "BlastRod remesh bench from checkpoint: t=" << std::setprecision(8)
              << cp->time << " step=" << cp->step
              << " remesh=" << cp->remeshCount
              << " h_far=" << h << " iter=" << meshIter << "\n"
              << "  checkpoint solid current quality: min_angle="
              << std::fixed << std::setprecision(3) << solidQ.minAngleDeg
              << " min_edge=" << std::scientific << std::setprecision(6)
              << solidQ.minEdge
              << " min_moving_boundary_edge=" << bStats.minMovingEdge
              << " inverted=" << solidQ.invertedElements
              << " edge_mid=(" << std::fixed << std::setprecision(6)
              << bStats.midpoint.x() << "," << bStats.midpoint.y() << ")"
              << std::defaultfloat << "\n";

    const double solidRemeshH = quick ? 0.0045 : 0.0030;
    const int solidRemeshIter = quick ? 140 : 260;
    int oldSolidNodes = solid.numNodes();
    int oldSolidElems = solid.numElements();
    Mesh newSolidCurrent =
        makeCurrentBlastRodSolidInteriorMesh(solid, solidRemeshH, solidRemeshIter, false);
    solid.remeshToCurrentMesh(newSolidCurrent);
    SolidMeshQuality remeshedSolidQ = solid.currentMeshQuality();
    SolidBoundaryStats remeshedBoundary = movingBoundaryStats(solid);
    std::cout << "  solid remesh bench: h=" << solidRemeshH
              << " iter=" << solidRemeshIter
              << " nodes=" << oldSolidNodes << "->" << solid.numNodes()
              << " elems=" << oldSolidElems << "->" << solid.numElements()
              << " min_angle=" << std::fixed << std::setprecision(3)
              << solidQ.minAngleDeg << "->" << remeshedSolidQ.minAngleDeg
              << " inverted=" << solidQ.invertedElements
              << "->" << remeshedSolidQ.invertedElements
              << " moving_edge=" << std::scientific << std::setprecision(6)
              << bStats.minMovingEdge << "->" << remeshedBoundary.minMovingEdge
              << " edge_mid=(" << std::fixed << std::setprecision(6)
              << remeshedBoundary.midpoint.x() << ","
              << remeshedBoundary.midpoint.y() << ")"
              << std::defaultfloat << "\n";

    auto start = std::chrono::steady_clock::now();
    Mesh mesh = makeCurrentSolidBlastRodMesh(geom, solid, h, meshIter, true);
    using seconds = std::chrono::duration<double>;
    double wall = std::chrono::duration_cast<seconds>(
        std::chrono::steady_clock::now() - start).count();
    double minAng = 0.0, meanAng = 0.0, minArea = 0.0, maxArea = 0.0;
    meshQuality(mesh, minAng, meanAng, minArea, maxArea);
    std::cout << "Remesh bench done: nodes=" << mesh.node.rows()
              << " elems=" << mesh.elem.rows()
              << " min_angle=" << minAng
              << " mean_angle=" << meanAng
              << " area=[" << minArea << "," << maxArea << "]"
              << " wall=" << wall << "s\n";
    return 0;
}

int runBlastRod(bool quick, bool freshStart) {
    namespace fs = std::filesystem;
    std::cout << std::unitbuf;

    int ord = 1;
    int nFrames = quick ? 36 : 900;
    double tEnd = quick ? 0.24 : 2.16;
    double cfl = quick ? 0.18 : 0.15;
    double h = quick ? 0.024 : 0.012;
    int meshIter = quick ? 120 : 260;
    int videoW = quick ? 900 : 1500;
    int renderH = quick ? 600 : 1000;
    int ssaa = quick ? 2 : 3;
    const std::string outputPrefix = quick ? "blast_rod_quick" : "blast_rod";

    BlastRodGeom geom;
    int renderW = fullDomainRenderWidth(videoW, geom);
    const int videoCropW = videoW;
    const double viewXa = geom.xa;
    const double viewXb = geom.xb;
    const double viewYa = geom.ya - 0.025;
    const double viewYb = geom.yb + 0.025;
    SolidMaterial material;
    material.density = 70.0;
    material.young = 2.8e2;
    material.poisson = 0.34;
    material.damping = 0.55;

    ElasticSolid2D solid;
    int solidNx = quick ? 8 : 12;
    int solidNy = quick ? 84 : 140;
    auto buildInitialSolid = [&]() {
        solid.buildRoundedRootBeam(geom.rodLeft(), geom.rodRight(), geom.rodBaseY,
                                   geom.rodTipY(), geom.rodRootRadius(),
                                   solidNx, solidNy, material);
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
    std::cout << "  physical channel=" << geom.flowXb() << "x" << geom.yb
              << " computational channel=" << geom.xb << "x" << geom.yb
              << " sponge=[" << geom.spongeStartX() << "," << geom.xb << "]"
              << " rod x=" << geom.rodX << " length=" << geom.rodL
              << " width=" << geom.rodW
              << " root_radius=" << geom.rodRootRadius() << "\n";
    std::cout << "  DG order=P" << ord << "\n";
    std::optional<RunCheckpoint> resumeCheckpoint;
    if (!freshStart) {
        resumeCheckpoint = loadLatestCheckpoint("blast_rod", quick, ord, nFrames, tEnd, h,
                                                solid.numNodes());
        if (resumeCheckpoint.has_value() &&
            !checkpointMatchesGeometry(*resumeCheckpoint, geom)) {
            std::cout << "  checkpoint geometry differs from current sponge domain; cold start\n";
            resumeCheckpoint.reset();
        }
    }

    Mesh base;
    if (resumeCheckpoint.has_value()) {
        base = resumeCheckpoint->referenceMesh;
        std::cout << "  using checkpoint reference mesh; skipping startup SDF mesh generation\n";
    } else {
        std::cout << "  generating common SDF mesh h_far=" << h << "...\n";
        base = makeBlastRodMesh(geom, h, meshIter, true);
    }
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
        if (tag == TAG_OUTFLOW) {
            return characteristicPressureOutletGhost(Um, nxn, nyn, wn,
                                                     Vector4d(1.0, 0.0, 0.0, pExt));
        }
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
    const Vector4d spongeReference = euler::primToCons(1.0, 0.0, 0.0, pExt);
    const double spongeSigmaMax =
        3.0 * (0.90 + std::sqrt(euler::GAMMA * pExt)) /
        std::max(geom.spongeWidth(), 1e-12);
    auto applyRightSponge = [&](MatrixXd& state, const Space& space, double dtStep) {
        if (geom.spongeWidth() <= 1e-14 || dtStep <= 0.0) return;
        MatrixXd dofLam = space.fem->lagrangeNodes();
        for (int elem = 0; elem < space.mesh.elem.rows(); ++elem) {
            Vector2d p0 = space.mesh.node.row(space.mesh.elem(elem, 0));
            Vector2d p1 = space.mesh.node.row(space.mesh.elem(elem, 1));
            Vector2d p2 = space.mesh.node.row(space.mesh.elem(elem, 2));
            for (int i = 0; i < space.fem->locDof; ++i) {
                Vector3d lam = dofLam.row(i).transpose();
                Vector2d p = lam(0) * p0 + lam(1) * p1 + lam(2) * p2;
                double s = geom.spongeCoordinate(p.x());
                if (s <= 0.0) continue;
                double ramp = s * s * (3.0 - 2.0 * s);
                double alpha = std::exp(-spongeSigmaMax * ramp * dtStep);
                int dof = space.e2d(elem, i);
                Vector4d Ui = state.row(dof).transpose();
                state.row(dof) = (spongeReference + alpha * (Ui - spongeReference)).transpose();
            }
        }
    };

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
                                        hiW, hiH, viewXa, viewXb, viewYa, viewYb, 0.55, 3.25,
                                        euler::CM_INFERNO);
        int outW = 0, outH = 0;
        std::vector<unsigned char> noMesh =
            downsampleImage(hiNoMesh, hiW, hiH, ssaa, outW, outH);
        if (!noMesh.empty()) {
            writePPM(noMeshPath, outW, outH, noMesh);
        }

        std::vector<unsigned char> hiMesh = hiNoMesh;
        overlayMesh(hiMesh, hiW, hiH, rsp.mesh, viewXa, viewXb, viewYa, viewYb);
        overlaySolidMesh(hiMesh, hiW, hiH, solid, viewXa, viewXb, viewYa, viewYb, true);
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
    SolidMeshQuality sq = solid.meshQuality();
    std::cout << "  solid beam FEM: nodes=" << solid.numNodes()
              << " tris=" << solid.numElements()
              << " min_angle=" << sq.minAngleDeg
              << " mean_angle=" << sq.meanAngleDeg
              << " min_edge=" << sq.minEdge
              << " inverted=" << sq.invertedElements
              << " area=[" << sq.minArea << "," << sq.maxArea << "]"
              << " mass=" << solid.totalMass()
              << " E=" << material.young
              << " nu=" << material.poisson
              << " damping=" << material.damping
              << " p_high=3.35\n";
    std::cout << "  render SSAA=" << ssaa << "x"
              << " final=" << renderW << "x" << renderH
              << " internal=" << renderW * ssaa << "x" << renderH * ssaa
              << " ppm_full_domain video_crop=" << videoCropW << "x" << renderH
              << " sponge_sigma_max=" << spongeSigmaMax
              << "\n";
    const double remeshSoftH = 0.50 * initialMinH;
    const double remeshHardH = 0.35 * initialMinH;
    const double remeshEmergencyH = 0.22 * initialMinH;
    const double remeshAcceptH = 0.58 * initialMinH;
    const double remeshAngleTrigger = 13.0;
    const double remeshAcceptAngle = 18.0;
    const int remeshCooldownSteps = quick ? 12 : 40;
    const int remeshAttemptCooldownSteps = quick ? 24 : 180;
    int lastRemeshStep = -1000000;
    int lastRemeshAttemptStep = -1000000;
    double lastFailedRemeshMinH = std::numeric_limits<double>::max();
    double lastAcceptedRemeshMinH = std::numeric_limits<double>::max();
    int remeshCount = 0;
    const double solidRemeshH = quick ? 0.0045 : 0.0030;
    const int solidRemeshIter = quick ? 140 : 260;
    const double solidRemeshAngleTrigger = 9.0;
    const double solidRemeshMovingEdgeTrigger =
        0.35 * std::max(sq.minEdge, 0.40 * solidRemeshH);
    const int solidRemeshCooldownSteps = quick ? 80 : 1400;
    int solidRemeshCount = 0;
    int lastSolidRemeshStep = -1000000;
    std::vector<CheckpointMilestone> checkpointPlan = checkpointSchedule(quick);
    std::vector<int> checkpointDone(checkpointPlan.size(), 0);
    double t = 0.0;
    double nextFrame = 0.0;
    double frameDt = tEnd / std::max(1, nFrames);
    int step = 0;
    int frame = 0;
    bool resumed = false;
    std::cout << "  remesh guard: initial_min_h=" << initialMinH
              << " soft_h=" << remeshSoftH
              << " hard_h=" << remeshHardH
              << " emergency_h=" << remeshEmergencyH
              << " accept_h=" << remeshAcceptH
              << " angle_trigger=" << remeshAngleTrigger
              << " accept_angle=" << remeshAcceptAngle
              << " cooldown_steps=" << remeshCooldownSteps
              << " failed_attempt_cooldown_steps=" << remeshAttemptCooldownSteps << "\n";
    std::cout << "  solid remesh guard: h=" << solidRemeshH
              << " iter=" << solidRemeshIter
              << " angle_trigger=" << solidRemeshAngleTrigger
              << " moving_edge_trigger=" << solidRemeshMovingEdgeTrigger
              << " cooldown_steps=" << solidRemeshCooldownSteps << "\n";

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
        if (resumeCheckpoint.has_value()) {
            const RunCheckpoint& cp = *resumeCheckpoint;
            if (cp.solidReferenceMesh.node.rows() > 0) {
                solid.resetReferenceMesh(cp.solidReferenceMesh, material);
                map.setSolid(&solid);
            }
            solid.setState(cp.solidNodes, cp.solidVelocities);
            aleReferenceNodes = cp.aleReferenceNodes;
            map.setReferenceNodes(aleReferenceNodes);
            map.setCurrent(cp.time, solid.currentNodes(), solid.velocities());
            forest = ALEAdaptiveForest(cp.referenceMesh, ord, 4);
            initializeStaticSpace(forest, ord, refMap, cp.time, tagger, stepSpace);
            initializeStaticSpace(forest, ord, refMap, cp.time, tagger, nextSpace);
            initializeStaticSpace(forest, ord, refMap, cp.time, tagger, renderSpace);
            if (cp.U.rows() == sp.nDof && cp.U.cols() == 4) {
                U = cp.U;
                applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
                t = cp.time;
                nextFrame = cp.nextFrame;
                step = cp.step;
                frame = cp.frame;
                remeshCount = cp.remeshCount;
                if (cp.milestoneDone.size() == checkpointDone.size()) {
                    checkpointDone = cp.milestoneDone;
                } else {
                    double frac = (tEnd > 0.0) ? cp.time / tEnd : 1.0;
                    for (int i = 0; i < static_cast<int>(checkpointPlan.size()); ++i)
                        checkpointDone[i] =
                            (frac + 1e-12 >= checkpointPlan[i].fraction) ? 1 : 0;
                }
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
        if (!quick) {
            pruneOldCheckpoints("blast_rod", quick, 0);
            std::cout << "  removed old full checkpoints for fresh run\n";
        }
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

    auto acceptRemeshCandidate = [&](double time, const std::string& reason,
                                     const Mesh& newBase,
                                     const MeshStats& stats,
                                     ALEAdaptiveForest&& candidateForest,
                                     StaticSpace&& candidateStep,
                                     const BlastRodMeshOptions& options,
                                     const std::string& qualityLabel) {
        Space oldSpace = std::move(stepSpace.space);
        MatrixXd oldU = U;
        forest = std::move(candidateForest);
        stepSpace = std::move(candidateStep);
        initializeStaticSpace(forest, ord, refMap, time, tagger, nextSpace);
        initializeStaticSpace(forest, ord, refMap, time, tagger, renderSpace);
        U = interpolateDGToSpace(oldSpace, oldU, stepSpace.space);
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
        ++remeshCount;
        lastRemeshStep = step;
        lastRemeshAttemptStep = step;
        lastFailedRemeshMinH = std::numeric_limits<double>::max();
        lastAcceptedRemeshMinH = stats.minH;
        SolidMeshQuality solidQ = solid.currentMeshQuality();
        SolidBoundaryStats bStats = movingBoundaryStats(solid);
        std::cout << "  remesh#" << remeshCount << " accepted at t=" << std::fixed
                  << std::setprecision(5) << time
                  << " reason=" << reason
                  << " quality=" << qualityLabel
                  << " elems=" << newBase.elem.rows()
                  << " min_angle=" << std::setprecision(3) << stats.minAngle
                  << " mean_angle=" << stats.meanAngle
                  << " minH=" << std::scientific << std::setprecision(6) << stats.minH
                  << " minH_xy=(" << std::fixed << std::setprecision(6)
                  << stats.minHCentroid.x() << "," << stats.minHCentroid.y() << ")"
                  << " minH_elem=" << stats.minHElem
                  << " minH_edges=(" << std::scientific << std::setprecision(6)
                  << stats.minHEdges[0] << "," << stats.minHEdges[1] << ","
                  << stats.minHEdges[2] << ")"
                  << " solid_min_angle=" << std::fixed << std::setprecision(3)
                  << solidQ.minAngleDeg
                  << " solid_inverted=" << solidQ.invertedElements
                  << " solid_min_moving_edge=" << std::scientific << std::setprecision(6)
                  << bStats.minMovingEdge
                  << " solid_edge_mid=(" << std::fixed << std::setprecision(6)
                  << bStats.midpoint.x() << "," << bStats.midpoint.y() << ")"
                  << " seed_offset=(" << options.seedOffsetX << "," << options.seedOffsetY << ")"
                  << " seed=" << options.randomSeed
                  << " h_near_factor=" << options.hNearFactor
                  << " grade_radius=" << options.gradeRadius
                  << std::defaultfloat << "\n";
    };

    auto remeshFluidAtCurrent = [&](double time, const std::string& reason,
                                    bool forceGeometryRefresh = false) -> bool {
        double oldMinH = sp.minH;
        MatrixXd previousReferenceNodes = aleReferenceNodes;
        aleReferenceNodes = solid.currentNodes();
        map.setReferenceNodes(aleReferenceNodes);
        map.setCurrent(time, solid.currentNodes(), solid.velocities());

        struct CandidateRecord {
            Mesh base;
            MeshStats stats;
            BlastRodMeshOptions options;
            int index = -1;
            double score = -1.0;
        };
        CandidateRecord best;
        std::vector<BlastRodMeshOptions> candidates = remeshCandidateOptions(h);
        double forcedMinH = std::max(0.35 * remeshAcceptH, 0.50 * remeshEmergencyH);
        double forcedMinAngle = 12.0;
        for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
            int iterBudget = meshIter + ((i >= 4) ? 80 : 0);
            Mesh newBase = makeCurrentSolidBlastRodMesh(geom, solid, h, iterBudget,
                                                        false, candidates[i]);
            ALEAdaptiveForest candidateForest(newBase, ord, 4);
            StaticSpace candidateStep;
            initializeStaticSpace(candidateForest, ord, refMap, time, tagger, candidateStep);
            MeshStats stats = currentMeshStats(candidateStep.space);
            double hScore = stats.minH / std::max(remeshAcceptH, 1e-14);
            double angleScore = stats.minAngle / std::max(remeshAcceptAngle, 1e-14);
            double score = std::min(hScore, 1.5) + std::min(angleScore, 1.5);
            if (score > best.score) {
                best.base = newBase;
                best.stats = stats;
                best.options = candidates[i];
                best.index = i;
                best.score = score;
            }

            bool accepted = stats.minH >= remeshAcceptH &&
                            stats.minAngle >= remeshAcceptAngle;
            bool forcedAccepted = forceGeometryRefresh &&
                                  stats.minH >= forcedMinH &&
                                  stats.minAngle >= forcedMinAngle;
            std::cout << "  remesh candidate " << i
                      << " reason=" << reason
                      << " elems=" << newBase.elem.rows()
                      << " min_angle=" << std::fixed << std::setprecision(3) << stats.minAngle
                      << " minH=" << std::scientific << std::setprecision(6) << stats.minH
                      << " h_ratio=" << stats.minH / std::max(oldMinH, 1e-300)
                      << " minH_xy=(" << std::fixed << std::setprecision(6)
                      << stats.minHCentroid.x() << "," << stats.minHCentroid.y() << ")"
                      << " seed_offset=(" << candidates[i].seedOffsetX << ","
                      << candidates[i].seedOffsetY << ")"
                      << " seed=" << candidates[i].randomSeed
                      << std::defaultfloat
                      << (accepted ? " accepted\n" :
                          (forcedAccepted ? " accepted_solid_geometry\n" : " rejected\n"));
            if (accepted || forcedAccepted) {
                acceptRemeshCandidate(time, reason, newBase, stats,
                                      std::move(candidateForest), std::move(candidateStep),
                                      candidates[i],
                                      accepted ? "strict" : "solid-geometry-limited");
                return true;
            }
        }

        bool softAngleRepair = reason.find("_soft_minH_angle") != std::string::npos &&
                               best.index >= 0 &&
                               best.stats.minAngle >= remeshAcceptAngle &&
                               best.stats.minH >= std::max(remeshHardH, 0.70 * oldMinH);
        if (softAngleRepair) {
            ALEAdaptiveForest bestForest(best.base, ord, 4);
            StaticSpace bestStep;
            initializeStaticSpace(bestForest, ord, refMap, time, tagger, bestStep);
            acceptRemeshCandidate(time, reason, best.base, best.stats,
                                  std::move(bestForest), std::move(bestStep),
                                  best.options, "angle-repair");
            return true;
        }

        bool marginal = best.index >= 0 &&
                        best.stats.minH >= std::max(remeshHardH * 1.05, oldMinH * 1.08) &&
                        best.stats.minAngle >= remeshAngleTrigger + 1.0;
        if (marginal) {
            ALEAdaptiveForest bestForest(best.base, ord, 4);
            StaticSpace bestStep;
            initializeStaticSpace(bestForest, ord, refMap, time, tagger, bestStep);
            acceptRemeshCandidate(time, reason, best.base, best.stats,
                                  std::move(bestForest), std::move(bestStep),
                                  best.options, "marginal");
            return true;
        }

        bool hardOrEmergency = reason.find("_hard_minH") != std::string::npos ||
                               reason.find("_emergency_minH") != std::string::npos;
        bool geometryLimited =
            best.index >= 0 &&
            ((hardOrEmergency &&
              best.stats.minAngle >= remeshAcceptAngle &&
              best.stats.minH >= 0.90 * oldMinH) ||
             (forceGeometryRefresh &&
              best.stats.minAngle >= forcedMinAngle &&
              best.stats.minH >= forcedMinH));
        if (geometryLimited) {
            ALEAdaptiveForest bestForest(best.base, ord, 4);
            StaticSpace bestStep;
            initializeStaticSpace(bestForest, ord, refMap, time, tagger, bestStep);
            acceptRemeshCandidate(time, reason, best.base, best.stats,
                                  std::move(bestForest), std::move(bestStep),
                                  best.options, "geometry-limited");
            return true;
        }

        SolidMeshQuality solidQ = solid.currentMeshQuality();
        SolidBoundaryStats bStats = movingBoundaryStats(solid);
        std::cout << "  remesh skipped at t=" << std::fixed << std::setprecision(5) << time
                  << " reason=" << reason
                  << " old_minH=" << std::scientific << std::setprecision(6) << oldMinH
                  << " best_candidate=" << best.index
                  << " best_minH=" << best.stats.minH
                  << " best_ratio=" << best.stats.minH / std::max(oldMinH, 1e-300)
                  << " best_min_angle=" << std::fixed << std::setprecision(3)
                  << best.stats.minAngle
                  << " best_minH_xy=(" << best.stats.minHCentroid.x()
                  << "," << best.stats.minHCentroid.y() << ")"
                  << " solid_min_angle=" << solidQ.minAngleDeg
                  << " solid_inverted=" << solidQ.invertedElements
                  << " solid_min_moving_edge=" << std::scientific << std::setprecision(6)
                  << bStats.minMovingEdge
                  << " solid_edge_mid=(" << std::fixed << std::setprecision(6)
                  << bStats.midpoint.x() << "," << bStats.midpoint.y() << ")"
                  << std::defaultfloat << " keeping_current_mesh\n";
        lastRemeshAttemptStep = step;
        lastFailedRemeshMinH = oldMinH;
        if (forceGeometryRefresh) {
            std::ostringstream oss;
            oss << "solid-triggered fluid remesh failed at t=" << std::setprecision(17)
                << time << " reason=" << reason
                << " best_candidate=" << best.index
                << " best_minH=" << best.stats.minH
                << " best_min_angle=" << best.stats.minAngle;
            throw std::runtime_error(oss.str());
        }
        aleReferenceNodes = previousReferenceNodes;
        map.setReferenceNodes(aleReferenceNodes);
        map.setCurrent(time, solid.currentNodes(), solid.velocities());
        return false;
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
        cp.solidReferenceMesh.node = solid.referenceNodes();
        cp.solidReferenceMesh.elem = solid.elements();
        cp.U = U;
        cp.solidNodes = solid.currentNodes();
        cp.solidVelocities = solid.velocities();
        cp.aleReferenceNodes = aleReferenceNodes;
        return cp;
    };

    auto writeMilestoneCheckpoints = [&]() {
        double frac = (tEnd > 0.0) ? t / tEnd : 1.0;
        for (int i = 0; i < static_cast<int>(checkpointPlan.size()); ++i) {
            if (checkpointDone[i]) continue;
            if (frac + 1e-12 < checkpointPlan[i].fraction) continue;
            checkpointDone[i] = 1;
            diag.flush();
            RunCheckpoint cp = makeCheckpoint();
            cp.milestoneDone = checkpointDone;
            fs::path path = checkpointPath("blast_rod", quick, checkpointPlan[i].label);
            if (writeCheckpointAtomic(path, cp)) {
                if (!quick) pruneOldCheckpoints("blast_rod", quick, 3);
                std::cout << "  checkpoint " << checkpointPlan[i].label << "% written at t="
                          << std::fixed << std::setprecision(5) << t
                          << " -> " << path.string() << "\n";
            }
        }
    };

    auto remeshNeeded = [&](const std::string& phase, std::string& reason) {
        if (!(sp.minH > 0.0) || !std::isfinite(sp.minH)) {
            reason = phase + "_invalid_minH";
            return true;
        }
        if (sp.minH < remeshEmergencyH) {
            bool recentAccepted = step - lastRemeshStep < remeshAttemptCooldownSteps;
            bool acceptedSubstantiallyWorse =
                lastAcceptedRemeshMinH < 0.5 * std::numeric_limits<double>::max() &&
                sp.minH < 0.85 * lastAcceptedRemeshMinH;
            if (recentAccepted && !acceptedSubstantiallyWorse) return false;
            bool recentFailure = step - lastRemeshAttemptStep < remeshAttemptCooldownSteps;
            bool substantiallyWorse =
                lastFailedRemeshMinH < 0.5 * std::numeric_limits<double>::max() &&
                sp.minH < 0.85 * lastFailedRemeshMinH;
            if (recentFailure && !substantiallyWorse) return false;
            reason = phase + "_emergency_minH";
            return true;
        }
        if (step - lastRemeshAttemptStep < remeshAttemptCooldownSteps) return false;
        if (step - lastRemeshStep < remeshCooldownSteps) return false;
        if (sp.minH < remeshHardH) {
            reason = phase + "_hard_minH";
            return true;
        }
        if (sp.minH < remeshSoftH) {
            MeshStats stats = currentMeshStats(sp);
            if (stats.minAngle < remeshAngleTrigger) {
                std::ostringstream oss;
                oss << phase << "_soft_minH_angle"
                    << "_angle=" << std::fixed << std::setprecision(3) << stats.minAngle;
                reason = oss.str();
                return true;
            }
        }
        return false;
    };

    auto solidRemeshNeeded = [&](std::string& reason) {
        if (step - lastSolidRemeshStep < solidRemeshCooldownSteps) return false;
        SolidMeshQuality q = solid.currentMeshQuality();
        SolidBoundaryStats bStats = movingBoundaryStats(solid);
        if (q.invertedElements > 0) {
            std::ostringstream oss;
            oss << "solid_inverted=" << q.invertedElements
                << "_angle=" << std::fixed << std::setprecision(3) << q.minAngleDeg;
            reason = oss.str();
            return true;
        }
        if (q.minAngleDeg > 0.0 && q.minAngleDeg < solidRemeshAngleTrigger) {
            std::ostringstream oss;
            oss << "solid_angle=" << std::fixed << std::setprecision(3) << q.minAngleDeg
                << "_edge=" << std::scientific << std::setprecision(6)
                << bStats.minMovingEdge;
            reason = oss.str();
            return true;
        }
        if (bStats.minMovingEdge > 0.0 &&
            bStats.minMovingEdge < solidRemeshMovingEdgeTrigger) {
            std::ostringstream oss;
            oss << "solid_edge=" << std::scientific << std::setprecision(6)
                << bStats.minMovingEdge
                << "_angle=" << std::fixed << std::setprecision(3) << q.minAngleDeg;
            reason = oss.str();
            return true;
        }
        return false;
    };

    auto remeshSolidAtCurrent = [&](double time, const std::string& reason) {
        SolidMeshQuality beforeQ = solid.currentMeshQuality();
        SolidBoundaryStats beforeB = movingBoundaryStats(solid);
        int oldNodes = solid.numNodes();
        int oldElems = solid.numElements();

        Mesh newCurrent = makeCurrentBlastRodSolidInteriorMesh(solid, solidRemeshH,
                                                               solidRemeshIter, false);
        solid.remeshToCurrentMesh(newCurrent);
        map.setSolid(&solid);
        aleReferenceNodes = solid.currentNodes();
        map.setReferenceNodes(aleReferenceNodes);
        map.setCurrent(time, solid.currentNodes(), solid.velocities());

        SolidMeshQuality referenceQ = solid.meshQuality();
        SolidMeshQuality afterQ = solid.currentMeshQuality();
        SolidBoundaryStats afterB = movingBoundaryStats(solid);
        ++solidRemeshCount;
        lastSolidRemeshStep = step;
        std::cout << "  solid_remesh#" << solidRemeshCount
                  << " at t=" << std::fixed << std::setprecision(5) << time
                  << " reason=" << reason
                  << " nodes=" << oldNodes << "->" << solid.numNodes()
                  << " elems=" << oldElems << "->" << solid.numElements()
                  << " ref_min_angle=" << std::setprecision(3)
                  << referenceQ.minAngleDeg
                  << " before_angle=" << beforeQ.minAngleDeg
                  << " after_angle=" << afterQ.minAngleDeg
                  << " before_inverted=" << beforeQ.invertedElements
                  << " after_inverted=" << afterQ.invertedElements
                  << " before_move_edge=" << std::scientific << std::setprecision(6)
                  << beforeB.minMovingEdge
                  << " after_move_edge=" << afterB.minMovingEdge
                  << " before_edge_mid=(" << std::fixed << std::setprecision(6)
                  << beforeB.midpoint.x() << "," << beforeB.midpoint.y() << ")"
                  << " after_edge_mid=(" << afterB.midpoint.x() << ","
                  << afterB.midpoint.y() << ")"
                  << std::defaultfloat << "\n";

        remeshFluidAtCurrent(time, "solid_remesh_" + reason, true);
    };

    while (t < tEnd - 1e-14) {
        setCurrentMap(t);
        updateStaticSpace(stepSpace, refMap, t);
        std::string solidRemeshReason;
        if (solidRemeshNeeded(solidRemeshReason)) {
            remeshSolidAtCurrent(t, solidRemeshReason);
        }
        std::string remeshReason;
        if (remeshNeeded("pre_step", remeshReason)) {
            remeshFluidAtCurrent(t, remeshReason);
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
        U = advanceOneStatic(stepSpace, nextSpace, refMap, t, dt, map, bc, U);
        t += dt;
        ++step;

        updateStaticSpace(stepSpace, refMap, t);
        applyRightSponge(U, sp, dt);
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
        remeshReason.clear();
        if (remeshNeeded("post_step", remeshReason)) {
            remeshFluidAtCurrent(t, remeshReason);
        }

        double rmin = U.col(0).minCoeff();
        double rmax = U.col(0).maxCoeff();
        diag << std::setprecision(12) << t << "," << solid.tipDisplacementX()
             << "," << solid.tipVelocityX() << ","
             << fFluid << "," << drag << "," << pMean << ","
             << sp.mesh.elem.rows() << "," << solid.numNodes() << ","
             << solid.numElements() << "," << sp.minH << "," << remeshCount
             << "," << rmin << "," << rmax << "\n";

        int frameWritten = -1;
        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            frameWritten = frame;
            writeFrame(frame++, t);
            nextFrame += frameDt;
        }

        SolidMeshQuality solidQStep = solid.currentMeshQuality();
        SolidBoundaryStats solidBStep = movingBoundaryStats(solid);
        std::ostringstream stepLog;
        stepLog << std::fixed << std::setprecision(6)
                << "  step_summary step=" << step
                << " t=" << t
                << " dt=" << dt
                << " frame=" << frameWritten
                << " tipX=" << solid.tipDisplacementX()
                << " tipV=" << solid.tipVelocityX()
                << " Fq=" << fFluid
                << " drag=" << drag
                << " pMean=" << pMean
                << " tris=" << sp.mesh.elem.rows()
                << " minH=" << sp.minH
                << " remesh=" << remeshCount
                << " solidRemesh=" << solidRemeshCount
                << " rho[" << rmin << "," << rmax << "]"
                << " solidMinAng=" << solidQStep.minAngleDeg
                << " solidInv=" << solidQStep.invertedElements
                << " solidMinMoveEdge=" << solidBStep.minMovingEdge
                << " solid_fem remeshable_mesh\n";
        std::cout << stepLog.str();
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
    const std::string cropFilter = "crop=" + std::to_string(videoCropW) + ":" +
                                   std::to_string(renderH) + ":0:0";
    runOutputCommand("ffmpeg -y -i " + still + " -frames:v 1 -update 1 " + stillPng);
    runOutputCommand("ffmpeg -y -i " + stillNoMesh + " -frames:v 1 -update 1 " + stillNoMeshPng);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(videoFps) + " -i " + dir +
                     "/frame_%05d.ppm -vf " + cropFilter +
                     " -c:v libx264 -pix_fmt yuv420p -crf 16 " + video);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(videoFps) + " -i " + dirNoMesh +
                     "/frame_%05d.ppm -vf " + cropFilter +
                     " -c:v libx264 -pix_fmt yuv420p -crf 16 " + videoNoMesh);
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    bool freshStart = false;
    bool remeshBench = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
        if (std::string(argv[i]) == "--fresh") freshStart = true;
        if (std::string(argv[i]) == "--remesh-bench") remeshBench = true;
    }
    if (remeshBench) return euler_ale::runBlastRodRemeshBench(quick);
    return euler_ale::runBlastRod(quick, freshStart);
}
