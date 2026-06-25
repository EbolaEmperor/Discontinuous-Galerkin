#include "BodyFittedMesh.h"
#include "SolidALEMap.h"
#include "Checkpoint.h"
#include "Core.h"
#include "DGState.h"
#include "Movie.h"
#include "NeoHookeanSolid.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace euler_ale {
namespace {

using namespace Eigen;

constexpr double PI = 3.141592653589793238462643383279502884;

struct DiskGeom {
    double xa = 0.0;
    double xb = 1.20;
    double ya = 0.0;
    double yb = 0.80;
    double cx = 0.40;
    double cy = 0.40;
    double radius = 0.125;
    double pinRadius = 0.017;
    double spongeFraction = 0.20;

    double spongeStartX() const {
        return xa + (1.0 - spongeFraction) * (xb - xa);
    }

    double spongeWidth() const {
        return std::max(0.0, xb - spongeStartX());
    }

    double spongeCoordinate(double x) const {
        double width = spongeWidth();
        if (width <= 1e-14) return 0.0;
        return std::clamp((x - spongeStartX()) / width, 0.0, 1.0);
    }
};

double smooth01(double s) {
    s = std::clamp(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
}

std::vector<Vector2d> diskBoundarySamples(const DiskGeom& g, int n) {
    std::vector<Vector2d> loop;
    n = std::max(32, n);
    loop.reserve(n);
    for (int i = 0; i < n; ++i) {
        double th = 2.0 * PI * static_cast<double>(i) / n;
        loop.emplace_back(g.cx + g.radius * std::cos(th),
                          g.cy + g.radius * std::sin(th));
    }
    return loop;
}

double diskSolidDistance(const DiskGeom& g, double x, double y) {
    return std::hypot(x - g.cx, y - g.cy) - g.radius;
}

void addOuterSegments(SolidBodyMeshSpec& spec, const DiskGeom& g) {
    spec.fixedSegments.push_back({Vector2d(g.xa, g.ya), Vector2d(g.xb, g.ya)});
    spec.fixedSegments.push_back({Vector2d(g.xb, g.ya), Vector2d(g.xb, g.yb)});
    spec.fixedSegments.push_back({Vector2d(g.xb, g.yb), Vector2d(g.xa, g.yb)});
    spec.fixedSegments.push_back({Vector2d(g.xa, g.yb), Vector2d(g.xa, g.ya)});
}

Mesh makeDiskSolidReferenceMesh(const DiskGeom& g, double h, int maxIter, bool verbose) {
    DistanceMeshSpec spec;
    spec.xa = g.cx - g.radius - 0.018;
    spec.xb = g.cx + g.radius + 0.018;
    spec.ya = g.cy - g.radius - 0.018;
    spec.yb = g.cy + g.radius + 0.018;
    spec.h0 = h;
    spec.seedH = 0.85 * h;
    spec.randomSeed = 20260623u;
    spec.signedDistance = [g](double x, double y) {
        return diskSolidDistance(g, x, y);
    };
    spec.targetSize = [h](double, double) { return h; };
    spec.fixedPoints = diskBoundarySamples(g, 96);
    spec.fixedPoints.emplace_back(g.cx, g.cy);
    for (int k = 0; k < 8; ++k) {
        double th = 2.0 * PI * static_cast<double>(k) / 8.0;
        spec.fixedPoints.emplace_back(g.cx + g.pinRadius * std::cos(th),
                                      g.cy + g.pinRadius * std::sin(th));
    }

    Mesh mesh;
    generateDistanceMesh(mesh, spec, maxIter, verbose);
    return mesh;
}

Mesh makeDiskSolidMeshFromCurrentBoundary(const DiskGeom& g,
                                          const ElasticSolid2D& solid,
                                          double h, int maxIter, bool verbose) {
    std::vector<Vector2d> loop = solidBoundaryLoop(solid);
    if (loop.size() < 8) return makeDiskSolidReferenceMesh(g, h, maxIter, verbose);

    Vector2d lo = loop.front();
    Vector2d hi = loop.front();
    for (const auto& p : loop) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    double pad = 0.02;
    DistanceMeshSpec spec;
    spec.xa = lo.x() - pad;
    spec.xb = hi.x() + pad;
    spec.ya = lo.y() - pad;
    spec.yb = hi.y() + pad;
    spec.h0 = h;
    spec.seedH = 0.85 * h;
    spec.randomSeed = 20260623u;
    spec.signedDistance = [loop](double x, double y) {
        return signedDistancePolygon(loop, x, y);
    };
    spec.targetSize = [h](double, double) { return h; };
    spec.fixedPoints = loop;
    spec.fixedPoints.emplace_back(g.cx, g.cy);
    for (int k = 0; k < 8; ++k) {
        double th = 2.0 * PI * static_cast<double>(k) / 8.0;
        spec.fixedPoints.emplace_back(g.cx + g.pinRadius * std::cos(th),
                                      g.cy + g.pinRadius * std::sin(th));
    }

    Mesh mesh;
    generateDistanceMesh(mesh, spec, maxIter, verbose);
    return mesh;
}

Mesh makeDiskFluidMesh(const DiskGeom& g, const ElasticSolid2D& solid,
                       double h, int maxIter, bool verbose) {
    std::vector<Vector2d> loop = solidBoundaryLoop(solid);
    SolidBodyMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.hFar = h;
    spec.hNearFactor = 0.22;
    spec.gradeRadius = 0.38;
    spec.randomSeed = 20260624u;
    spec.maxIter = maxIter;
    spec.verbose = verbose;
    spec.solidDistance = [loop](double x, double y) {
        return signedDistancePolygon(loop, x, y);
    };
    spec.solid = &solid;
    spec.addMovingSolidBoundary = true;
    addOuterSegments(spec, g);
    return makeSolidBodyFittedMesh(spec);
}

int diskBoundaryTag(double x, double y, double time,
                    const SolidALEMap& map, const DiskGeom& g) {
    double dLeft = std::abs(x - g.xa);
    double dRight = std::abs(x - g.xb);
    double dBottom = std::abs(y - g.ya);
    double dTop = std::abs(y - g.yb);
    double dSolid = map.distanceToBoundary(x, y, time);

    double best = dLeft;
    int tag = TAG_EXACT;
    if (dRight < best) {
        best = dRight;
        tag = TAG_OUTFLOW;
    }
    if (dBottom < best) {
        best = dBottom;
        tag = TAG_SLIP_WALL;
    }
    if (dTop < best) {
        best = dTop;
        tag = TAG_SLIP_WALL;
    }
    if (dSolid < best) tag = TAG_MOVING_WALL;
    return tag;
}

Vector2d movingCentroid(const ElasticSolid2D& solid) {
    const MatrixXd& x = solid.currentNodes();
    const VectorXi& fixed = solid.fixedMask();
    Vector2d c = Vector2d::Zero();
    int n = 0;
    for (int i = 0; i < x.rows(); ++i) {
        if (fixed(i)) continue;
        c += x.row(i).transpose();
        ++n;
    }
    if (n == 0) return Vector2d::Zero();
    return c / static_cast<double>(n);
}

double maxDisplacement(const ElasticSolid2D& solid) {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    double best = 0.0;
    for (int i = 0; i < x.rows(); ++i) {
        best = std::max(best, (x.row(i) - X.row(i)).norm());
    }
    return best;
}

double extremeReferenceDisplacement(const ElasticSolid2D& solid,
                                    int referenceAxis, bool maximumSide,
                                    int displacementComponent) {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    int best = 0;
    for (int i = 1; i < X.rows(); ++i) {
        if (maximumSide) {
            if (X(i, referenceAxis) > X(best, referenceAxis)) best = i;
        } else {
            if (X(i, referenceAxis) < X(best, referenceAxis)) best = i;
        }
    }
    return x(best, displacementComponent) - X(best, displacementComponent);
}

double maxFixedDisplacement(const ElasticSolid2D& solid) {
    const MatrixXd& X = solid.referenceNodes();
    const MatrixXd& x = solid.currentNodes();
    const VectorXi& fixed = solid.fixedMask();
    double best = 0.0;
    for (int i = 0; i < x.rows(); ++i) {
        if (!fixed(i)) continue;
        best = std::max(best, (x.row(i) - X.row(i)).norm());
    }
    return best;
}

void copyFrameStill(const std::string& framePath, const std::string& stillPath) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::copy_file(framePath, stillPath, fs::copy_options::overwrite_existing, ec);
    if (ec) {
        std::cerr << "Warning: cannot write still " << stillPath << ": "
                  << ec.message() << "\n";
    }
}

std::string frameName(const std::string& folder, int frame) {
    char fn[512];
    std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", folder.c_str(), frame);
    return std::string(fn);
}

} // namespace

int runShockDisk(bool quick, bool freshStart = false) {
    namespace fs = std::filesystem;

    const int ord = 1;
    const int nFrames = quick ? 28 : 1800;
    const double tEnd = quick ? 0.22 : 7.20;
    const double hFluid = quick ? 0.052 : 0.021;
    const double hSolid = quick ? 0.021 : 0.010;
    const int fluidIter = quick ? 20 : 30;
    const int solidIter = quick ? 25 : 35;
    const double cfl = quick ? 0.16 : 0.13;
    const double rhoFloor = 0.08;
    const double pFloor = 0.06;
    const double speedMax = 4.5;
    const double pExt = 1.0;
    const double remeshFluidMinAngle = quick ? 4.0 : 5.0;
    const double remeshSolidMinAngle = quick ? 6.0 : 8.0;
    const int forceSmoothPasses = quick ? 2 : 4;
    const double forceSmoothBlend = 0.55;

    DiskGeom geom;
    SolidMaterial material;
    material.density = 8.0;
    material.thickness = 1.0;
    material.young = 12.0;
    material.poisson = 0.33;
    material.damping = 0.35;
    NeoHookeanSolidModel solidModel(NeoHookeanMaterial::fromSolidMaterial(material));

    std::cout << "Shock-loaded centre-pinned elastic disk: Euler ALE + NeoHookeanSolidModel\n";
    std::cout << "  domain=[" << geom.xa << "," << geom.xb << "]x["
              << geom.ya << "," << geom.yb << "] aspect="
              << (geom.xb - geom.xa) / (geom.yb - geom.ya)
              << " sponge=[" << geom.spongeStartX() << "," << geom.xb << "]\n";
    std::cout << "  generating solid FEM mesh h=" << hSolid << "...\n";
    Mesh solidMesh = makeDiskSolidReferenceMesh(geom, hSolid, solidIter, true);
    ElasticSolid2D solid;
    solid.resetReferenceMesh(solidMesh, material);
    solid.setFixedNodesInDisk(Vector2d(geom.cx, geom.cy), geom.pinRadius, true);
    solid.setAllBoundarySegmentsMoving();
    SolidMeshQuality solidQ = solid.meshQuality();
    std::cout << "  solid nodes=" << solid.numNodes()
              << " elems=" << solid.numElements()
              << " min_angle=" << solidQ.minAngleDeg
              << " fixed_center_radius=" << geom.pinRadius
              << " mass=" << solid.totalMass()
              << " force_smooth_passes=" << forceSmoothPasses
              << " force_smooth_blend=" << forceSmoothBlend << "\n";

    SolidALEMap map;
    map.setSolid(&solid);
    map.setDomain(geom.xa, geom.xb, geom.ya, geom.yb);
    map.setInfluence(0.16, 0.12);
    map.setCurrent(0.0, solid.currentNodes(), solid.velocities());

    std::cout << "  generating body-fitted fluid mesh h=" << hFluid << "...\n";
    Mesh base = makeDiskFluidMesh(geom, solid, hFluid, fluidIter, true);
    double minAng = 0.0;
    double meanAng = 0.0;
    double minArea = 0.0;
    double maxArea = 0.0;
    meshQuality(base, minAng, meanAng, minArea, maxArea);
    std::cout << "  fluid nodes=" << base.node.rows()
              << " elems=" << base.elem.rows()
              << " min_angle=" << minAng
              << " mean_angle=" << meanAng << "\n";

    RefMapFn refMap = [&](const Vector2d& X, double time) {
        return map.refToPhys(X, time);
    };
    MaxMeshSpeedFn maxSpeed = [&](double time) {
        return map.maxMeshSpeed(time);
    };
    MeshVelocityFn meshVel = [&](double x, double y, double time) {
        return map.velocityAt(x, y, time);
    };
    Tagger tagger = [&](double x, double y, double time) {
        return diskBoundaryTag(x, y, time, map, geom);
    };
    auto driverPrim = [](double time) {
        double period = 1.3;
        double onDuration = 0.8;
        double phase = std::fmod(time, period);
        double active = (phase < onDuration) ? 1.0 : 0.0;
        double rampUp = smooth01(std::min(phase, 0.030) / 0.030);
        double rampDown = smooth01(std::clamp((onDuration - phase) / 0.030, 0.0, 1.0));
        double envelope = active * std::min(rampUp, rampDown);
        int pulseIndex = static_cast<int>(std::floor(time / period));
        double intensity = 1.0 + 0.5 * std::max(0, pulseIndex - 2);
        double rho = 1.0 + 0.80 * intensity * envelope;
        double u = 0.98 * intensity * envelope;
        double p = 1.0 + 2.35 * intensity * envelope;
        return Vector4d(rho, u, 0.0, p);
    };
    ALEBCFn bc = [&](double, double, double time, const Vector4d& Um,
                     double nx, double ny, int tag, double wn) {
        if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL) {
            return movingWallGhost(Um, nx, ny, (tag == TAG_MOVING_WALL) ? wn : 0.0);
        }
        if (tag == TAG_EXACT) {
            Vector4d pr = driverPrim(time);
            return euler::primToCons(pr(0), pr(1), pr(2), pr(3));
        }
        if (tag == TAG_OUTFLOW) {
            return characteristicPressureOutletGhost(Um, nx, ny, wn,
                                                     Vector4d(1.0, 0.0, 0.0, pExt));
        }
        return Um;
    };

    ALEAdaptiveForest forest(base, ord, 4);
    Space sp;
    rebuildSpace(forest, ord, refMap, 0.0, tagger, sp);
    MatrixXd U = euler::projectInitial(*sp.fem, sp.mesh, sp.e2d,
        [&](double x, double) {
            if (x < 0.16) return Vector4d(1.74, 0.76, 0.0, 3.22);
            return Vector4d(1.0, 0.0, 0.0, 1.0);
        });
    applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
    const Vector4d spongeReference = euler::primToCons(1.0, 0.0, 0.0, pExt);
    const double spongeSigmaMax =
        3.0 * (0.98 + std::sqrt(euler::GAMMA * pExt)) /
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

    double t = 0.0;
    double nextFrame = 0.0;
    double frameDt = tEnd / std::max(1, nFrames);
    int step = 0;
    int frame = 0;
    int remeshCount = 0;

    const std::string checkpointPrefix = "shock_disk";
    std::string prefix = quick ? "shock_disk_quick" : "shock_disk";
    std::string dir = "out/" + prefix + "_frames";
    std::string dirMesh = "out/" + prefix + "_mesh_frames";
    std::string dirSch = "out/" + prefix + "_schlieren_frames";
    std::string csvPath = "out/" + prefix + "_diagnostics.csv";
    fs::create_directories("out");
    fs::create_directories(dir);
    fs::create_directories(dirMesh);
    fs::create_directories(dirSch);

    std::vector<CheckpointMilestone> checkpointPlan = checkpointSchedule(quick);
    std::vector<int> checkpointDone(checkpointPlan.size(), 0);
    bool resumed = false;

    if (!freshStart) {
        std::optional<RunCheckpoint> resumeCP =
            loadLatestCheckpoint(checkpointPrefix, quick, ord, nFrames, tEnd,
                                 hFluid, solid.numNodes(), true);
        if (resumeCP.has_value()) {
            const RunCheckpoint& cp = *resumeCP;
            if (cp.solidReferenceMesh.node.rows() > 0) {
                solid.resetReferenceMesh(cp.solidReferenceMesh, material);
                solid.setFixedNodesInDisk(Vector2d(geom.cx, geom.cy), geom.pinRadius, true);
                solid.setAllBoundarySegmentsMoving();
                map.setSolid(&solid);
            }
            solid.setState(cp.solidNodes, cp.solidVelocities);
            map.setCurrent(cp.time, solid.currentNodes(), solid.velocities());
            forest = ALEAdaptiveForest(cp.referenceMesh, ord, 4);
            rebuildSpace(forest, ord, refMap, cp.time, tagger, sp);
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
                pruneFramesFrom(dirMesh, frame);
                pruneFramesFrom(dirSch, frame);
                trimDiagnosticsToTime(csvPath, t);
                resumed = true;
                std::cout << "  resumed from checkpoint at t=" << std::fixed
                          << std::setprecision(5) << t
                          << " step=" << step
                          << " frame=" << frame
                          << " remesh=" << remeshCount << "\n";
            } else {
                std::cerr << "Warning: checkpoint DG state incompatible; cold start\n";
            }
        }
    } else {
        std::cout << "  fresh start; checkpoint auto-resume disabled\n";
        if (!quick) {
            pruneOldCheckpoints(checkpointPrefix, quick, 0);
        }
    }

    if (!resumed) {
        clearFrameDirectory(dir);
        clearFrameDirectory(dirMesh);
        clearFrameDirectory(dirSch);
    }

    std::ofstream diag;
    if (resumed) {
        diag.open(csvPath, std::ios::app);
    } else {
        diag.open(csvPath);
        diag << "time,centroid_x,centroid_y,left_dx,right_dx,top_dy,bottom_dy,"
                "max_displacement,pin_max_displacement,max_speed,"
                "drag,lift,mean_pressure,fluid_triangles,solid_nodes,solid_elements,"
                "solid_min_angle,solid_inverted,rho_min,rho_max,min_h,strain_energy,"
                "kinetic_energy\n";
    }

    const int W = quick ? 960 : 1440;
    const int H = quick ? 640 : 960;
    const int ssaa = quick ? 1 : 2;
    const double viewXa = geom.xa;
    const double viewXb = geom.xb;
    const double viewYa = geom.ya;
    const double viewYb = geom.yb;

    auto writeFrame = [&](int idx, double time) {
        rebuildSpace(forest, ord, refMap, time, tagger, sp);
        int hiW = W * ssaa;
        int hiH = H * ssaa;
        std::vector<unsigned char> hi =
            euler::renderScalarPPMImage(*sp.fem, sp.mesh, sp.e2d, U.col(0),
                                        hiW, hiH, viewXa, viewXb, viewYa, viewYb,
                                        0.55, 3.25, euler::CM_INFERNO);
        overlaySolidMesh(hi, hiW, hiH, solid, viewXa, viewXb, viewYa, viewYb, true);
        int outW = hiW;
        int outH = hiH;
        std::vector<unsigned char> img = hi;
        if (ssaa > 1) img = downsampleImage(hi, hiW, hiH, ssaa, outW, outH);
        writePPM(frameName(dir, idx), outW, outH, img);

        std::vector<unsigned char> meshImg = img;
        overlayMesh(meshImg, outW, outH, sp.mesh, viewXa, viewXb, viewYa, viewYb);
        overlaySolidMesh(meshImg, outW, outH, solid, viewXa, viewXb, viewYa, viewYb, true);
        writePPM(frameName(dirMesh, idx), outW, outH, meshImg);

        std::string schPath = frameName(dirSch, idx);
        euler::writeSchlierenPPM(schPath, *sp.fem, sp.mesh, sp.e2d, U.col(0),
                                 outW, outH, viewXa, viewXb, viewYa, viewYb, 9.0);
        overlaySolidMesh(schPath, solid, viewXa, viewXb, viewYa, viewYb, true);
    };

    if (!resumed) {
        writeFrame(frame++, t);
        nextFrame += frameDt;
    }

    while (t < tEnd - 1e-14) {
        map.setCurrent(t, solid.currentNodes(), solid.velocities());
        rebuildSpace(forest, ord, refMap, t, tagger, sp);

        double pMean = pExt;
        double drag = 0.0;
        double lift = 0.0;
        solid.clearExternalForces();
        loadSolidFromFluidPressure(sp, U, solid, pExt, &pMean, &drag, &lift);
        solid.smoothMovingBoundaryForces(forceSmoothPasses, forceSmoothBlend);

        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        double dtSolid = solidModel.stableTimeStep(solid, 0.30);
        dt = std::min(dt, dtSolid);
        dt = std::min(dt, quick ? 0.00115 : 0.00085);
        if (!(dt > 1e-15)) {
            std::cerr << "DT-DEBUG t=" << t << " step=" << step
                      << " estDt=" << estimateDt(sp, U, maxSpeed, t, ord, cfl)
                      << " solidDt=" << dtSolid
                      << " minH=" << sp.minH
                      << " elems=" << sp.mesh.elem.rows()
                      << " Urows=" << U.rows() << " nDof=" << sp.nDof
                      << " rho_min=" << U.col(0).minCoeff()
                      << " rho_max=" << U.col(0).maxCoeff() << std::endl;
            break;
        }

        MatrixXd solidX0 = solid.currentNodes();
        MatrixXd solidV0 = solid.velocities();
        solidModel.advanceExplicit(solid, dt);
        MatrixXd solidX1 = solid.currentNodes();
        MatrixXd solidV1 = solid.velocities();
        map.setMotion(t, t + dt, solidX0, solidV0, solidX1, solidV1);
        U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
        t += dt;
        ++step;

        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        applyRightSponge(U, sp, dt);
        applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);

        bool needRemesh = false;
        SolidMeshQuality solidQ2 = solid.currentMeshQuality();
        if (solidQ2.invertedElements > 0 || solidQ2.minAngleDeg < remeshSolidMinAngle) {
            needRemesh = true;
        }
        if (!needRemesh) {
            double fluidMinAng = 0.0, fluidMeanAng = 0.0, fluidMinA = 0.0, fluidMaxA = 0.0;
            meshQuality(sp.mesh, fluidMinAng, fluidMeanAng, fluidMinA, fluidMaxA);
            if (fluidMinAng < remeshFluidMinAngle) needRemesh = true;
        }

        if (needRemesh && t < tEnd - 1e-14) {
            double fluidMinAngBefore = 0.0, fluidMeanAngBefore = 0.0, fluidMinABefore = 0.0, fluidMaxABefore = 0.0;
            meshQuality(sp.mesh, fluidMinAngBefore, fluidMeanAngBefore, fluidMinABefore, fluidMaxABefore);
            std::cout << "  coupled remesh at t=" << std::fixed << std::setprecision(4) << t
                      << " solid_min_angle=" << solidQ2.minAngleDeg
                      << " solid_inverted=" << solidQ2.invertedElements
                      << " fluid_min_angle=" << fluidMinAngBefore << "\n";
            Mesh newSolidMesh = makeDiskSolidMeshFromCurrentBoundary(geom, solid, hSolid, solidIter, false);
            solid.remeshToCurrentMesh(newSolidMesh);
            solid.setFixedNodesInDisk(Vector2d(geom.cx, geom.cy), geom.pinRadius, true);
            solid.setAllBoundarySegmentsMoving();
            map.setSolid(&solid);
            map.setReferenceNodes(solid.currentNodes());
            map.setCurrent(t, solid.currentNodes(), solid.velocities());

            Space oldSp = std::move(sp);
            MatrixXd oldU = U;
            Mesh newFluidBase = makeDiskFluidMesh(geom, solid, hFluid, fluidIter, false);
            base = newFluidBase;
            forest = ALEAdaptiveForest(base, ord, 4);
            rebuildSpace(forest, ord, refMap, t, tagger, sp);
            U = interpolateDGToSpace(oldSp, oldU, sp);
            applyPrimitiveBounds(U, rhoFloor, pFloor, speedMax);
            ++remeshCount;
            double fluidMinAngAfter = 0.0, fluidMeanAngAfter = 0.0, fluidMinAAfter = 0.0, fluidMaxAAfter = 0.0;
            meshQuality(sp.mesh, fluidMinAngAfter, fluidMeanAngAfter, fluidMinAAfter, fluidMaxAAfter);
            std::cout << "    new fluid mesh: nodes=" << sp.mesh.node.rows()
                      << " elems=" << sp.mesh.elem.rows()
                      << " fluid_min_angle=" << fluidMinAngAfter
                      << " solid_min_angle=" << solid.currentMeshQuality().minAngleDeg
                      << " solid_dt=" << solidModel.stableTimeStep(solid, 0.30) << "\n";
        }

        double rmin = U.col(0).minCoeff();
        double rmax = U.col(0).maxCoeff();
        SolidMeshQuality sq = solid.currentMeshQuality();
        Vector2d c = movingCentroid(solid);
        diag << std::setprecision(12) << t << "," << c.x() << "," << c.y()
             << "," << extremeReferenceDisplacement(solid, 0, false, 0)
             << "," << extremeReferenceDisplacement(solid, 0, true, 0)
             << "," << extremeReferenceDisplacement(solid, 1, true, 1)
             << "," << extremeReferenceDisplacement(solid, 1, false, 1)
             << "," << maxDisplacement(solid)
             << "," << maxFixedDisplacement(solid)
             << "," << solid.maxNodeSpeed()
             << "," << drag << "," << lift << "," << pMean
             << "," << sp.mesh.elem.rows()
             << "," << solid.numNodes()
             << "," << solid.numElements()
             << "," << sq.minAngleDeg
             << "," << sq.invertedElements
             << "," << rmin << "," << rmax
             << "," << sp.minH
             << "," << solidModel.strainEnergy(solid)
             << "," << solidModel.kineticEnergy(solid) << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step
                      << " leftDx=" << extremeReferenceDisplacement(solid, 0, false, 0)
                      << " rightDx=" << extremeReferenceDisplacement(solid, 0, true, 0)
                      << " maxDisp=" << maxDisplacement(solid)
                      << " pinDisp=" << maxFixedDisplacement(solid)
                      << " speed=" << solid.maxNodeSpeed()
                      << " F=(" << drag << "," << lift << ")"
                      << " pMean=" << pMean
                      << " fluidTris=" << sp.mesh.elem.rows()
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " solidMinAng=" << sq.minAngleDeg
                      << " inv=" << sq.invertedElements
                      << " remesh=" << remeshCount << "\n";
            nextFrame += frameDt;
        }

        double frac = (tEnd > 0.0) ? t / tEnd : 1.0;
        for (int i = 0; i < static_cast<int>(checkpointPlan.size()); ++i) {
            if (checkpointDone[i]) continue;
            if (frac + 1e-12 < checkpointPlan[i].fraction) continue;
            checkpointDone[i] = 1;
            diag.flush();
            RunCheckpoint cp;
            cp.quick = quick;
            cp.ord = ord;
            cp.nFrames = nFrames;
            cp.tEnd = tEnd;
            cp.h = hFluid;
            cp.time = t;
            cp.nextFrame = nextFrame;
            cp.step = step;
            cp.frame = frame;
            cp.remeshCount = remeshCount;
            cp.milestoneDone = checkpointDone;
            cp.referenceMesh = base;
            cp.solidReferenceMesh.node = solid.referenceNodes();
            cp.solidReferenceMesh.elem = solid.elements();
            cp.U = U;
            cp.solidNodes = solid.currentNodes();
            cp.solidVelocities = solid.velocities();
            cp.aleReferenceNodes = solid.referenceNodes();
            fs::path cpPath = checkpointPath(checkpointPrefix, quick,
                                             checkpointPlan[i].label);
            if (writeCheckpointAtomic(cpPath, cp)) {
                if (!quick) pruneOldCheckpoints(checkpointPrefix, quick, 3);
                std::cout << "  checkpoint " << checkpointPlan[i].label
                          << "% written at t=" << std::fixed << std::setprecision(5)
                          << t << " -> " << cpPath.string() << "\n";
            }
        }
    }

    std::string still = "out/" + prefix + ".ppm";
    std::string stillPng = "out/" + prefix + ".png";
    std::string stillMesh = "out/" + prefix + "_mesh.ppm";
    std::string stillMeshPng = "out/" + prefix + "_mesh.png";
    std::string stillSch = "out/" + prefix + "_schlieren.ppm";
    std::string stillSchPng = "out/" + prefix + "_schlieren.png";
    copyFrameStill(frameName(dir, frame - 1), still);
    copyFrameStill(frameName(dirMesh, frame - 1), stillMesh);
    copyFrameStill(frameName(dirSch, frame - 1), stillSch);

    int fps = quick ? 30 : 60;
    std::string video = "out/" + prefix + ".mp4";
    std::string videoMesh = "out/" + prefix + "_mesh.mp4";
    std::string videoSch = "out/" + prefix + "_schlieren.mp4";
    runOutputCommand("ffmpeg -y -i " + still + " -frames:v 1 -update 1 " + stillPng);
    runOutputCommand("ffmpeg -y -i " + stillMesh + " -frames:v 1 -update 1 " + stillMeshPng);
    runOutputCommand("ffmpeg -y -i " + stillSch + " -frames:v 1 -update 1 " + stillSchPng);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(fps) + " -i " + dir +
                     "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 " + video);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(fps) + " -i " + dirMesh +
                     "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 17 " + videoMesh);
    runOutputCommand("ffmpeg -y -framerate " + std::to_string(fps) + " -i " + dirSch +
                     "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 " + videoSch);

    std::cout << "Done. frames=" << frame
              << " density=" << video
              << " mesh=" << videoMesh
              << " schlieren=" << videoSch
              << " diagnostics=" << csvPath << "\n";
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    bool fresh = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
        if (std::string(argv[i]) == "--fresh") fresh = true;
    }
    return euler_ale::runShockDisk(quick, fresh);
}
