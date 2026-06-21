#include "CGCase.h"

#include "BodyFittedMesh.h"
#include "Core.h"
#include "Domain.h"
#include "ElasticSolid.h"
#include "Output.h"

#include "DG.h"
#include "FEM.h"
#include "MeshGen.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace euler_ale {
namespace {

struct CGSpace {
    Mesh mesh;
    MatrixXi edge;
    MatrixXi edge2side;
    VectorXi tag;
    VectorXd mass;
    double minH = 0.0;
};

int fullDomainRenderWidth(int physicalWidth, const BlastRodGeom& geom) {
    double physicalLength = std::max(geom.flowXb() - geom.xa, 1e-14);
    double computationalLength = std::max(geom.xb - geom.xa, physicalLength);
    int w = static_cast<int>(std::lround(physicalWidth * computationalLength / physicalLength));
    w = std::max(physicalWidth, w);
    if (w % 2) --w;
    return std::max(2, w);
}

Vector4d inflowPrimitive(double time) {
    double r = blastRodSmoothRamp(time / 0.035);
    double rho = 1.0 + 0.72 * r;
    double u = 0.90 * r;
    double p = 1.0 + 2.35 * r;
    return Vector4d(rho, u, 0.0, p);
}

Vector4d initialPrimitive(const BlastRodGeom& geom, double x) {
    if (x < geom.rodLeft() - 0.055) return Vector4d(1.72, 0.10, 0.0, 3.35);
    return Vector4d(1.0, 0.0, 0.0, 1.0);
}

void applyNodeBounds(MatrixXd& U, double rhoFloor, double pFloor, double speedMax) {
    for (int i = 0; i < U.rows(); ++i) {
        Vector4d Ui = U.row(i).transpose();
        double rho = std::max(Ui(0), rhoFloor);
        double u = Ui(1) / std::max(Ui(0), rhoFloor);
        double v = Ui(2) / std::max(Ui(0), rhoFloor);
        double speed = std::hypot(u, v);
        if (speed > speedMax) {
            double s = speedMax / std::max(speed, 1e-14);
            u *= s;
            v *= s;
        }
        double p = std::max(euler::pressure(Ui, rhoFloor, pFloor), pFloor);
        U.row(i) = euler::primToCons(rho, u, v, p).transpose();
    }
}

Mesh mappedMesh(const Mesh& referenceMesh, const RefMapFn& refMap, double time) {
    Mesh mesh = referenceMesh;
    for (int i = 0; i < mesh.node.rows(); ++i) {
        Vector2d X = referenceMesh.node.row(i).transpose();
        mesh.node.row(i) = refMap(X, time).transpose();
    }
    return mesh;
}

CGSpace buildCGSpace(const Mesh& referenceMesh, const RefMapFn& refMap,
                     double time, const SolidALEMap& map,
                     const BlastRodGeom& geom) {
    CGSpace sp;
    sp.mesh = mappedMesh(referenceMesh, refMap, time);
    sp.mesh.getEdge2Side(sp.edge, sp.edge2side);
    sp.tag = VectorXi::Constant(sp.edge.rows(), TAG_INTERIOR);
    sp.mass = VectorXd::Zero(sp.mesh.node.rows());
    sp.minH = 1e300;

    for (int e = 0; e < sp.edge.rows(); ++e) {
        bool boundary = (sp.edge2side(e, 0) < 0 || sp.edge2side(e, 1) < 0);
        if (!boundary) continue;
        Vector2d mid = 0.5 * (sp.mesh.node.row(sp.edge(e, 0)).transpose() +
                              sp.mesh.node.row(sp.edge(e, 1)).transpose());
        sp.tag(e) = blastRodBoundaryTag(mid.x(), mid.y(), time, map, geom);
    }

    for (int t = 0; t < sp.mesh.elem.rows(); ++t) {
        double area = triArea(sp.mesh, t);
        sp.minH = std::min(sp.minH, hCFL(sp.mesh, t));
        for (int k = 0; k < 3; ++k) sp.mass(sp.mesh.elem(t, k)) += area / 3.0;
    }
    if (!std::isfinite(sp.minH)) sp.minH = 0.0;
    return sp;
}

std::array<Vector2d, 3> p1Gradients(const Mesh& mesh, int elem) {
    Vector2d p0 = mesh.node.row(mesh.elem(elem, 0)).transpose();
    Vector2d p1 = mesh.node.row(mesh.elem(elem, 1)).transpose();
    Vector2d p2 = mesh.node.row(mesh.elem(elem, 2)).transpose();
    double twoA = (p1.x() - p0.x()) * (p2.y() - p0.y()) -
                  (p1.y() - p0.y()) * (p2.x() - p0.x());
    if (std::abs(twoA) < 1e-300) twoA = (twoA < 0.0) ? -1e-300 : 1e-300;
    return {
        Vector2d((p1.y() - p2.y()) / twoA, (p2.x() - p1.x()) / twoA),
        Vector2d((p2.y() - p0.y()) / twoA, (p0.x() - p2.x()) / twoA),
        Vector2d((p0.y() - p1.y()) / twoA, (p1.x() - p0.x()) / twoA)
    };
}

double nodeSignalSpeed(const Vector4d& U) {
    double rho = std::max(U(0), 1e-14);
    double u = U(1) / rho;
    double v = U(2) / rho;
    return std::hypot(u, v) + euler::soundSpeed(U);
}

double estimateCGDt(const CGSpace& sp, const MatrixXd& U,
                    const MaxMeshSpeedFn& maxMeshSpeed, double time, double cfl) {
    double maxMesh = maxMeshSpeed(time);
    double dt = 1e300;
    for (int t = 0; t < sp.mesh.elem.rows(); ++t) {
        double a = maxMesh;
        for (int k = 0; k < 3; ++k) {
            int n = sp.mesh.elem(t, k);
            a = std::max(a, nodeSignalSpeed(U.row(n).transpose()) + maxMesh);
        }
        dt = std::min(dt, cfl * hCFL(sp.mesh, t) / std::max(a, 1e-12));
    }
    return std::isfinite(dt) ? dt : 1e-5;
}

Vector4d boundaryGhost(const Vector4d& Um, double nx, double ny, double wn,
                       int tag, double time) {
    if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL)
        return movingWallGhost(Um, nx, ny, (tag == TAG_MOVING_WALL) ? wn : 0.0);
    if (tag == TAG_OUTFLOW)
        return characteristicPressureOutletGhost(Um, nx, ny, wn,
                                                 Vector4d(1.0, 0.0, 0.0, 1.0));
    if (tag == TAG_EXACT) {
        Vector4d pr = inflowPrimitive(time);
        return euler::primToCons(pr(0), pr(1), pr(2), pr(3));
    }
    return Um;
}

MatrixXd cgResidual(const CGSpace& sp, const MatrixXd& U,
                    const MeshVelocityFn& meshVelocity, double time,
                    double graphViscosity) {
    MatrixXd R = MatrixXd::Zero(U.rows(), 4);

    for (int t = 0; t < sp.mesh.elem.rows(); ++t) {
        double area = triArea(sp.mesh, t);
        auto grad = p1Gradients(sp.mesh, t);
        Vector4d Uc = Vector4d::Zero();
        Vector2d xc = Vector2d::Zero();
        for (int k = 0; k < 3; ++k) {
            int n = sp.mesh.elem(t, k);
            Uc += U.row(n).transpose();
            xc += sp.mesh.node.row(n).transpose();
        }
        Uc /= 3.0;
        xc /= 3.0;
        Vector4d Fx, Fy;
        euler::fluxes(Uc, Fx, Fy);
        Vector2d w = meshVelocity(xc.x(), xc.y(), time);
        Fx -= w.x() * Uc;
        Fy -= w.y() * Uc;
        for (int k = 0; k < 3; ++k) {
            int n = sp.mesh.elem(t, k);
            Vector4d r = area * (grad[k].x() * Fx + grad[k].y() * Fy);
            R.row(n) += r.transpose();
        }
    }

    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) == TAG_INTERIOR) continue;
        int tt = (sp.edge2side(e, 0) >= 0) ? sp.edge2side(e, 0) : sp.edge2side(e, 1);
        if (tt < 0) continue;
        int n1 = sp.edge(e, 0);
        int n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, tt, n1, n2);
        Vector4d Um = 0.5 * (U.row(n1).transpose() + U.row(n2).transpose());
        Vector2d mid = 0.5 * (sp.mesh.node.row(n1).transpose() +
                              sp.mesh.node.row(n2).transpose());
        double wn = meshVelocity(mid.x(), mid.y(), time).dot(eo.nout);
        Vector4d Up = boundaryGhost(Um, eo.nout.x(), eo.nout.y(), wn, sp.tag(e), time);
        Vector4d flux = aleRusanov(Um, Up, eo.nout.x(), eo.nout.y(), wn);
        R.row(n1) -= (0.5 * eo.he * flux).transpose();
        R.row(n2) -= (0.5 * eo.he * flux).transpose();
    }

    for (int e = 0; e < sp.edge.rows(); ++e) {
        int n1 = sp.edge(e, 0);
        int n2 = sp.edge(e, 1);
        Vector4d U1 = U.row(n1).transpose();
        Vector4d U2 = U.row(n2).transpose();
        double len = (sp.mesh.node.row(n1) - sp.mesh.node.row(n2)).norm();
        double a = std::max(nodeSignalSpeed(U1), nodeSignalSpeed(U2));
        double coeff = graphViscosity * len * a;
        R.row(n1) += (coeff * (U2 - U1)).transpose();
        R.row(n2) += (coeff * (U1 - U2)).transpose();
    }

    return R;
}

MatrixXd advanceCG(const CGSpace& sp, const MatrixXd& U,
                   const MeshVelocityFn& meshVelocity, double time, double dt,
                   double graphViscosity, double rhoFloor, double pFloor,
                   double speedMax) {
    MatrixXd R0 = cgResidual(sp, U, meshVelocity, time, graphViscosity);
    MatrixXd U1 = U;
    for (int i = 0; i < U.rows(); ++i) U1.row(i) += (dt / sp.mass(i)) * R0.row(i);
    applyNodeBounds(U1, rhoFloor, pFloor, speedMax);

    MatrixXd R1 = cgResidual(sp, U1, meshVelocity, time + dt, graphViscosity);
    MatrixXd Un = U;
    for (int i = 0; i < U.rows(); ++i) {
        Un.row(i) = 0.5 * (U.row(i) + U1.row(i) + (dt / sp.mass(i)) * R1.row(i));
    }
    applyNodeBounds(Un, rhoFloor, pFloor, speedMax);
    return Un;
}

double loadSolidFromCG(const CGSpace& sp, const MatrixXd& U, ElasticSolid2D& solid,
                       double pExt, double* meanPressure, double* drag) {
    double fx = 0.0;
    double area = 0.0;
    double pInt = 0.0;
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) != TAG_MOVING_WALL) continue;
        int tt = (sp.edge2side(e, 0) >= 0) ? sp.edge2side(e, 0) : sp.edge2side(e, 1);
        if (tt < 0) continue;
        int n1 = sp.edge(e, 0);
        int n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, tt, n1, n2);
        Vector4d Um = 0.5 * (U.row(n1).transpose() + U.row(n2).transpose());
        double p = euler::pressure(Um, 1e-8, 1e-8);
        Vector2d mid = 0.5 * (sp.mesh.node.row(n1).transpose() +
                              sp.mesh.node.row(n2).transpose());
        Vector2d traction = (p - pExt) * eo.nout;
        solid.addBoundaryTractionAt(mid, traction, eo.he);
        fx += traction.x() * eo.he;
        pInt += p * eo.he;
        area += eo.he;
    }
    if (meanPressure) *meanPressure = (area > 0.0) ? pInt / area : pExt;
    if (drag) *drag = fx;
    return fx;
}

void applyRightSpongeToNodes(MatrixXd& U, const CGSpace& sp, const BlastRodGeom& geom,
                             double dt, double sigmaMax, const Vector4d& referenceState) {
    if (geom.spongeWidth() <= 1e-14 || dt <= 0.0) return;
    for (int i = 0; i < sp.mesh.node.rows(); ++i) {
        double s = geom.spongeCoordinate(sp.mesh.node(i, 0));
        if (s <= 0.0) continue;
        double ramp = s * s * (3.0 - 2.0 * s);
        double alpha = std::exp(-sigmaMax * ramp * dt);
        Vector4d Ui = U.row(i).transpose();
        U.row(i) = (referenceState + alpha * (Ui - referenceState)).transpose();
    }
}

void writeCGFrame(const CGSpace& sp, const MatrixXd& U, ElasticSolid2D& solid,
                  const BlastRodGeom& geom, int renderW, int renderH, int ssaa,
                  const std::string& noMeshPath, const std::string& meshPath) {
    Mesh renderMesh = sp.mesh;
    FEM fem(1, renderMesh, false);
    MatrixXi elem2dof = sp.mesh.elem;
    int hiW = renderW * ssaa;
    int hiH = renderH * ssaa;
    std::vector<unsigned char> hiNoMesh =
        euler::renderScalarPPMImage(fem, sp.mesh, elem2dof, U.col(0), hiW, hiH,
                                    geom.xa, geom.xb,
                                    geom.ya - 0.025, geom.yb + 0.025,
                                    0.55, 3.25, euler::CM_INFERNO);
    int outW = 0, outH = 0;
    std::vector<unsigned char> noMesh =
        downsampleImage(hiNoMesh, hiW, hiH, ssaa, outW, outH);
    if (!noMesh.empty()) writePPM(noMeshPath, outW, outH, noMesh);

    std::vector<unsigned char> hiMesh = hiNoMesh;
    overlayMesh(hiMesh, hiW, hiH, sp.mesh, geom.xa, geom.xb,
                geom.ya - 0.025, geom.yb + 0.025);
    overlaySolidMesh(hiMesh, hiW, hiH, solid, geom.xa, geom.xb,
                     geom.ya - 0.025, geom.yb + 0.025, true);
    outW = 0;
    outH = 0;
    std::vector<unsigned char> meshImage =
        downsampleImage(hiMesh, hiW, hiH, ssaa, outW, outH);
    if (!meshImage.empty()) writePPM(meshPath, outW, outH, meshImage);
}

} // namespace

int runBlastRodCG(bool quick, bool freshStart) {
    (void)freshStart;
    namespace fs = std::filesystem;
    std::cout << std::unitbuf;

    const int nFrames = quick ? 36 : 900;
    const double tEnd = quick ? 0.24 : 2.16;
    const double cfl = quick ? 0.080 : 0.055;
    const double h = quick ? 0.024 : 0.012;
    const int meshIter = quick ? 120 : 260;
    const int videoW = quick ? 900 : 1500;
    const int renderH = quick ? 600 : 1000;
    const int ssaa = quick ? 2 : 3;
    const double graphViscosity = quick ? 0.34 : 0.38;
    const std::string outputPrefix = quick ? "blast_rod_cg_quick" : "blast_rod_cg";

    BlastRodGeom geom;
    const int renderW = fullDomainRenderWidth(videoW, geom);
    const int videoCropW = videoW;
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

    SolidALEMap map;
    map.setSolid(&solid);
    map.setDomain(geom.xa, geom.xb, geom.ya, geom.yb);
    map.setInfluence(0.070, 0.075);
    MatrixXd aleReferenceNodes = solid.referenceNodes();
    map.setReferenceNodes(aleReferenceNodes);
    map.setCurrent(0.0, solid.currentNodes(), solid.velocities());
    const double pExt = 1.0;

    std::cout << "Body-fitted ALE continuous-P1 CG high-pressure gas over elastic beam\n";
    std::cout << "  physical channel=" << geom.flowXb() << "x" << geom.yb
              << " computational channel=" << geom.xb << "x" << geom.yb
              << " sponge=[" << geom.spongeStartX() << "," << geom.xb << "]"
              << " rod x=" << geom.rodX << " length=" << geom.rodL
              << " width=" << geom.rodW
              << " root_radius=" << geom.rodRootRadius() << "\n";
    std::cout << "  fluid space=P1 continuous Lagrange, lumped mass, graph viscosity="
              << graphViscosity << "\n";
    std::cout << "  generating common SDF mesh h_far=" << h << "...\n";
    Mesh referenceMesh = makeBlastRodMesh(geom, h, meshIter, true);

    RefMapFn refMap = [&](const Vector2d& X, double time) { return map.refToPhys(X, time); };
    MeshVelocityFn meshVel = [&](double x, double y, double t) { return map.velocityAt(x, y, t); };
    MaxMeshSpeedFn maxSpeed = [&](double time) { return map.maxMeshSpeed(time); };

    CGSpace sp = buildCGSpace(referenceMesh, refMap, 0.0, map, geom);
    MatrixXd U(sp.mesh.node.rows(), 4);
    for (int i = 0; i < sp.mesh.node.rows(); ++i) {
        Vector4d pr = initialPrimitive(geom, sp.mesh.node(i, 0));
        U.row(i) = euler::primToCons(pr(0), pr(1), pr(2), pr(3)).transpose();
    }
    const double rhoFloor = quick ? 0.05 : 0.08;
    const double pFloor = quick ? 0.04 : 0.06;
    const double speedMax = 4.0;
    applyNodeBounds(U, rhoFloor, pFloor, speedMax);
    const Vector4d spongeReference = euler::primToCons(1.0, 0.0, 0.0, pExt);
    const double spongeSigmaMax =
        3.0 * (0.90 + std::sqrt(euler::GAMMA * pExt)) /
        std::max(geom.spongeWidth(), 1e-12);

    std::string dir = "out/" + outputPrefix + "_frames";
    std::string dirNoMesh = "out/" + outputPrefix + "_nomesh_frames";
    std::string csvPath = "out/" + outputPrefix + "_diagnostics.csv";
    fs::create_directories("out");
    fs::create_directories(dir);
    fs::create_directories(dirNoMesh);
    clearFrameDirectory(dir);
    clearFrameDirectory(dirNoMesh);

    std::ofstream diag(csvPath, std::ios::trunc);
    diag << "time,tip_displacement,tip_velocity,fluid_force_x,drag,mean_rod_pressure,"
            "fluid_nodes,fluid_triangles,solid_nodes,solid_triangles,min_h,rho_min,rho_max\n";

    std::cout << "  base nodes=" << referenceMesh.node.rows()
              << " elems=" << referenceMesh.elem.rows()
              << " minH=" << sp.minH << "\n";
    SolidMeshQuality sq = solid.meshQuality();
    std::cout << "  solid beam FEM: nodes=" << solid.numNodes()
              << " tris=" << solid.numElements()
              << " min_angle=" << sq.minAngleDeg
              << " mean_angle=" << sq.meanAngleDeg
              << " min_edge=" << sq.minEdge
              << " area=[" << sq.minArea << "," << sq.maxArea << "]"
              << " mass=" << solid.totalMass()
              << " E=" << material.young
              << " nu=" << material.poisson
              << " damping=" << material.damping << "\n";
    std::cout << "  render SSAA=" << ssaa << "x"
              << " final=" << renderW << "x" << renderH
              << " internal=" << renderW * ssaa << "x" << renderH * ssaa
              << " ppm_full_domain video_crop=" << videoCropW << "x" << renderH
              << " sponge_sigma_max=" << spongeSigmaMax << "\n";

    double t = 0.0;
    double nextFrame = 0.0;
    double frameDt = tEnd / std::max(1, nFrames);
    int step = 0;
    int frame = 0;

    auto writeFrame = [&](int idx, double time) {
        CGSpace renderSpace = buildCGSpace(referenceMesh, refMap, time, map, geom);
        char fn[512];
        char fnNoMesh[512];
        std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", dir.c_str(), idx);
        std::snprintf(fnNoMesh, sizeof(fnNoMesh), "%s/frame_%05d.ppm",
                      dirNoMesh.c_str(), idx);
        writeCGFrame(renderSpace, U, solid, geom, renderW, renderH, ssaa, fnNoMesh, fn);
    };

    writeFrame(frame++, t);
    nextFrame += frameDt;

    while (t < tEnd - 1e-14) {
        map.setCurrent(t, solid.currentNodes(), solid.velocities());
        sp = buildCGSpace(referenceMesh, refMap, t, map, geom);

        solid.clearExternalForces();
        double pMean = pExt;
        double drag = 0.0;
        double fFluid = loadSolidFromCG(sp, U, solid, pExt, &pMean, &drag);

        double dt = std::min(estimateCGDt(sp, U, maxSpeed, t, cfl), tEnd - t);
        dt = std::min(dt, solid.stableTimeStep(0.38));
        dt = std::min(dt, quick ? 0.0012 : 0.0010);

        MatrixXd solidX0 = solid.currentNodes();
        MatrixXd solidV0 = solid.velocities();
        solid.advanceExplicit(dt);
        MatrixXd solidX1 = solid.currentNodes();
        MatrixXd solidV1 = solid.velocities();
        map.setMotion(t, t + dt, solidX0, solidV0, solidX1, solidV1);

        U = advanceCG(sp, U, meshVel, t, dt, graphViscosity, rhoFloor, pFloor, speedMax);
        t += dt;
        ++step;

        map.setCurrent(t, solid.currentNodes(), solid.velocities());
        sp = buildCGSpace(referenceMesh, refMap, t, map, geom);
        applyRightSpongeToNodes(U, sp, geom, dt, spongeSigmaMax, spongeReference);
        applyNodeBounds(U, rhoFloor, pFloor, speedMax);

        double rmin = U.col(0).minCoeff();
        double rmax = U.col(0).maxCoeff();
        diag << std::setprecision(12) << t << "," << solid.tipDisplacementX()
             << "," << solid.tipVelocityX() << ","
             << fFluid << "," << drag << "," << pMean << ","
             << U.rows() << "," << sp.mesh.elem.rows() << ","
             << solid.numNodes() << "," << solid.numElements() << ","
             << sp.minH << "," << rmin << "," << rmax << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step
                      << " tipX=" << solid.tipDisplacementX()
                      << " tipV=" << solid.tipVelocityX()
                      << " Fq=" << fFluid
                      << " drag=" << drag
                      << " nodes=" << U.rows()
                      << " tris=" << sp.mesh.elem.rows()
                      << " minH=" << sp.minH
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " cg_p1 static_mesh\n";
            nextFrame += frameDt;
        }
    }

    std::string still = "out/" + outputPrefix + ".ppm";
    std::string stillNoMesh = "out/" + outputPrefix + "_nomesh.ppm";
    std::string stillPng = "out/" + outputPrefix + ".png";
    std::string stillNoMeshPng = "out/" + outputPrefix + "_nomesh.png";
    std::string video = "out/" + outputPrefix + ".mp4";
    std::string videoNoMesh = "out/" + outputPrefix + "_nomesh.mp4";
    writeCGFrame(sp, U, solid, geom, renderW, renderH, ssaa, stillNoMesh, still);
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
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
        if (std::string(argv[i]) == "--fresh") freshStart = true;
    }
    return euler_ale::runBlastRodCG(quick, freshStart);
}
