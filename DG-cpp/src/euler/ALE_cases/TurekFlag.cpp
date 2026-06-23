#include "TurekFlag.h"

#include "MeshGen.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace euler_ale {
namespace {

void addFixedSegment(DistanceMeshSpec& spec, const Vector2d& a, const Vector2d& b, double h) {
    int n = std::max(1, static_cast<int>(std::ceil((b - a).norm() / std::max(h, 1e-12))));
    for (int i = 0; i <= n; ++i) {
        double s = static_cast<double>(i) / n;
        spec.fixedPoints.push_back((1.0 - s) * a + s * b);
    }
}

} // namespace

double TurekFlagGeom::flagX1() const { return flagX0 + flagL; }
double TurekFlagGeom::flagY0() const { return cy - 0.5 * flagT; }
double TurekFlagGeom::flagY1() const { return cy + 0.5 * flagT; }

double signedBox(double x, double y, double x0, double x1, double y0, double y1) {
    double dx = std::max({x0 - x, 0.0, x - x1});
    double dy = std::max({y0 - y, 0.0, y - y1});
    double outside = std::hypot(dx, dy);
    if (dx > 0.0 || dy > 0.0) return outside;
    return -std::min({x - x0, x1 - x, y - y0, y1 - y});
}

double turekFluidSdf(const TurekFlagGeom& g, double x, double y) {
    double dRect = -std::min({x - g.xa, g.xb - x, y - g.ya, g.yb - y});
    double dCyl = std::hypot(x - g.cx, y - g.cy) - g.r;
    double dFlag = signedBox(x, y, g.flagX0, g.flagX1(), g.flagY0(), g.flagY1());
    double dObstacle = std::min(dCyl, dFlag);
    return std::max(dRect, -dObstacle);
}

double TurekFlagMap::smooth01(double s) {
    s = std::clamp(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
}

double TurekFlagMap::q(double t) const {
    double a = std::clamp((t - t0) / std::max(1e-14, t1 - t0), 0.0, 1.0);
    return (1.0 - a) * q0 + a * q1;
}

double TurekFlagMap::qd(double) const {
    return (q1 - q0) / std::max(1e-14, t1 - t0);
}

double TurekFlagMap::mode(double x) const {
    if (x <= geom.flagX0) return 0.0;
    if (x <= geom.flagX1()) {
        double s = (x - geom.flagX0) / geom.flagL;
        return smooth01(s);
    }
    double z = (x - geom.flagX1()) / wakeDecay;
    return std::exp(-z * z);
}

double TurekFlagMap::influence(double x, double y) const {
    double dBand = std::max(0.0, std::abs(y - geom.cy) - 0.5 * geom.flagT);
    double band = std::exp(-(dBand / normalDecay) * (dBand / normalDecay));
    double dWall = std::min(y - geom.ya, geom.yb - y);
    double wall = smooth01(dWall / wallMargin);
    return mode(x) * band * wall;
}

double TurekFlagMap::flagMode(double x) const {
    double s = std::clamp((x - geom.flagX0) / geom.flagL, 0.0, 1.0);
    return smooth01(s);
}

double TurekFlagMap::flagY(double x, double yRef, double t) const {
    return yRef + q(t) * flagMode(x);
}

Vector2d TurekFlagMap::refToPhys(const Vector2d& X, double t) const {
    return Vector2d(X.x(), X.y() + q(t) * influence(X.x(), X.y()));
}

Vector2d TurekFlagMap::velocityAt(double x, double y, double t) const {
    return Vector2d(0.0, qd(t) * influence(x, y));
}

void TurekFlagMap::setLinearMotion(double a, double b, double qa, double qb) {
    t0 = a;
    t1 = b;
    q0 = qa;
    q1 = qb;
}

double TurekFlagMap::maxMeshSpeed(double t) const { return std::abs(qd(t)); }

ModalState advanceFlagSymplectic(ModalState state, double fluidForce,
                                 const FlagModalParams& params, double dt) {
    double a = (fluidForce - params.stiffness * state.q - params.damping * state.qd) / params.mass;
    state.qd += dt * a;
    state.q += dt * state.qd;
    if (state.q < params.qMin) {
        state.q = params.qMin;
        state.qd = std::max(0.0, state.qd);
    }
    if (state.q > params.qMax) {
        state.q = params.qMax;
        state.qd = std::min(0.0, state.qd);
    }
    return state;
}

Mesh makeTurekFlagMesh(const TurekFlagGeom& g, double h, int maxIter, bool verbose) {
    DistanceMeshSpec spec;
    spec.xa = g.xa;
    spec.xb = g.xb;
    spec.ya = g.ya;
    spec.yb = g.yb;
    spec.h0 = h;
    spec.signedDistance = [g](double x, double y) { return turekFluidSdf(g, x, y); };
    spec.targetSize = [g, h](double x, double y) {
        double dCyl = std::max(0.0, std::hypot(x - g.cx, y - g.cy) - g.r);
        double dFlag = std::max(0.0, signedBox(x, y, g.flagX0, g.flagX1(), g.flagY0(), g.flagY1()));
        double dNear = std::min(dCyl, dFlag);
        double wake = std::max(0.0, x - g.flagX0);
        double wakeBand = std::abs(y - g.cy);
        double wakeDist = (x <= g.flagX1() + 0.55) ? std::max(0.0, wakeBand - 0.10) : 1.0;
        double local = h * (1.0 + 3.5 * dNear / h);
        local = std::min(local, h * (1.0 + 2.0 * wakeDist / h));
        if (wake > 0.0 && x < g.flagX1() + 0.45) local = std::min(local, 1.6 * h);
        return std::min(3.2 * h, std::max(0.45 * h, local));
    };

    double hWall = 1.9 * h;
    addFixedSegment(spec, Vector2d(g.xa, g.ya), Vector2d(g.xb, g.ya), hWall);
    addFixedSegment(spec, Vector2d(g.xb, g.ya), Vector2d(g.xb, g.yb), hWall);
    addFixedSegment(spec, Vector2d(g.xb, g.yb), Vector2d(g.xa, g.yb), hWall);
    addFixedSegment(spec, Vector2d(g.xa, g.yb), Vector2d(g.xa, g.ya), hWall);

    constexpr double pi = 3.14159265358979323846;
    int nCyl = std::max(24, static_cast<int>(std::round(2.0 * pi * g.r / (0.65 * h))));
    for (int k = 0; k < nCyl; ++k) {
        double th = 2.0 * pi * k / nCyl;
        spec.fixedPoints.emplace_back(g.cx + g.r * std::cos(th), g.cy + g.r * std::sin(th));
    }

    double hf = std::min(0.85 * h, 0.026);
    addFixedSegment(spec, Vector2d(g.flagX0, g.flagY0()), Vector2d(g.flagX1(), g.flagY0()), hf);
    addFixedSegment(spec, Vector2d(g.flagX1(), g.flagY0()), Vector2d(g.flagX1(), g.flagY1()), hf);
    addFixedSegment(spec, Vector2d(g.flagX1(), g.flagY1()), Vector2d(g.flagX0, g.flagY1()), hf);
    addFixedSegment(spec, Vector2d(g.flagX0, g.flagY1()), Vector2d(g.flagX0, g.flagY0()), hf);

    Mesh mesh;
    generateDistanceMesh(mesh, spec, maxIter, verbose);
    return mesh;
}

int turekBoundaryTag(double x, double y, double t, const TurekFlagMap& map) {
    const TurekFlagGeom& g = map.geom;
    double dLeft = std::abs(x - g.xa);
    double dRight = std::abs(x - g.xb);
    double dBottom = std::abs(y - g.ya);
    double dTop = std::abs(y - g.yb);
    double dCyl = std::abs(std::hypot(x - g.cx, y - g.cy) - g.r);

    double yLo = map.flagY(std::clamp(x, g.flagX0, g.flagX1()), g.flagY0(), t);
    double yHi = map.flagY(std::clamp(x, g.flagX0, g.flagX1()), g.flagY1(), t);
    double dFlag = 1e9;
    if (x >= g.flagX0 - 0.025 && x <= g.flagX1() + 0.025) {
        dFlag = std::min(dFlag, std::abs(y - yLo));
        dFlag = std::min(dFlag, std::abs(y - yHi));
    }
    if (y >= yLo - 0.025 && y <= yHi + 0.025) {
        dFlag = std::min(dFlag, std::abs(x - g.flagX0));
        dFlag = std::min(dFlag, std::abs(x - g.flagX1()));
    }

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
    if (dCyl < best) {
        best = dCyl;
        tag = TAG_SLIP_WALL;
    }
    if (dFlag < best) tag = TAG_MOVING_WALL;
    return tag;
}

double flagGeneralizedForce(const Space& sp, const MatrixXd& U, const TurekFlagMap& map,
                            double time, double pExt, double* meanPressure,
                            double* lift, double* drag) {
    (void)time;
    MatrixXd q1d;
    VectorXd w1d;
    sp.fem->quad1d(q1d, w1d);
    int nqe = static_cast<int>(w1d.size());
    int locDof = sp.fem->locDof;
    double force = 0.0;
    double area = 0.0;
    double pInt = 0.0;
    double fy = 0.0;
    double fx = 0.0;
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) != TAG_MOVING_WALL) continue;
        int tt = (sp.e2s(e, 0) >= 0) ? sp.e2s(e, 0) : sp.e2s(e, 1);
        if (tt < 0) continue;
        int n1 = sp.edge(e, 0);
        int n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, tt, n1, n2);
        Matrix<double, Dynamic, 4> Ue(locDof, 4);
        for (int i = 0; i < locDof; ++i) Ue.row(i) = U.row(sp.e2d(tt, i));
        for (int q = 0; q < nqe; ++q) {
            double l1 = q1d(q, 0);
            double l2 = q1d(q, 1);
            Vector3d lam = Vector3d::Zero();
            if (eo.et == 0) {
                lam(0) = (eo.dir == 0) ? l1 : l2;
                lam(1) = (eo.dir == 0) ? l2 : l1;
            } else if (eo.et == 1) {
                lam(1) = (eo.dir == 0) ? l1 : l2;
                lam(2) = (eo.dir == 0) ? l2 : l1;
            } else {
                lam(2) = (eo.dir == 0) ? l1 : l2;
                lam(0) = (eo.dir == 0) ? l2 : l1;
            }
            RowVectorXd phi = sp.fem->computeBasisValue_all(lam.transpose()).row(0);
            Vector4d Uq = (phi * Ue).transpose();
            double p = euler::pressure(Uq);
            Vector2d Pp = l1 * sp.mesh.node.row(n1).transpose() +
                          l2 * sp.mesh.node.row(n2).transpose();
            double ds = w1d(q) * eo.he;
            double dp = p - pExt;
            double mode = map.flagMode(Pp.x());
            force += dp * eo.nout.y() * mode * ds;
            fx += dp * eo.nout.x() * ds;
            fy += dp * eo.nout.y() * ds;
            pInt += p * ds;
            area += ds;
        }
    }
    if (meanPressure) *meanPressure = (area > 0.0) ? pInt / area : pExt;
    if (lift) *lift = fy;
    if (drag) *drag = fx;
    return force;
}

int runTurek(bool quick) {
    namespace fs = std::filesystem;
    int ord = 1;
    int maxGen = quick ? 1 : 2;
    int remeshEvery = quick ? 10 : 8;
    int nFrames = quick ? 30 : 80;
    double tEnd = quick ? 0.18 : 0.32;
    double cfl = 0.22;
    double h = quick ? 0.065 : 0.055;
    int meshIter = quick ? 70 : 100;
    double uMax = 0.30;
    double thRef = 0.010;
    double thCrs = 0.003;

    TurekFlagGeom geom;
    TurekFlagMap map;
    map.geom = geom;
    FlagModalParams params;
    ModalState flag;
    flag.q = 0.004;
    flag.qd = 0.0;

    std::cout << "Turek-Hron-style cylinder + flexible flag ALE-AMR benchmark surrogate\n";
    std::cout << "  channel=" << geom.xb << "x" << geom.yb
              << " cylinder r=" << geom.r << " center=(" << geom.cx << "," << geom.cy << ")"
              << " flag=" << geom.flagL << "x" << geom.flagT << "\n";
    std::cout << "  generating common SDF mesh h=" << h << "...\n";
    Mesh base = makeTurekFlagMesh(geom, h, meshIter, true);
    double minAng = 0.0;
    double meanAng = 0.0;
    double minArea = 0.0;
    double maxArea = 0.0;
    meshQuality(base, minAng, meanAng, minArea, maxArea);

    RefMapFn refMap = [&](const Vector2d& X, double time) { return map.refToPhys(X, time); };
    MaxMeshSpeedFn maxSpeed = [&](double time) { return map.maxMeshSpeed(time); };
    auto inflowU = [=](double y) {
        double yy = std::clamp((y - geom.ya) / (geom.yb - geom.ya), 0.0, 1.0);
        return uMax * 4.0 * yy * (1.0 - yy);
    };
    Tagger tagger = [&](double x, double y, double t) { return turekBoundaryTag(x, y, t, map); };
    MeshVelocityFn meshVel = [&](double x, double y, double t) { return map.velocityAt(x, y, t); };
    ALEBCFn bc = [&](double, double y, double, const Vector4d& Um,
                     double nxn, double nyn, int tag, double wn) {
        if (tag == TAG_MOVING_WALL || tag == TAG_SLIP_WALL)
            return movingWallGhost(Um, nxn, nyn, (tag == TAG_MOVING_WALL) ? wn : 0.0);
        if (tag == TAG_EXACT)
            return euler::primToCons(1.0, inflowU(y), 0.0, 1.0);
        return Um;
    };

    map.setLinearMotion(0.0, 1e-12, flag.q, flag.q);
    ALEAdaptiveForest forest(base, ord, 4);
    Space sp;
    rebuildSpace(forest, ord, refMap, 0.0, tagger, sp);
    MatrixXd U = euler::projectInitial(*sp.fem, sp.mesh, sp.e2d,
                                       [&](double, double y) {
                                           return Vector4d(1.0, inflowU(y), 0.0, 1.0);
                                       });

    std::string dir = quick ? "out/turek_flag_quick_frames" : "out/turek_flag_frames";
    std::string csvPath = quick ? "out/turek_flag_quick_diagnostics.csv" :
                                  "out/turek_flag_diagnostics.csv";
    fs::create_directories("out");
    fs::create_directories(dir);
    for (const auto& e : fs::directory_iterator(dir))
        if (e.path().extension() == ".ppm") fs::remove(e.path());
    std::ofstream diag(csvPath);
    diag << "time,tip_displacement,tip_velocity,generalized_force,lift,drag,mean_flag_pressure,"
            "triangles,rho_min,rho_max\n";

    auto setCurrentMap = [&](double time) {
        map.setLinearMotion(time, time + 1e-12, flag.q, flag.q);
    };
    auto writeFrame = [&](int idx, double time) {
        rebuildSpace(forest, ord, refMap, time, tagger, sp);
        char fn[512];
        std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm", dir.c_str(), idx);
        int W = quick ? 1200 : 1500;
        int H = quick ? 260 : 320;
        euler::writeScalarPPM(fn, *sp.fem, sp.mesh, sp.e2d, U.col(0), W, H,
                              geom.xa - 0.02, geom.xb + 0.02,
                              geom.ya - 0.025, geom.yb + 0.025,
                              0.965, 1.035, euler::CM_COOLWARM);
        overlayMesh(fn, sp.mesh, geom.xa - 0.02, geom.xb + 0.02,
                    geom.ya - 0.025, geom.yb + 0.025);
    };

    std::cout << "  base nodes=" << base.node.rows() << " elems=" << base.elem.rows()
              << " min_angle=" << std::setprecision(3) << minAng
              << " mean_angle=" << meanAng << " max_gen=" << maxGen << "\n";
    std::cout << "  modal flag: m=" << params.mass << " k=" << params.stiffness
              << " c=" << params.damping << " q0=" << flag.q << " uMax=" << uMax << "\n";

    double t = 0.0;
    double nextFrame = 0.0;
    double frameDt = tEnd / std::max(1, nFrames);
    int step = 0;
    int frame = 0;
    int lastRef = 0;
    int lastCrs = 0;
    setCurrentMap(t);
    writeFrame(frame++, t);
    nextFrame += frameDt;

    while (t < tEnd - 1e-14) {
        setCurrentMap(t);
        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        double pMean = params.pExt;
        double lift = 0.0;
        double drag = 0.0;
        double fFluid = flagGeneralizedForce(sp, U, map, t, params.pExt, &pMean, &lift, &drag);
        double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
        dt = std::min(dt, quick ? 0.0015 : 0.0012);

        ModalState nextFlag = advanceFlagSymplectic(flag, fFluid, params, dt);
        map.setLinearMotion(t, t + dt, flag.q, nextFlag.q);
        U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
        flag = nextFlag;
        t += dt;
        ++step;

        rebuildSpace(forest, ord, refMap, t, tagger, sp);
        if (step % remeshEvery == 0 && t < tEnd - 1e-14) {
            std::vector<int> flagAmr = computeAMRFlags(sp, U, forest, maxGen, thRef, thCrs, 2, true);
            forest.syncFromState(U, sp.e2d, sp.fem->locDof);
            auto rc = forest.adapt(flagAmr, maxGen);
            lastRef = rc.first;
            lastCrs = rc.second;
            rebuildSpace(forest, ord, refMap, t, tagger, sp);
            U = forest.gatherState(sp.e2d, sp.fem->locDof, sp.nDof);
        }

        double rmin = U.col(0).minCoeff();
        double rmax = U.col(0).maxCoeff();
        diag << std::setprecision(12) << t << "," << flag.q << "," << flag.qd << ","
             << fFluid << "," << lift << "," << drag << "," << pMean << ","
             << sp.mesh.elem.rows() << "," << rmin << "," << rmax << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++, t);
            int gmax = 0;
            for (int k = 0; k < sp.mesh.elem.rows(); ++k) gmax = std::max(gmax, forest.gen(k));
            std::cout << "  t=" << std::fixed << std::setprecision(4) << t
                      << " step=" << step << " qTip=" << flag.q << " qd=" << flag.qd
                      << " Fq=" << fFluid << " lift=" << lift << " drag=" << drag
                      << " tris=" << sp.mesh.elem.rows()
                      << " rho[" << std::setprecision(3) << rmin << "," << rmax << "]"
                      << " gmax=" << gmax << " remesh(+" << lastRef << ",-" << lastCrs << ")\n";
            nextFrame += frameDt;
        }
    }

    std::string still = quick ? "out/turek_flag_quick.ppm" : "out/turek_flag.ppm";
    rebuildSpace(forest, ord, refMap, t, tagger, sp);
    euler::writeScalarPPM(still, *sp.fem, sp.mesh, sp.e2d, U.col(0), quick ? 1200 : 1500,
                          quick ? 260 : 320,
                          geom.xa - 0.02, geom.xb + 0.02,
                          geom.ya - 0.025, geom.yb + 0.025,
                          0.965, 1.035, euler::CM_COOLWARM);
    overlayMesh(still, sp.mesh, geom.xa - 0.02, geom.xb + 0.02,
                geom.ya - 0.025, geom.yb + 0.025);
    std::cout << "Done. frames=" << frame << ", still=" << still
              << ", diagnostics=" << csvPath << "\n";
    std::cout << "ffmpeg -y -framerate 25 -i " << dir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 "
              << (quick ? "out/turek_flag_quick.mp4" : "out/turek_flag.mp4") << "\n";
    return 0;
}

} // namespace euler_ale

int main(int argc, char** argv) {
    bool quick = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
    }
    return euler_ale::runTurek(quick);
}
