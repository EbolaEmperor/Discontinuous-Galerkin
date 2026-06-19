#include "Convergence.h"

#include "Core.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace euler_ale {
namespace {

struct ConvergenceMap {
    double amp = 0.05;
    double period = 1.0;

    double s(double t) const {
        constexpr double pi = 3.14159265358979323846;
        return amp * std::sin(2.0 * pi * t / period);
    }

    double sd(double t) const {
        constexpr double pi = 3.14159265358979323846;
        return amp * (2.0 * pi / period) * std::cos(2.0 * pi * t / period);
    }

    Vector2d refToPhys(const Vector2d& X, double t) const {
        double st = s(t);
        return Vector2d(st + (1.0 - st) * X.x(), X.y());
    }

    Vector2d velocityAt(double x, double, double t) const {
        double st = s(t);
        double X = (x - st) / std::max(1e-12, 1.0 - st);
        return Vector2d(sd(t) * (1.0 - X), 0.0);
    }

    double maxMeshSpeed(double t) const { return std::abs(sd(t)); }
};

Vector4d exactVortexPrim(double x, double y, double t) {
    constexpr double pi = 3.14159265358979323846;
    const double gamma = euler::GAMMA;
    const double beta = 2.0;
    const double u0 = 0.45;
    const double v0 = 0.10;
    const double xc = 0.38 + u0 * t;
    const double yc = 0.50 + v0 * t;
    double X = x - xc;
    double Y = y - yc;
    double r2 = X * X + Y * Y;
    double e = std::exp(0.5 * (1.0 - r2));
    double du = -(beta / (2.0 * pi)) * e * Y;
    double dv = (beta / (2.0 * pi)) * e * X;
    double temp = 1.0 - (gamma - 1.0) * beta * beta / (8.0 * gamma * pi * pi) *
                         std::exp(1.0 - r2);
    double rho = std::pow(temp, 1.0 / (gamma - 1.0));
    double p = std::pow(temp, gamma / (gamma - 1.0));
    return Vector4d(rho, u0 + du, v0 + dv, p);
}

Vector4d exactVortexCons(double x, double y, double t) {
    Vector4d pr = exactVortexPrim(x, y, t);
    return euler::primToCons(pr(0), pr(1), pr(2), pr(3));
}

double densityL2(FEM& fem, const Mesh& mesh, const MatrixXi& e2d,
                 const MatrixXd& U, double t) {
    return euler::l2Error(fem, mesh, e2d, U, 0,
                          [=](double x, double y) { return exactVortexCons(x, y, t); });
}

} // namespace

int runConvergence() {
    int ord = 1;
    double tEnd = 0.05;
    double cfl = 0.08;
    ConvergenceMap map;

    RefMapFn refMap = [&](const Vector2d& X, double time) { return map.refToPhys(X, time); };
    MaxMeshSpeedFn maxSpeed = [&](double time) { return map.maxMeshSpeed(time); };
    Tagger tagger = [](double, double, double) { return TAG_EXACT; };
    MeshVelocityFn meshVel = [&](double x, double y, double t) { return map.velocityAt(x, y, t); };
    ALEBCFn bc = [](double x, double y, double t, const Vector4d&, double, double, int, double) {
        return exactVortexCons(x, y, t);
    };

    std::vector<int> Ns = {8, 12, 16, 24, 32};
    std::vector<double> hs, errs;
    std::cout << "ALE-DG moving-mesh convergence: isentropic vortex, dP" << ord << "\n";
    std::cout << "   N        h          L2(rho)       order       steps\n";
    double prevErr = 0.0;
    double prevH = 0.0;
    for (int N : Ns) {
        Mesh base;
        euler::makeRectMesh(base, 0.0, 1.0, 0.0, 1.0, N, N);
        ALEAdaptiveForest forest(base, ord, 4);
        Space sp;
        rebuildSpace(forest, ord, refMap, 0.0, tagger, sp);
        MatrixXd U = euler::projectInitial(*sp.fem, sp.mesh, sp.e2d,
                                           [](double x, double y) {
                                               return exactVortexPrim(x, y, 0.0);
                                           });
        double t = 0.0;
        int steps = 0;
        while (t < tEnd - 1e-14) {
            rebuildSpace(forest, ord, refMap, t, tagger, sp);
            double dt = std::min(estimateDt(sp, U, maxSpeed, t, ord, cfl), tEnd - t);
            U = advanceOne(forest, ord, refMap, t, dt, tagger, meshVel, bc, U);
            t += dt;
            ++steps;
        }
        rebuildSpace(forest, ord, refMap, tEnd, tagger, sp);
        double err = densityL2(*sp.fem, sp.mesh, sp.e2d, U, tEnd);
        double h = 1.0 / N;
        double order = (prevErr > 0.0) ? std::log(prevErr / err) / std::log(prevH / h) : 0.0;
        std::cout << std::setw(4) << N << "  "
                  << std::scientific << std::setprecision(3) << h << "  "
                  << err << "  ";
        if (prevErr > 0.0) std::cout << std::fixed << std::setprecision(2) << order;
        else std::cout << "  --";
        std::cout << "      " << steps << "\n";
        prevErr = err;
        prevH = h;
        hs.push_back(h);
        errs.push_back(err);
    }
    return 0;
}

} // namespace euler_ale

int main() {
    return euler_ale::runConvergence();
}
