#include "EulerALE1D.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace euler_fsi;

namespace {

constexpr double L = 1.0;
constexpr double U0 = 0.45;
constexpr double P0 = 1.0;
constexpr double AMP = 0.18;
constexpr double MESH_AMP = 0.055;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double OMEGA = 2.0 * PI;

double wrap(double x) {
    x = std::fmod(x, L);
    if (x < 0.0) x += L;
    return x;
}

Prim1D exactPrim(double x, double t) {
    double xp = wrap(x - U0 * t);
    double rho = 1.0 + AMP * std::sin(2.0 * PI * xp / L);
    return {rho, U0, P0};
}

double exactRho(double x, double t) {
    return exactPrim(x, t).rho;
}

double minDx(const Grid1D& g) {
    double h = 1e300;
    for (int i = 0; i < g.cells(); ++i) h = std::min(h, g.dx(i));
    return h;
}

double maxFaceSpeed(const Grid1D& g) {
    double m = 0.0;
    for (double w : g.wFace) m = std::max(m, std::abs(w));
    return m;
}

double runEntropyWave(int n, double dtMax, double T, bool movingMesh) {
    auto gridFn = [=](double t) {
        return movingMesh ? Grid1D::oscillatory(n, L, t, MESH_AMP, OMEGA)
                          : Grid1D::affinePiston(n, L, 0.0);
    };
    auto bcFn = [](double) {
        BoundaryCondition1D bc;
        bc.left = BoundaryKind::Periodic;
        bc.right = BoundaryKind::Periodic;
        return bc;
    };
    EulerALE1D solver(n);
    solver.setCellAverages([](double x) { return exactPrim(x, 0.0); }, gridFn(0.0));
    solver.setGrid(gridFn(0.0), 0.0);
    double t = 0.0;
    while (t < T - 1e-14) {
        Grid1D g = gridFn(t);
        double cflDt = 0.38 * minDx(g) / (solver.maxWaveSpeed() + maxFaceSpeed(g) + 1e-12);
        double dt = std::min({dtMax, cflDt, T - t});
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Unlimited)) {
            std::cerr << "entropy-wave step failed at t=" << t << "\n";
            return 1e300;
        }
        t += dt;
    }
    return densityL2Error(solver, exactRho, T);
}

std::vector<double> runEntropyDensity(int n, double dtMax, double T) {
    auto gridFn = [=](double) { return Grid1D::affinePiston(n, L, 0.0); };
    auto bcFn = [](double) {
        BoundaryCondition1D bc;
        bc.left = BoundaryKind::Periodic;
        bc.right = BoundaryKind::Periodic;
        return bc;
    };
    EulerALE1D solver(n);
    solver.setCellAverages([](double x) { return exactPrim(x, 0.0); }, gridFn(0.0));
    solver.setGrid(gridFn(0.0), 0.0);
    double t = 0.0;
    while (t < T - 1e-14) {
        double dt = std::min(dtMax, T - t);
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Unlimited)) {
            std::cerr << "entropy-density step failed at t=" << t << "\n";
            break;
        }
        t += dt;
    }
    std::vector<double> rho(n);
    for (int i = 0; i < n; ++i) rho[i] = solver.cellAverage(i).rho;
    return rho;
}

double densityDifference(const std::vector<double>& a, const std::vector<double>& b) {
    double dx = L / static_cast<double>(a.size());
    double e = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        e += dx * d * d;
    }
    return std::sqrt(e);
}

double gclUniformError(int n) {
    const double T = 1.0;
    auto gridFn = [=](double t) { return Grid1D::oscillatory(n, L, t, 0.08, OMEGA); };
    auto bcFn = [](double) {
        BoundaryCondition1D bc;
        bc.left = BoundaryKind::Periodic;
        bc.right = BoundaryKind::Periodic;
        return bc;
    };
    Prim1D q0{1.0, 0.2, 1.0};
    EulerALE1D solver(n);
    solver.setUniform(q0, gridFn(0.0));
    solver.setGrid(gridFn(0.0), 0.0);
    double t = 0.0;
    while (t < T - 1e-14) {
        Grid1D g = gridFn(t);
        double dt = std::min(0.30 * minDx(g) / (solver.maxWaveSpeed() + maxFaceSpeed(g) + 1e-12),
                             T - t);
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Unlimited)) return 1e300;
        t += dt;
    }
    double err = 0.0;
    for (int i = 0; i < solver.cells(); ++i) {
        Prim1D q = solver.cellPrimitive(i);
        err = std::max(err, std::abs(q.rho - q0.rho));
        err = std::max(err, std::abs(q.u - q0.u));
        err = std::max(err, std::abs(q.p - q0.p));
    }
    return err;
}

} // namespace

int main() {
    std::cout << "============================================================\n";
    std::cout << " Euler-ALE/FSI infrastructure verification\n";
    std::cout << " 1D conservative ALE finite-volume core, moving mesh + Euler\n";
    std::cout << "============================================================\n";

    std::cout << "\n[1] Geometric conservation law / uniform-flow preservation\n";
    double maxGcl = 0.0;
    for (int n : {40, 80, 160, 320}) {
        double e = gclUniformError(n);
        maxGcl = std::max(maxGcl, e);
        std::cout << "  n=" << std::setw(4) << n << "  max primitive drift="
                  << std::scientific << std::setprecision(3) << e << "\n";
    }

    std::cout << "\n[2] Spatial convergence on a moving mesh\n";
    const double T = 0.2;
    double prevE = 0.0, prevH = 0.0;
    double minSpatialRate = 1e300;
    double finalSpatialErr = 0.0;
    for (int n : {50, 100, 200, 400}) {
        double h = L / n;
        double err = runEntropyWave(n, 0.08 * h, T, true);
        double rate = prevE > 0.0 ? std::log(prevE / err) / std::log(prevH / h) : 0.0;
        if (prevE > 0.0) minSpatialRate = std::min(minSpatialRate, rate);
        finalSpatialErr = err;
        std::cout << "  n=" << std::setw(4) << n << "  h=" << std::fixed << std::setprecision(5) << h
                  << "  L2(rho)=" << std::scientific << std::setprecision(3) << err
                  << "  rate=" << std::fixed << std::setprecision(2) << rate << "\n";
        prevE = err; prevH = h;
    }
    std::cout << "  expected asymptotic rate: 2\n";

    std::cout << "\n[3] Temporal convergence on a fixed mesh (successive differences)\n";
    std::vector<double> dtList = {0.0010, 0.0005, 0.00025, 0.000125, 0.0000625};
    std::vector<std::vector<double>> sols;
    for (double dt : dtList) sols.push_back(runEntropyDensity(220, dt, T));
    prevE = 0.0;
    double minTemporalRate = 1e300;
    for (size_t i = 0; i + 1 < sols.size(); ++i) {
        double err = densityDifference(sols[i], sols[i + 1]);
        double rate = prevE > 0.0 ? std::log(prevE / err) / std::log(2.0) : 0.0;
        if (prevE > 0.0) minTemporalRate = std::min(minTemporalRate, rate);
        std::cout << "  dt=" << std::scientific << std::setprecision(3) << dtList[i]
                  << " -> " << dtList[i + 1]
                  << "  ||rho(dt)-rho(dt/2)||_L2=" << err
                  << "  rate=" << std::fixed << std::setprecision(2) << rate << "\n";
        prevE = err;
    }
    std::cout << "  expected asymptotic rate: 2\n";

    bool ok = true;
    ok = ok && std::isfinite(maxGcl) && maxGcl < 1e-10;
    ok = ok && std::isfinite(finalSpatialErr) && finalSpatialErr < 4e-6;
    ok = ok && std::isfinite(minSpatialRate) && minSpatialRate > 1.85;
    ok = ok && std::isfinite(minTemporalRate) && minTemporalRate > 1.90;

    std::cout << "\nVerification " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) {
        std::cerr << "Euler-ALE convergence gate failed: maxGcl=" << maxGcl
                  << " finalSpatialErr=" << finalSpatialErr
                  << " minSpatialRate=" << minSpatialRate
                  << " minTemporalRate=" << minTemporalRate << "\n";
        return 2;
    }
    return 0;
}
