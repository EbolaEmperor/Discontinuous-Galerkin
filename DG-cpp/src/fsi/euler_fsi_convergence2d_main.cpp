#include "EulerALE2D.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace euler_fsi;

namespace {

constexpr double LX = 1.0;
constexpr double LY = 0.75;
constexpr double U0 = 0.37;
constexpr double V0 = -0.21;
constexpr double P0 = 1.0;
constexpr double AMP = 0.10;
constexpr double MESH_AMP = 0.045;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double OMEGA = 2.0 * PI;

double wrap(double x, double L) {
    x = std::fmod(x, L);
    if (x < 0.0) x += L;
    return x;
}

Prim2D exactPrim(double x, double y, double t) {
    double xp = wrap(x - U0 * t, LX);
    double yp = wrap(y - V0 * t, LY);
    double rho = 1.0 + AMP * std::sin(2.0 * PI * xp / LX) *
                         std::cos(2.0 * PI * yp / LY);
    return {rho, U0, V0, P0};
}

double exactRho(double x, double y, double t) {
    return exactPrim(x, y, t).rho;
}

double minSpacing(const Grid2D& g) {
    double h = 1e300;
    for (int i = 0; i < g.nx; ++i) h = std::min(h, g.dx(i));
    for (int j = 0; j < g.ny; ++j) h = std::min(h, g.dy(j));
    return h;
}

double maxFaceSpeed(const Grid2D& g) {
    double m = 0.0;
    for (double w : g.wXFace) m = std::max(m, std::abs(w));
    for (double w : g.wYFace) m = std::max(m, std::abs(w));
    return m;
}

BoundaryCondition2D periodicBC() {
    BoundaryCondition2D bc;
    bc.left = BoundaryKind::Periodic;
    bc.right = BoundaryKind::Periodic;
    bc.bottom = BoundaryKind::Periodic;
    bc.top = BoundaryKind::Periodic;
    return bc;
}

double runEntropyWave2D(int nx, int ny, double dtMax, double T, bool movingMesh) {
    auto gridFn = [=](double t) {
        return movingMesh ? Grid2D::oscillatoryX(nx, ny, LX, LY, t, MESH_AMP, OMEGA)
                          : Grid2D::affinePiston(nx, ny, LX, LY, 0.0);
    };
    auto bcFn = [](double) { return periodicBC(); };
    EulerALE2D solver(nx, ny);
    solver.setCellAverages([](double x, double y) { return exactPrim(x, y, 0.0); }, gridFn(0.0));
    solver.setGrid(gridFn(0.0), 0.0);
    double t = 0.0;
    while (t < T - 1e-14) {
        Grid2D g = gridFn(t);
        double cflDt = 0.32 * minSpacing(g) / (solver.maxWaveSpeed() + maxFaceSpeed(g) + 1e-12);
        double dt = std::min({dtMax, cflDt, T - t});
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Unlimited)) return 1e300;
        t += dt;
    }
    return densityL2Error(solver, exactRho, T);
}

std::vector<double> runEntropyDensity2D(int nx, int ny, double dtMax, double T) {
    auto gridFn = [=](double) { return Grid2D::affinePiston(nx, ny, LX, LY, 0.0); };
    auto bcFn = [](double) { return periodicBC(); };
    EulerALE2D solver(nx, ny);
    solver.setCellAverages([](double x, double y) { return exactPrim(x, y, 0.0); }, gridFn(0.0));
    solver.setGrid(gridFn(0.0), 0.0);
    double t = 0.0;
    while (t < T - 1e-14) {
        double dt = std::min(dtMax, T - t);
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Unlimited)) break;
        t += dt;
    }
    std::vector<double> rho(static_cast<size_t>(nx) * ny);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) rho[static_cast<size_t>(j) * nx + i] = solver.cellAverage(i, j).rho;
    return rho;
}

double densityDifference2D(const std::vector<double>& a, const std::vector<double>& b,
                           int nx, int ny) {
    double dx = LX / nx, dy = LY / ny;
    double e = 0.0;
    for (size_t k = 0; k < a.size(); ++k) {
        double d = a[k] - b[k];
        e += dx * dy * d * d;
    }
    return std::sqrt(e);
}

double gclUniformError2D(int nx, int ny) {
    const double T = 1.0;
    auto gridFn = [=](double t) { return Grid2D::oscillatoryX(nx, ny, LX, LY, t, 0.07, OMEGA); };
    auto bcFn = [](double) { return periodicBC(); };
    Prim2D q0{1.0, 0.18, -0.11, 1.0};
    EulerALE2D solver(nx, ny);
    solver.setUniform(q0, gridFn(0.0));
    solver.setGrid(gridFn(0.0), 0.0);
    double t = 0.0;
    while (t < T - 1e-14) {
        Grid2D g = gridFn(t);
        double dt = std::min(0.30 * minSpacing(g) / (solver.maxWaveSpeed() + maxFaceSpeed(g) + 1e-12),
                             T - t);
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Unlimited)) return 1e300;
        t += dt;
    }
    double err = 0.0;
    for (int j = 0; j < solver.ny(); ++j)
        for (int i = 0; i < solver.nx(); ++i) {
            Prim2D q = solver.cellPrimitive(i, j);
            err = std::max(err, std::abs(q.rho - q0.rho));
            err = std::max(err, std::abs(q.u - q0.u));
            err = std::max(err, std::abs(q.v - q0.v));
            err = std::max(err, std::abs(q.p - q0.p));
        }
    return err;
}

} // namespace

int main() {
    std::cout << "============================================================\n";
    std::cout << " Euler-ALE/FSI 2D infrastructure verification\n";
    std::cout << " structured 2D conservative ALE FV core\n";
    std::cout << "============================================================\n";

    std::cout << "\n[1] 2D GCL / uniform-flow preservation\n";
    double maxGcl = 0.0;
    for (auto dims : {std::pair<int, int>{24, 18}, {48, 36}, {96, 72}}) {
        double e = gclUniformError2D(dims.first, dims.second);
        maxGcl = std::max(maxGcl, e);
        std::cout << "  nx=" << std::setw(4) << dims.first
                  << " ny=" << std::setw(4) << dims.second
                  << "  max primitive drift=" << std::scientific << std::setprecision(3) << e << "\n";
    }

    std::cout << "\n[2] 2D spatial convergence on a moving mesh\n";
    const double T = 0.16;
    double prevE = 0.0, prevH = 0.0, minSpatialRate = 1e300, finalSpatialErr = 0.0;
    for (int nx : {24, 48, 96, 192}) {
        int ny = 3 * nx / 4;
        double h = LX / nx;
        double err = runEntropyWave2D(nx, ny, 0.07 * h, T, true);
        double rate = prevE > 0.0 ? std::log(prevE / err) / std::log(prevH / h) : 0.0;
        if (prevE > 0.0) minSpatialRate = std::min(minSpatialRate, rate);
        finalSpatialErr = err;
        std::cout << "  nx=" << std::setw(4) << nx << " ny=" << std::setw(4) << ny
                  << "  h=" << std::fixed << std::setprecision(5) << h
                  << "  L2(rho)=" << std::scientific << std::setprecision(3) << err
                  << "  rate=" << std::fixed << std::setprecision(2) << rate << "\n";
        prevE = err;
        prevH = h;
    }
    std::cout << "  expected asymptotic rate: 2\n";

    std::cout << "\n[3] 2D temporal convergence on a fixed mesh (successive differences)\n";
    const int nxT = 80, nyT = 60;
    std::vector<double> dtList = {0.0010, 0.0005, 0.00025, 0.000125};
    std::vector<std::vector<double>> sols;
    for (double dt : dtList) sols.push_back(runEntropyDensity2D(nxT, nyT, dt, T));
    prevE = 0.0;
    double minTemporalRate = 1e300;
    for (size_t k = 0; k + 1 < sols.size(); ++k) {
        double err = densityDifference2D(sols[k], sols[k + 1], nxT, nyT);
        double rate = prevE > 0.0 ? std::log(prevE / err) / std::log(2.0) : 0.0;
        if (prevE > 0.0) minTemporalRate = std::min(minTemporalRate, rate);
        std::cout << "  dt=" << std::scientific << std::setprecision(3) << dtList[k]
                  << " -> " << dtList[k + 1]
                  << "  ||rho(dt)-rho(dt/2)||_L2=" << err
                  << "  rate=" << std::fixed << std::setprecision(2) << rate << "\n";
        prevE = err;
    }
    std::cout << "  expected asymptotic rate: 2\n";

    bool ok = true;
    ok = ok && std::isfinite(maxGcl) && maxGcl < 1e-10;
    ok = ok && std::isfinite(finalSpatialErr) && finalSpatialErr < 5e-5;
    ok = ok && std::isfinite(minSpatialRate) && minSpatialRate > 1.75;
    ok = ok && std::isfinite(minTemporalRate) && minTemporalRate > 1.85;
    std::cout << "\nVerification " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) {
        std::cerr << "Euler-ALE 2D convergence gate failed: maxGcl=" << maxGcl
                  << " finalSpatialErr=" << finalSpatialErr
                  << " minSpatialRate=" << minSpatialRate
                  << " minTemporalRate=" << minTemporalRate << "\n";
        return 2;
    }
    return 0;
}
