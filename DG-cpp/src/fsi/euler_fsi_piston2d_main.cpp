#include "EulerALE2D.h"
#include "FSIDiagnostics.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace euler_fsi;

namespace {

double minSpacing(const Grid2D& g) {
    double h = 1e300;
    for (int i = 0; i < g.nx; ++i) h = std::min(h, g.dx(i));
    for (int j = 0; j < g.ny; ++j) h = std::min(h, g.dy(j));
    return h;
}

Prim2D shockChannelIC(double x, double L) {
    if (x < 0.28 * L) return {1.0, 0.0, 0.0, 8.0};
    return {1.0, 0.0, 0.0, 1.0};
}

double rightWallPressure(const EulerALE2D& solver) {
    const Grid2D& g = solver.grid();
    double num = 0.0, den = 0.0;
    int i = solver.nx() - 1;
    for (int j = 0; j < solver.ny(); ++j) {
        double w = g.dy(j);
        num += w * solver.cellPrimitive(i, j).p;
        den += w;
    }
    return num / std::max(den, 1e-300);
}

} // namespace

int main(int argc, char** argv) {
    const int nx = 240;
    const int ny = 48;
    const double Ly = 0.25;
    const double tEnd = 0.50;
    const double cfl = 0.42;
    const int frames = 150;
    const int W = 1024, H = 256;
    std::string outDir = "out/euler_fsi_piston2d_frames";
    if (argc > 1) outDir = argv[1];

    std::filesystem::create_directories("out");
    std::filesystem::remove_all(outDir);
    std::filesystem::create_directories(outDir);

    SpringPiston piston;
    piston.x = 1.0;
    piston.v = 0.0;
    piston.xRest = 1.0;
    piston.mass = 7.5;
    piston.stiffness = 180.0;
    piston.damping = 6.0;
    piston.area = Ly;
    piston.externalPressure = 1.0;
    piston.minX = 0.72;
    piston.maxX = 1.26;
    const double initialPistonX = piston.x;

    EulerALE2D solver(nx, ny);
    Grid2D g0 = Grid2D::affinePiston(nx, ny, piston.x, Ly, piston.v);
    solver.setCellAverages([&](double x, double) { return shockChannelIC(x, piston.x); }, g0);
    solver.setGrid(g0, 0.0);
    const auto initialTotals = solver.conservedTotals();
    FSIEnergyBudget energyBudget;
    energyBudget.reset(initialTotals[3], piston);

    double minRhoSeen = solver.minDensity();
    double minPSeen = solver.minPressure();
    double maxMassDrift = 0.0;

    std::ofstream diag("out/euler_fsi_piston2d_diagnostics.csv");
    diag << "step,t,piston_x,piston_v,p_wall,rho_min,p_min,mass,momentum_x,momentum_y,gas_energy,"
            "piston_ke,spring_energy,external_pressure_potential,"
         << FSIEnergyBudget::csvHeader() << "\n";

    auto writeFrame = [&](int frame) {
        char path[512];
        std::snprintf(path, sizeof(path), "%s/frame_%05d.ppm", outDir.c_str(), frame);
        writePistonFramePPM2D(path, solver, piston, W, H, 0.65, 3.2);
    };

    std::cout << "2D Euler-ALE FSI shock channel + spring piston demo\n";
    std::cout << "  cells=" << nx << "x" << ny << "  t_end=" << tEnd
              << "  frames=" << frames << "\n";
    std::cout << "  output frames: " << outDir << "\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    double t = 0.0, nextFrame = 0.0, frameDt = tEnd / frames;
    int step = 0, frame = 0;
    writeFrame(frame++);
    nextFrame += frameDt;

    while (t < tEnd - 1e-14) {
        double pWall = rightWallPressure(solver);
        double dt = cfl * minSpacing(solver.grid()) /
                    (solver.maxWaveSpeed() + std::abs(piston.v) + 1e-12);
        dt = std::min(dt, tEnd - t);
        dt = std::min(dt, std::max(1e-8, nextFrame - t));

        double xOld = piston.x;
        double vOld = piston.v;
        piston.advanceExplicit(dt, pWall);
        double xNew = piston.x;
        double wallVel = (xNew - xOld) / dt;

        auto gridFn = [=](double tau) {
            double a = (tau - t) / dt;
            a = std::min(1.0, std::max(0.0, a));
            return Grid2D::affinePiston(nx, ny, xOld + (xNew - xOld) * a, Ly, wallVel);
        };
        auto bcFn = [](double) {
            BoundaryCondition2D bc;
            bc.left = BoundaryKind::FixedWall;
            bc.right = BoundaryKind::MovingWall;
            bc.bottom = BoundaryKind::FixedWall;
            bc.top = BoundaryKind::FixedWall;
            return bc;
        };
        if (!solver.step(dt, gridFn, bcFn, SlopeMode::Minmod)) {
            std::cerr << "ERROR: non-positive 2D state at step " << step << ", t=" << t << "\n";
            return 1;
        }
        t += dt;
        ++step;

        auto tot = solver.conservedTotals();
        FSIEnergySnapshot energy = energyBudget.advance(dt, tot[3], piston, pWall, xOld, vOld);
        maxMassDrift = std::max(maxMassDrift, std::abs(tot[0] - initialTotals[0]));
        minRhoSeen = std::min(minRhoSeen, solver.minDensity());
        minPSeen = std::min(minPSeen, solver.minPressure());

        diag << step << "," << std::setprecision(16) << t << "," << piston.x << "," << piston.v
             << "," << pWall << "," << solver.minDensity() << "," << solver.minPressure()
             << "," << tot[0] << "," << tot[1] << "," << tot[2] << "," << tot[3]
             << "," << piston.kineticEnergy() << "," << piston.springEnergy()
             << "," << piston.externalPressurePotential() << "," << energyBudget.csvValues(energy) << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++);
            nextFrame += frameDt;
            std::cout << "  frame " << std::setw(4) << frame
                      << "  t=" << std::fixed << std::setprecision(4) << t
                      << "  piston x=" << piston.x
                      << "  v=" << piston.v
                      << "  p_wall=" << pWall << "\n";
        }
    }

    writePistonFramePPM2D("out/euler_fsi_piston2d_final.ppm", solver, piston, W, H, 0.65, 3.2);
    double wall = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << "\nDone. steps=" << step << " frames=" << frame
              << " wall=" << std::fixed << std::setprecision(2) << wall << "s\n";
    std::cout << "  final still: out/euler_fsi_piston2d_final.ppm\n";
    std::cout << "  diagnostics: out/euler_fsi_piston2d_diagnostics.csv\n";
    std::cout << "  ffmpeg -y -framerate 30 -i " << outDir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_piston2d.mp4\n";

    bool ok = true;
    ok = ok && std::isfinite(piston.x) && std::isfinite(piston.v);
    ok = ok && std::abs(piston.x - initialPistonX) > 2e-3;
    ok = ok && minRhoSeen > 0.0 && minPSeen > 0.0;
    ok = ok && maxMassDrift < 1e-8;
    ok = ok && energyBudget.maxGasWorkDrift() < 1e-3;
    ok = ok && energyBudget.maxCoupledEnergyDrift() < 1e-3;
    ok = ok && frame == frames + 1;
    ok = ok && std::filesystem::exists("out/euler_fsi_piston2d_final.ppm");
    std::cout << "  verification: min_rho=" << std::scientific << minRhoSeen
              << " min_p=" << minPSeen
              << " max_mass_drift=" << maxMassDrift
              << " max_gas_work_drift=" << energyBudget.maxGasWorkDrift()
              << " max_coupled_energy_drift=" << energyBudget.maxCoupledEnergyDrift()
              << " piston_dx=" << std::abs(piston.x - initialPistonX)
              << " frames=" << std::fixed << frame << "\n";
    std::cout << "Verification " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) return 2;
    return 0;
}
