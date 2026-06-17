#include "EulerALE2D.h"
#include "FSIDiagnostics.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace euler_fsi;

namespace {

double clamp01(double x) {
    return std::max(0.0, std::min(1.0, x));
}

double minSpacing(const Grid2D& g) {
    double h = 1e300;
    for (int i = 0; i < g.nx; ++i) h = std::min(h, g.dx(i));
    for (int j = 0; j < g.ny; ++j) h = std::min(h, g.dy(j));
    return h;
}

Prim2D shockChannelIC(double x, double L) {
    if (x < 0.20 * L) return {1.0, 0.0, 0.0, 28.0};
    return {0.75, 0.0, 0.0, 0.8};
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

void showcaseColour(double s, unsigned char& r, unsigned char& g, unsigned char& b) {
    s = clamp01(s);
    const double stops[5][3] = {
        {12.0, 10.0, 32.0},
        {34.0, 70.0, 150.0},
        {30.0, 180.0, 190.0},
        {246.0, 190.0, 58.0},
        {255.0, 82.0, 38.0}
    };
    double x = s * 4.0;
    int k = std::max(0, std::min(3, static_cast<int>(std::floor(x))));
    double a = x - k;
    auto interp = [&](int c) {
        return static_cast<unsigned char>(std::lround((1.0 - a) * stops[k][c] + a * stops[k + 1][c]));
    };
    r = interp(0);
    g = interp(1);
    b = interp(2);
}

void writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img) {
    std::ofstream os(path, std::ios::binary);
    if (!os) {
        std::cerr << "showcase writePPM: cannot open " << path << "\n";
        return;
    }
    os << "P6\n" << W << " " << H << "\n255\n";
    os.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
}

void renderShowcase(const std::string& path, const EulerALE2D& solver,
                    const SpringPiston& piston, const std::vector<double>& pistonTrail,
                    int W, int H) {
    std::vector<unsigned char> img(static_cast<size_t>(W) * H * 3, 12);
    const Grid2D& g = solver.grid();
    double xmin = 0.0;
    double xmax = piston.maxX;
    double ymin = g.yFace.front();
    double ymax = g.yFace.back();
    double xspan = std::max(xmax - xmin, 1e-300);
    double yspan = std::max(ymax - ymin, 1e-300);

    auto setPixel = [&](int col, int row, unsigned char r, unsigned char gg, unsigned char b) {
        if (col < 0 || col >= W || row < 0 || row >= H) return;
        size_t k = (static_cast<size_t>(row) * W + col) * 3;
        img[k] = r;
        img[k + 1] = gg;
        img[k + 2] = b;
    };
    auto blendPixel = [&](int col, int row, unsigned char r, unsigned char gg,
                          unsigned char b, double alpha) {
        if (col < 0 || col >= W || row < 0 || row >= H) return;
        alpha = clamp01(alpha);
        if (alpha <= 0.0) return;
        size_t k = (static_cast<size_t>(row) * W + col) * 3;
        double beta = 1.0 - alpha;
        img[k] = static_cast<unsigned char>(std::lround(beta * img[k] + alpha * r));
        img[k + 1] = static_cast<unsigned char>(std::lround(beta * img[k + 1] + alpha * gg));
        img[k + 2] = static_cast<unsigned char>(std::lround(beta * img[k + 2] + alpha * b));
    };
    auto xToPx = [&](double x) { return (x - xmin) / xspan * (W - 1); };
    auto yToPx = [&](double y) { return (ymax - y) / yspan * (H - 1); };
    auto coverage1D = [](double c, double a, double b) {
        if (a > b) std::swap(a, b);
        return std::max(0.0, std::min(c + 0.5, b) - std::max(c - 0.5, a));
    };
    auto aaRect = [&](double x0, double y0, double x1, double y1,
                      unsigned char r, unsigned char gg, unsigned char b,
                      double opacity = 1.0) {
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);
        int c0 = std::max(0, static_cast<int>(std::floor(x0 - 1.0)));
        int c1 = std::min(W - 1, static_cast<int>(std::ceil(x1 + 1.0)));
        int r0 = std::max(0, static_cast<int>(std::floor(y0 - 1.0)));
        int r1 = std::min(H - 1, static_cast<int>(std::ceil(y1 + 1.0)));
        for (int row = r0; row <= r1; ++row) {
            double ay = coverage1D(row, y0, y1);
            for (int col = c0; col <= c1; ++col)
                blendPixel(col, row, r, gg, b, opacity * ay * coverage1D(col, x0, x1));
        }
    };
    auto aaLine = [&](double x0, double y0, double x1, double y1,
                      unsigned char r, unsigned char gg, unsigned char b,
                      double halfWidth = 0.9, double opacity = 1.0) {
        int c0 = std::max(0, static_cast<int>(std::floor(std::min(x0, x1) - halfWidth - 2.0)));
        int c1 = std::min(W - 1, static_cast<int>(std::ceil(std::max(x0, x1) + halfWidth + 2.0)));
        int r0 = std::max(0, static_cast<int>(std::floor(std::min(y0, y1) - halfWidth - 2.0)));
        int r1 = std::min(H - 1, static_cast<int>(std::ceil(std::max(y0, y1) + halfWidth + 2.0)));
        double dx = x1 - x0, dy = y1 - y0;
        double len2 = dx * dx + dy * dy;
        for (int row = r0; row <= r1; ++row)
            for (int col = c0; col <= c1; ++col) {
                double a = len2 > 1e-300 ? ((col - x0) * dx + (row - y0) * dy) / len2 : 0.0;
                a = clamp01(a);
                double px = x0 + a * dx, py = y0 + a * dy;
                double alpha = clamp01(halfWidth + 0.5 - std::hypot(col - px, row - py));
                blendPixel(col, row, r, gg, b, opacity * alpha);
            }
    };
    auto aaCircle = [&](double cx, double cy, double radius,
                        unsigned char r, unsigned char gg, unsigned char b,
                        double opacity = 1.0) {
        int c0 = std::max(0, static_cast<int>(std::floor(cx - radius - 1.0)));
        int c1 = std::min(W - 1, static_cast<int>(std::ceil(cx + radius + 1.0)));
        int r0 = std::max(0, static_cast<int>(std::floor(cy - radius - 1.0)));
        int r1 = std::min(H - 1, static_cast<int>(std::ceil(cy + radius + 1.0)));
        for (int row = r0; row <= r1; ++row)
            for (int col = c0; col <= c1; ++col) {
                double alpha = clamp01(radius + 0.5 - std::hypot(col - cx, row - cy));
                blendPixel(col, row, r, gg, b, opacity * alpha);
            }
    };

    for (int row = 0; row < H; ++row) {
        double y = ymax - yspan * (row + 0.5) / H;
        auto jy = std::upper_bound(g.yFace.begin(), g.yFace.end(), y);
        int j = std::max(0, std::min(solver.ny() - 1, static_cast<int>(jy - g.yFace.begin()) - 1));
        for (int col = 0; col < W; ++col) {
            double x = xmin + xspan * (col + 0.5) / W;
            unsigned char r = 20, gg = 22, b = 32;
            if (x >= g.xFace.front() && x <= g.xFace.back() && y >= ymin && y <= ymax) {
                auto ix = std::upper_bound(g.xFace.begin(), g.xFace.end(), x);
                int i = std::max(0, std::min(solver.nx() - 1, static_cast<int>(ix - g.xFace.begin()) - 1));
                Prim2D q = solver.cellPrimitive(i, j);
                double rhoS = clamp01((q.rho - 0.55) / 5.2);
                double pS = clamp01((q.p - 0.6) / 28.0);
                double s = clamp01(0.72 * rhoS + 0.28 * pS);
                showcaseColour(std::pow(s, 0.72), r, gg, b);
            }
            setPixel(col, row, r, gg, b);
        }
    }

    double left = xToPx(g.xFace.front());
    double wall = xToPx(g.xFace.back());
    double rest = xToPx(piston.xRest);
    double top = yToPx(g.yFace.back());
    double bottom = yToPx(g.yFace.front());
    if (top > bottom) std::swap(top, bottom);
    aaRect(left - 1.0, top, left + 1.0, bottom, 230, 235, 244, 0.65);
    aaRect(left, top - 1.0, wall, top + 1.0, 230, 235, 244, 0.6);
    aaRect(left, bottom - 1.0, wall, bottom + 1.0, 230, 235, 244, 0.6);

    double slabHalf = 8.5;
    aaRect(wall - slabHalf - 4.0, top - 6.0, wall + slabHalf + 4.0, bottom + 6.0,
           235, 242, 255, 0.72);
    aaRect(wall - slabHalf, top - 6.0, wall + slabHalf, bottom + 6.0,
           12, 16, 26, 1.0);
    aaRect(rest - 0.8, top, rest + 0.8, bottom, 160, 170, 190, 0.36);

    if (pistonTrail.size() > 1) {
        double y = std::max(5.0, top - H / 11.0);
        for (size_t k = 1; k < pistonTrail.size(); ++k) {
            double a = static_cast<double>(k) / static_cast<double>(pistonTrail.size() - 1);
            double x0 = xToPx(pistonTrail[k - 1]);
            double x1 = xToPx(pistonTrail[k]);
            aaLine(x0, y, x1, y, 255, 70, 50, 1.0, 0.18 + 0.52 * a);
        }
    }
    aaCircle(wall, std::max(5.0, top - H / 11.0), 5.0, 255, 66, 48, 1.0);

    double springY = std::min(H - 5.0, bottom + H / 9.0);
    double start = wall + slabHalf + 5.0;
    double anchor = xToPx(xmax);
    double end = std::max(start, anchor - 10.0);
    double amp = std::max(4.0, H / 56.0);
    double prevX = start, prevY = springY;
    for (int k = 1; k <= 18; ++k) {
        double x = start + (end - start) * k / 18.0;
        double y = springY + ((k % 2) ? amp : -amp);
        aaLine(prevX, prevY, x, y, 230, 235, 244, 1.15, 0.9);
        prevX = x;
        prevY = y;
    }
    aaRect(anchor - 2.4, springY - 24.0, anchor + 2.4, springY + 24.0, 230, 235, 244, 0.9);
    aaLine(wall + slabHalf, springY, start, springY, 230, 235, 244, 1.0, 0.8);

    writePPM(path, W, H, img);
}

} // namespace

int main(int argc, char** argv) {
    const int nx = 360;
    const int ny = 72;
    const double Ly = 0.32;
    const double tEnd = 0.95;
    const double cfl = 0.40;
    const int frames = 240;
    const int W = 1280, H = 360;
    std::string outDir = "out/euler_fsi_piston2d_showcase_frames";
    if (argc > 1) outDir = argv[1];

    std::filesystem::create_directories("out");
    std::filesystem::remove_all(outDir);
    std::filesystem::create_directories(outDir);

    SpringPiston piston;
    piston.x = 1.0;
    piston.v = 0.0;
    piston.xRest = 1.0;
    piston.mass = 0.75;
    piston.stiffness = 10.0;
    piston.damping = 0.85;
    piston.area = Ly;
    piston.externalPressure = 0.8;
    piston.minX = 0.62;
    piston.maxX = 1.85;
    const double initialPistonX = piston.x;

    EulerALE2D solver(nx, ny);
    Grid2D g0 = Grid2D::affinePiston(nx, ny, piston.x, Ly, piston.v);
    solver.setCellAverages([&](double x, double) { return shockChannelIC(x, piston.x); }, g0);
    solver.setGrid(g0, 0.0);
    const auto initialTotals = solver.conservedTotals();
    FSIEnergyBudget energyBudget;
    energyBudget.reset(initialTotals[3], piston);

    std::ofstream diag("out/euler_fsi_piston2d_showcase_diagnostics.csv");
    diag << "step,t,piston_x,piston_v,p_wall,rho_min,p_min,mass,momentum_x,momentum_y,gas_energy,"
            "piston_ke,spring_energy,external_pressure_potential,"
         << FSIEnergyBudget::csvHeader() << "\n";

    std::vector<double> trail;
    auto writeFrame = [&](int frame) {
        trail.push_back(piston.x);
        if (trail.size() > 96) trail.erase(trail.begin());
        char path[512];
        std::snprintf(path, sizeof(path), "%s/frame_%05d.ppm", outDir.c_str(), frame);
        renderShowcase(path, solver, piston, trail, W, H);
    };

    std::cout << "2D Euler-ALE FSI showcase: shock-loaded spring piston\n";
    std::cout << "  cells=" << nx << "x" << ny << "  t_end=" << tEnd
              << "  frames=" << frames << "\n";
    std::cout << "  output frames: " << outDir << "\n";

    auto clock0 = std::chrono::high_resolution_clock::now();
    double t = 0.0, nextFrame = 0.0, frameDt = tEnd / frames;
    int step = 0, frame = 0;
    double minRhoSeen = solver.minDensity();
    double minPSeen = solver.minPressure();
    double maxMassDrift = 0.0;
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
            std::cerr << "ERROR: non-positive showcase state at step " << step << ", t=" << t << "\n";
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
            if (frame % 20 == 0 || t >= tEnd - 1e-14) {
                std::cout << "  frame " << std::setw(4) << frame
                          << "  t=" << std::fixed << std::setprecision(4) << t
                          << "  piston x=" << piston.x
                          << "  v=" << piston.v
                          << "  p_wall=" << pWall << "\n";
            }
        }
    }

    renderShowcase("out/euler_fsi_piston2d_showcase_final.ppm", solver, piston, trail, W, H);
    double elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock0).count();
    std::cout << "\nDone. steps=" << step << " frames=" << frame
              << " wall=" << std::fixed << std::setprecision(2) << elapsed << "s\n";
    std::cout << "  final still: out/euler_fsi_piston2d_showcase_final.ppm\n";
    std::cout << "  diagnostics: out/euler_fsi_piston2d_showcase_diagnostics.csv\n";
    std::cout << "  ffmpeg -y -framerate 30 -i " << outDir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_piston2d_showcase.mp4\n";

    bool ok = true;
    ok = ok && std::isfinite(piston.x) && std::isfinite(piston.v);
    ok = ok && std::abs(piston.x - initialPistonX) > 0.12;
    ok = ok && minRhoSeen > 0.0 && minPSeen > 0.0;
    ok = ok && maxMassDrift < 1e-8;
    ok = ok && frame == frames + 1;
    ok = ok && std::filesystem::exists("out/euler_fsi_piston2d_showcase_final.ppm");
    std::cout << "  verification: min_rho=" << std::scientific << minRhoSeen
              << " min_p=" << minPSeen
              << " max_mass_drift=" << maxMassDrift
              << " max_gas_work_drift=" << energyBudget.maxGasWorkDrift()
              << " max_coupled_energy_drift=" << energyBudget.maxCoupledEnergyDrift()
              << " piston_dx=" << std::abs(piston.x - initialPistonX)
              << " frames=" << std::fixed << frame << "\n";
    std::cout << "Verification " << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 2;
}
