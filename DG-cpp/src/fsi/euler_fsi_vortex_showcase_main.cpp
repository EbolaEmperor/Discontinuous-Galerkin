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

constexpr double PI = 3.141592653589793238462643383279502884;

double clamp01(double x) {
    return std::max(0.0, std::min(1.0, x));
}

double minSpacing(const Grid2D& g) {
    double h = 1e300;
    for (int i = 0; i < g.nx; ++i) h = std::min(h, g.dx(i));
    for (int j = 0; j < g.ny; ++j) h = std::min(h, g.dy(j));
    return h;
}

Prim2D shockRmiIC(double x, double y, double Ly) {
    if (x < 0.10) return {1.2, 0.0, 0.0, 28.0};

    double iface = 0.42
                 + 0.038 * std::sin(2.0 * PI * y / Ly)
                 + 0.014 * std::sin(6.0 * PI * y / Ly + 0.45);
    double h = 0.5 * (1.0 + std::tanh((x - iface) / 0.010));
    double rho = (1.0 - h) * 3.2 + h * 0.55;
    double vSeed = 0.018 * std::sin(4.0 * PI * y / Ly)
                 * std::exp(-std::pow((x - iface) / 0.045, 2.0));
    return {rho, 0.0, vSeed, 0.80};
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

std::vector<double> computeVorticity(const EulerALE2D& solver, double& maxAbs) {
    const Grid2D& g = solver.grid();
    std::vector<double> omega(static_cast<size_t>(solver.nx()) * solver.ny(), 0.0);
    maxAbs = 0.0;
    for (int j = 0; j < solver.ny(); ++j) {
        int jm = std::max(0, j - 1);
        int jp = std::min(solver.ny() - 1, j + 1);
        double ym = g.centerY(jm);
        double yp = g.centerY(jp);
        for (int i = 0; i < solver.nx(); ++i) {
            int im = std::max(0, i - 1);
            int ip = std::min(solver.nx() - 1, i + 1);
            double xm = g.centerX(im);
            double xp = g.centerX(ip);
            Prim2D qL = solver.cellPrimitive(im, j);
            Prim2D qR = solver.cellPrimitive(ip, j);
            Prim2D qB = solver.cellPrimitive(i, jm);
            Prim2D qT = solver.cellPrimitive(i, jp);
            double dvdx = (qR.v - qL.v) / std::max(xp - xm, 1e-300);
            double dudy = (qT.u - qB.u) / std::max(yp - ym, 1e-300);
            double w = dvdx - dudy;
            omega[solver.index(i, j)] = w;
            maxAbs = std::max(maxAbs, std::abs(w));
        }
    }
    return omega;
}

void falseColour(double s, unsigned char& r, unsigned char& g, unsigned char& b) {
    s = clamp01(s);
    const double c[5][3] = {
        {8.0, 9.0, 24.0},
        {25.0, 60.0, 150.0},
        {20.0, 180.0, 190.0},
        {248.0, 200.0, 74.0},
        {255.0, 78.0, 42.0}
    };
    double x = 4.0 * s;
    int k = std::max(0, std::min(3, static_cast<int>(std::floor(x))));
    double a = x - k;
    auto mix = [&](int ch) {
        return static_cast<unsigned char>(std::lround((1.0 - a) * c[k][ch] + a * c[k + 1][ch]));
    };
    r = mix(0);
    g = mix(1);
    b = mix(2);
}

void writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream os(path, std::ios::binary);
    if (!os) {
        std::cerr << "vortex showcase: cannot open " << path << "\n";
        return;
    }
    os << "P6\n" << W << " " << H << "\n255\n";
    os.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
}

void renderVortexShowcase(const std::string& path, const EulerALE2D& solver,
                          const SpringPiston& piston, const std::vector<double>& pistonTrail,
                          int W, int H) {
    std::vector<unsigned char> img(static_cast<size_t>(W) * H * 3, 10);
    const Grid2D& g = solver.grid();
    double xmin = 0.0;
    double xmax = piston.maxX;
    double ymin = g.yFace.front();
    double ymax = g.yFace.back();
    double xspan = std::max(xmax - xmin, 1e-300);
    double yspan = std::max(ymax - ymin, 1e-300);
    double maxOmega = 0.0;
    std::vector<double> omega = computeVorticity(solver, maxOmega);
    double omegaScale = std::max(18.0, 0.38 * maxOmega);

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
            unsigned char r = 14, gg = 16, b = 28;
            if (x >= g.xFace.front() && x <= g.xFace.back() && y >= ymin && y <= ymax) {
                auto ix = std::upper_bound(g.xFace.begin(), g.xFace.end(), x);
                int i = std::max(0, std::min(solver.nx() - 1, static_cast<int>(ix - g.xFace.begin()) - 1));
                Prim2D q = solver.cellPrimitive(i, j);
                double shock = clamp01((q.p - 0.75) / 22.0);
                double rho = clamp01((q.rho - 0.35) / 3.1);
                falseColour(std::pow(clamp01(0.52 * rho + 0.48 * shock), 0.70), r, gg, b);
                double w = omega[solver.index(i, j)];
                double vort = std::sqrt(clamp01(std::abs(w) / omegaScale));
                if (vort > 0.02) {
                    unsigned char vr = (w >= 0.0) ? 255 : 40;
                    unsigned char vg = (w >= 0.0) ? 70 : 220;
                    unsigned char vb = (w >= 0.0) ? 45 : 255;
                    double alpha = 0.82 * vort;
                    r = static_cast<unsigned char>(std::lround((1.0 - alpha) * r + alpha * vr));
                    gg = static_cast<unsigned char>(std::lround((1.0 - alpha) * gg + alpha * vg));
                    b = static_cast<unsigned char>(std::lround((1.0 - alpha) * b + alpha * vb));
                }
            }
            setPixel(col, row, r, gg, b);
        }
    }

    double left = xToPx(g.xFace.front());
    double wall = xToPx(g.xFace.back());
    double top = yToPx(g.yFace.back());
    double bottom = yToPx(g.yFace.front());
    if (top > bottom) std::swap(top, bottom);
    aaRect(left - 1.0, top, left + 1.0, bottom, 235, 240, 250, 0.62);
    aaRect(left, top - 1.0, wall, top + 1.0, 235, 240, 250, 0.56);
    aaRect(left, bottom - 1.0, wall, bottom + 1.0, 235, 240, 250, 0.56);

    double slabHalf = 8.0;
    aaRect(wall - slabHalf - 4.0, top - 5.0, wall + slabHalf + 4.0, bottom + 5.0,
           225, 235, 255, 0.70);
    aaRect(wall - slabHalf, top - 5.0, wall + slabHalf, bottom + 5.0,
           10, 14, 24, 1.0);

    if (pistonTrail.size() > 1) {
        double y = std::max(6.0, top - H / 10.0);
        for (size_t k = 1; k < pistonTrail.size(); ++k) {
            double a = static_cast<double>(k) / static_cast<double>(pistonTrail.size() - 1);
            aaLine(xToPx(pistonTrail[k - 1]), y, xToPx(pistonTrail[k]), y,
                   255, 64, 46, 1.0, 0.18 + 0.55 * a);
        }
    }
    aaCircle(wall, std::max(6.0, top - H / 10.0), 5.2, 255, 66, 48, 1.0);

    double springY = std::min(H - 6.0, bottom + H / 9.0);
    double anchor = xToPx(xmax);
    double start = wall + slabHalf + 5.0;
    double end = std::max(start, anchor - 10.0);
    double amp = std::max(4.0, H / 58.0);
    double px = start, py = springY;
    for (int k = 1; k <= 18; ++k) {
        double x = start + (end - start) * k / 18.0;
        double y = springY + ((k % 2) ? amp : -amp);
        aaLine(px, py, x, y, 230, 235, 246, 1.1, 0.88);
        px = x;
        py = y;
    }
    aaRect(anchor - 2.2, springY - 23.0, anchor + 2.2, springY + 23.0, 230, 235, 246, 0.88);
    aaLine(wall + slabHalf, springY, start, springY, 230, 235, 246, 1.0, 0.78);

    writePPM(path, W, H, img);
}

} // namespace

int main(int argc, char** argv) {
    const int nx = 420;
    const int ny = 126;
    const double Ly = 0.36;
    const double tEnd = 0.76;
    const double cfl = 0.36;
    const int frames = 210;
    const int W = 1280, H = 420;
    std::string outDir = "out/euler_fsi_vortex_showcase_frames";
    if (argc > 1) outDir = argv[1];

    std::filesystem::create_directories("out");
    std::filesystem::remove_all(outDir);
    std::filesystem::create_directories(outDir);

    SpringPiston piston;
    piston.x = 1.0;
    piston.v = 0.0;
    piston.xRest = 1.0;
    piston.mass = 0.70;
    piston.stiffness = 13.0;
    piston.damping = 0.70;
    piston.area = Ly;
    piston.externalPressure = 0.80;
    piston.minX = 0.58;
    piston.maxX = 1.70;
    const double initialPistonX = piston.x;

    EulerALE2D solver(nx, ny);
    Grid2D g0 = Grid2D::affinePiston(nx, ny, piston.x, Ly, piston.v);
    solver.setCellAverages([&](double x, double y) { return shockRmiIC(x, y, Ly); }, g0);
    solver.setGrid(g0, 0.0);
    const auto initialTotals = solver.conservedTotals();
    FSIEnergyBudget energyBudget;
    energyBudget.reset(initialTotals[3], piston);

    std::ofstream diag("out/euler_fsi_vortex_showcase_diagnostics.csv");
    diag << "step,t,piston_x,piston_v,p_wall,rho_min,p_min,mass,momentum_x,momentum_y,gas_energy,"
            "max_abs_vorticity,piston_ke,spring_energy,external_pressure_potential,"
         << FSIEnergyBudget::csvHeader() << "\n";

    double minRhoSeen = solver.minDensity();
    double minPSeen = solver.minPressure();
    double maxMassDrift = 0.0;
    double maxVorticitySeen = 0.0;
    std::vector<double> trail;
    auto writeFrame = [&](int frame) {
        trail.push_back(piston.x);
        if (trail.size() > 110) trail.erase(trail.begin());
        char path[512];
        std::snprintf(path, sizeof(path), "%s/frame_%05d.ppm", outDir.c_str(), frame);
        renderVortexShowcase(path, solver, piston, trail, W, H);
    };

    std::cout << "2D Euler-ALE FSI vortex showcase: shock + perturbed density interface + spring piston\n";
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
            std::cerr << "ERROR: non-positive vortex showcase state at step " << step << ", t=" << t << "\n";
            return 1;
        }
        t += dt;
        ++step;

        auto tot = solver.conservedTotals();
        FSIEnergySnapshot energy = energyBudget.advance(dt, tot[3], piston, pWall, xOld, vOld);
        maxMassDrift = std::max(maxMassDrift, std::abs(tot[0] - initialTotals[0]));
        minRhoSeen = std::min(minRhoSeen, solver.minDensity());
        minPSeen = std::min(minPSeen, solver.minPressure());
        double maxW = 0.0;
        if (step % 12 == 0 || t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            std::vector<double> omega = computeVorticity(solver, maxW);
            (void)omega;
            maxVorticitySeen = std::max(maxVorticitySeen, maxW);
        }
        diag << step << "," << std::setprecision(16) << t << "," << piston.x << "," << piston.v
             << "," << pWall << "," << solver.minDensity() << "," << solver.minPressure()
             << "," << tot[0] << "," << tot[1] << "," << tot[2] << "," << tot[3]
             << "," << maxVorticitySeen
             << "," << piston.kineticEnergy() << "," << piston.springEnergy()
             << "," << piston.externalPressurePotential() << "," << energyBudget.csvValues(energy) << "\n";

        if (t >= nextFrame - 1e-14 || t >= tEnd - 1e-14) {
            writeFrame(frame++);
            nextFrame += frameDt;
            if (frame % 18 == 0 || t >= tEnd - 1e-14) {
                std::cout << "  frame " << std::setw(4) << frame
                          << "  t=" << std::fixed << std::setprecision(4) << t
                          << "  piston x=" << piston.x
                          << "  v=" << piston.v
                          << "  p_wall=" << pWall
                          << "  max|omega|=" << std::scientific << maxVorticitySeen << std::fixed << "\n";
            }
        }
    }

    renderVortexShowcase("out/euler_fsi_vortex_showcase_final.ppm", solver, piston, trail, W, H);
    double elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << "\nDone. steps=" << step << " frames=" << frame
              << " wall=" << std::fixed << std::setprecision(2) << elapsed << "s\n";
    std::cout << "  final still: out/euler_fsi_vortex_showcase_final.ppm\n";
    std::cout << "  diagnostics: out/euler_fsi_vortex_showcase_diagnostics.csv\n";
    std::cout << "  ffmpeg -y -framerate 30 -i " << outDir
              << "/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_vortex_showcase.mp4\n";

    bool ok = true;
    ok = ok && std::isfinite(piston.x) && std::isfinite(piston.v);
    ok = ok && std::abs(piston.x - initialPistonX) > 0.06;
    ok = ok && maxVorticitySeen > 10.0;
    ok = ok && minRhoSeen > 0.0 && minPSeen > 0.0;
    ok = ok && maxMassDrift < 1e-8;
    ok = ok && frame == frames + 1;
    ok = ok && std::filesystem::exists("out/euler_fsi_vortex_showcase_final.ppm");
    std::cout << "  verification: min_rho=" << std::scientific << minRhoSeen
              << " min_p=" << minPSeen
              << " max_mass_drift=" << maxMassDrift
              << " max_abs_vorticity=" << maxVorticitySeen
              << " max_gas_work_drift=" << energyBudget.maxGasWorkDrift()
              << " max_coupled_energy_drift=" << energyBudget.maxCoupledEnergyDrift()
              << " piston_dx=" << std::abs(piston.x - initialPistonX)
              << " frames=" << std::fixed << frame << "\n";
    std::cout << "Verification " << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 2;
}
