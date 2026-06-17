#include "EulerALE1D.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstdlib>

namespace euler_fsi {

namespace {

constexpr double RHO_FLOOR = 1e-12;
constexpr double P_FLOOR = 1e-12;
constexpr double PI = 3.141592653589793238462643383279502884;

double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

State1D minmodState(const State1D& a, const State1D& b) {
    return {minmod(a.rho, b.rho), minmod(a.mom, b.mom), minmod(a.ene, b.ene)};
}

double clamp01(double x) {
    return std::min(1.0, std::max(0.0, x));
}

void colour(double s, unsigned char& r, unsigned char& g, unsigned char& b) {
    s = clamp01(s);
    const double c0[3] = {12, 7, 30};
    const double c1[3] = {50, 90, 170};
    const double c2[3] = {40, 170, 145};
    const double c3[3] = {245, 205, 70};
    const double* A = c0;
    const double* B = c1;
    double t = s * 3.0;
    if (t >= 2.0) { A = c2; B = c3; t -= 2.0; }
    else if (t >= 1.0) { A = c1; B = c2; t -= 1.0; }
    r = static_cast<unsigned char>(std::lround(A[0] + t * (B[0] - A[0])));
    g = static_cast<unsigned char>(std::lround(A[1] + t * (B[1] - A[1])));
    b = static_cast<unsigned char>(std::lround(A[2] + t * (B[2] - A[2])));
}

void ensureParent(const std::string& path) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
}

void writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img) {
    ensureParent(path);
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        std::cerr << "writePPM: cannot open " << path << "\n";
        return;
    }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
}

Grid1D withSweptFaceVelocity(const Grid1D& grid0, const Grid1D& grid1, double dt) {
    Grid1D grid = grid0;
    if (grid0.xFace.size() != grid1.xFace.size() || dt <= 0.0) return grid;
    grid.wFace.resize(grid.xFace.size());
    for (size_t i = 0; i < grid.xFace.size(); ++i) {
        grid.wFace[i] = (grid1.xFace[i] - grid0.xFace[i]) / dt;
    }
    return grid;
}

} // namespace

State1D operator+(const State1D& a, const State1D& b) {
    return {a.rho + b.rho, a.mom + b.mom, a.ene + b.ene};
}
State1D operator-(const State1D& a, const State1D& b) {
    return {a.rho - b.rho, a.mom - b.mom, a.ene - b.ene};
}
State1D operator*(double s, const State1D& a) {
    return {s * a.rho, s * a.mom, s * a.ene};
}
State1D operator*(const State1D& a, double s) { return s * a; }
State1D operator/(const State1D& a, double s) {
    return {a.rho / s, a.mom / s, a.ene / s};
}

double pressure(const State1D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    double u = U.mom / rho;
    return (GAMMA - 1.0) * (U.ene - 0.5 * rho * u * u);
}

double soundSpeed(const State1D& U) {
    return std::sqrt(GAMMA * std::max(pressure(U), P_FLOOR) / std::max(U.rho, RHO_FLOOR));
}

Prim1D consToPrim(const State1D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    return {rho, U.mom / rho, std::max(pressure(U), P_FLOOR)};
}

State1D primToCons(const Prim1D& q) {
    double rho = std::max(q.rho, RHO_FLOOR);
    double p = std::max(q.p, P_FLOOR);
    return {rho, rho * q.u, p / (GAMMA - 1.0) + 0.5 * rho * q.u * q.u};
}

State1D flux(const State1D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    double u = U.mom / rho;
    double p = pressure(U);
    return {U.mom, U.mom * u + p, (U.ene + p) * u};
}

State1D aleRusanovFlux(const State1D& UL, const State1D& UR, double wFace) {
    State1D FL = flux(UL) - wFace * UL;
    State1D FR = flux(UR) - wFace * UR;
    double uL = UL.mom / std::max(UL.rho, RHO_FLOOR);
    double uR = UR.mom / std::max(UR.rho, RHO_FLOOR);
    double lam = std::max(std::abs(uL - wFace) + soundSpeed(UL),
                          std::abs(uR - wFace) + soundSpeed(UR));
    return 0.5 * (FL + FR) - 0.5 * lam * (UR - UL);
}

State1D movingWallFlux(const State1D& Uinside, double wallVelocity) {
    double p = std::max(pressure(Uinside), P_FLOOR);
    return {0.0, p, p * wallVelocity};
}

Grid1D Grid1D::affinePiston(int n, double length, double velocity) {
    Grid1D g;
    g.xFace.resize(n + 1);
    g.wFace.resize(n + 1);
    for (int i = 0; i <= n; ++i) {
        double xi = static_cast<double>(i) / n;
        g.xFace[i] = length * xi;
        g.wFace[i] = velocity * xi;
    }
    return g;
}

Grid1D Grid1D::oscillatory(int n, double length, double t, double amplitude, double omega) {
    Grid1D g;
    g.xFace.resize(n + 1);
    g.wFace.resize(n + 1);
    for (int i = 0; i <= n; ++i) {
        double eta = static_cast<double>(i) / n;
        double shape = std::sin(2.0 * PI * eta);
        g.xFace[i] = length * (eta + amplitude * shape * std::sin(omega * t));
        g.wFace[i] = length * amplitude * shape * omega * std::cos(omega * t);
    }
    return g;
}

EulerALE1D::EulerALE1D(int cells) : n_(cells), q_(cells) {}

void EulerALE1D::setGrid(const Grid1D& grid, double time) {
    grid_ = grid;
    time_ = time;
}

void EulerALE1D::setCellAverages(const std::function<Prim1D(double)>& prim, const Grid1D& grid) {
    grid_ = grid;
    q_.assign(n_, {});
    const double gp[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    for (int i = 0; i < n_; ++i) {
        double xa = grid.xFace[i], xb = grid.xFace[i + 1], dx = xb - xa;
        State1D integral{};
        for (double s : gp) {
            double x = 0.5 * (xa + xb) + 0.5 * dx * s;
            integral = integral + (0.5 * dx) * primToCons(prim(x));
        }
        q_[i] = integral;
    }
}

void EulerALE1D::setUniform(const Prim1D& prim, const Grid1D& grid) {
    setCellAverages([&](double) { return prim; }, grid);
}

State1D EulerALE1D::cellAverage(int i) const {
    return q_[i] / grid_.dx(i);
}

std::vector<State1D> EulerALE1D::reconstructSlopes(const std::vector<State1D>& U,
                                                   const Grid1D& grid,
                                                   const BoundaryCondition1D& bc,
                                                   SlopeMode mode) const {
    std::vector<State1D> slope(n_);
    if (mode == SlopeMode::FirstOrder) return slope;
    bool periodic = (bc.left == BoundaryKind::Periodic && bc.right == BoundaryKind::Periodic);
    for (int i = 0; i < n_; ++i) {
        int im = i - 1, ip = i + 1;
        if (periodic) {
            im = (i + n_ - 1) % n_;
            ip = (i + 1) % n_;
        } else {
            im = std::max(0, im);
            ip = std::min(n_ - 1, ip);
        }
        if (!periodic && (im == i || ip == i)) {
            slope[i] = {};
            continue;
        }
        double xim = grid.center(im), xi = grid.center(i), xip = grid.center(ip);
        if (periodic) {
            double L = grid.length();
            if (im > i) xim -= L;
            if (ip < i) xip += L;
        }
        State1D left = (U[i] - U[im]) / std::max(xi - xim, 1e-300);
        State1D right = (U[ip] - U[i]) / std::max(xip - xi, 1e-300);
        if (mode == SlopeMode::Minmod) slope[i] = minmodState(left, right);
        else slope[i] = (U[ip] - U[im]) / std::max(xip - xim, 1e-300);
    }
    return slope;
}

State1D EulerALE1D::boundaryFluxLeft(const State1D& U0, const Grid1D& grid,
                                     const BoundaryCondition1D& bc) const {
    double w = grid.wFace.front();
    switch (bc.left) {
        case BoundaryKind::FixedWall:
            return movingWallFlux(U0, 0.0);
        case BoundaryKind::MovingWall:
            return movingWallFlux(U0, w);
        case BoundaryKind::Dirichlet:
            return aleRusanovFlux(primToCons(bc.leftPrimitive), U0, w);
        case BoundaryKind::Transmissive:
        default:
            return aleRusanovFlux(U0, U0, w);
    }
}

State1D EulerALE1D::boundaryFluxRight(const State1D& UN, const Grid1D& grid,
                                      const BoundaryCondition1D& bc) const {
    double w = grid.wFace.back();
    switch (bc.right) {
        case BoundaryKind::FixedWall:
            return movingWallFlux(UN, 0.0);
        case BoundaryKind::MovingWall:
            return movingWallFlux(UN, w);
        case BoundaryKind::Dirichlet:
            return aleRusanovFlux(UN, primToCons(bc.rightPrimitive), w);
        case BoundaryKind::Transmissive:
        default:
            return aleRusanovFlux(UN, UN, w);
    }
}

std::vector<State1D> EulerALE1D::rhs(const std::vector<State1D>& q, const Grid1D& grid,
                                    const BoundaryCondition1D& bc, SlopeMode slopes) const {
    std::vector<State1D> U(n_);
    for (int i = 0; i < n_; ++i) U[i] = q[i] / grid.dx(i);
    std::vector<State1D> s = reconstructSlopes(U, grid, bc, slopes);

    std::vector<State1D> faceFlux(n_ + 1);
    bool periodic = (bc.left == BoundaryKind::Periodic && bc.right == BoundaryKind::Periodic);
    for (int f = 1; f < n_; ++f) {
        int il = f - 1, ir = f;
        State1D UL = U[il] + s[il] * (grid.xFace[f] - grid.center(il));
        State1D UR = U[ir] + s[ir] * (grid.xFace[f] - grid.center(ir));
        faceFlux[f] = aleRusanovFlux(UL, UR, grid.wFace[f]);
    }
    if (periodic) {
        State1D UL = U[n_ - 1] + s[n_ - 1] * (grid.xFace.back() - grid.center(n_ - 1));
        State1D UR = U[0] + s[0] * (grid.xFace.front() - grid.center(0));
        faceFlux[0] = aleRusanovFlux(UL, UR, grid.wFace.front());
        faceFlux[n_] = faceFlux[0];
    } else {
        State1D Uleft = U[0] + s[0] * (grid.xFace.front() - grid.center(0));
        State1D Uright = U[n_ - 1] + s[n_ - 1] * (grid.xFace.back() - grid.center(n_ - 1));
        faceFlux[0] = boundaryFluxLeft(Uleft, grid, bc);
        faceFlux[n_] = boundaryFluxRight(Uright, grid, bc);
    }

    std::vector<State1D> r(n_);
    for (int i = 0; i < n_; ++i) r[i] = faceFlux[i] - faceFlux[i + 1];
    return r;
}

bool EulerALE1D::step(double dt, const GridFn& gridFn, const BoundaryFn& bcFn,
                      SlopeMode slopes) {
    Grid1D g0 = gridFn(time_);
    Grid1D g1 = gridFn(time_ + dt);
    Grid1D g0Stage = withSweptFaceVelocity(g0, g1, dt);
    Grid1D g1Stage = g1;
    g1Stage.wFace = g0Stage.wFace;
    auto r0 = rhs(q_, g0Stage, bcFn(time_), slopes);
    std::vector<State1D> qPred(n_);
    for (int i = 0; i < n_; ++i) qPred[i] = q_[i] + dt * r0[i];

    auto r1 = rhs(qPred, g1Stage, bcFn(time_ + dt), slopes);
    std::vector<State1D> qNew(n_);
    for (int i = 0; i < n_; ++i) qNew[i] = 0.5 * q_[i] + 0.5 * (qPred[i] + dt * r1[i]);

    q_ = qNew;
    grid_ = g1;
    time_ += dt;
    bool ok = minDensity() > 0.0 && minPressure() > 0.0;
    if (!ok && std::getenv("EULER_FSI_DEBUG")) {
        std::cerr << "EulerALE1D non-positive state: t=" << time_
                  << " dt=" << dt << " minRho=" << minDensity()
                  << " minP=" << minPressure() << "\n";
    }
    return ok;
}

std::array<double, 3> EulerALE1D::conservedTotals() const {
    std::array<double, 3> total{0.0, 0.0, 0.0};
    for (const State1D& q : q_) {
        total[0] += q.rho;
        total[1] += q.mom;
        total[2] += q.ene;
    }
    return total;
}

double EulerALE1D::maxWaveSpeed() const {
    double m = 0.0;
    for (int i = 0; i < n_; ++i) {
        State1D U = cellAverage(i);
        double u = U.mom / std::max(U.rho, RHO_FLOOR);
        m = std::max(m, std::abs(u) + soundSpeed(U));
    }
    return m;
}

double EulerALE1D::minDensity() const {
    double m = 1e300;
    for (int i = 0; i < n_; ++i) m = std::min(m, cellAverage(i).rho);
    return m;
}

double EulerALE1D::minPressure() const {
    double m = 1e300;
    for (int i = 0; i < n_; ++i) m = std::min(m, pressure(cellAverage(i)));
    return m;
}

double SpringPiston::acceleration(double gasPressure) const {
    double force = area * (gasPressure - externalPressure)
                 - stiffness * (x - xRest) - damping * v;
    return force / std::max(mass, 1e-300);
}

double SpringPiston::kineticEnergy() const {
    return 0.5 * mass * v * v;
}

double SpringPiston::springEnergy() const {
    double dx = x - xRest;
    return 0.5 * stiffness * dx * dx;
}

double SpringPiston::externalPressurePotential() const {
    return area * externalPressure * x;
}

void SpringPiston::advanceExplicit(double dt, double gasPressure) {
    double a = acceleration(gasPressure);
    v += dt * a;
    x += dt * v;
    if (x < minX) { x = minX; v = std::max(0.0, v); }
    if (x > maxX) { x = maxX; v = std::min(0.0, v); }
}

double densityL2Error(const EulerALE1D& solver,
                      const std::function<double(double, double)>& exactRho,
                      double t) {
    const Grid1D& g = solver.grid();
    double e = 0.0;
    const double gp[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    for (int i = 0; i < solver.cells(); ++i) {
        double xa = g.xFace[i], xb = g.xFace[i + 1], dx = xb - xa;
        double avg = 0.0;
        for (double s : gp) {
            double x = 0.5 * (xa + xb) + 0.5 * dx * s;
            avg += 0.5 * exactRho(x, t);
        }
        double d = solver.cellAverage(i).rho - avg;
        e += dx * d * d;
    }
    return std::sqrt(e);
}

void writePistonFramePPM(const std::string& path, const EulerALE1D& solver,
                         const SpringPiston& piston, int width, int height,
                         double rhoMin, double rhoMax) {
    std::vector<unsigned char> img(static_cast<size_t>(width) * height * 3, 245);
    const Grid1D& g = solver.grid();
    double xmax = std::max(piston.maxX, g.xFace.back());
    int bandTop = height / 5, bandBot = 4 * height / 5;
    for (int col = 0; col < width; ++col) {
        double x = xmax * (col + 0.5) / width;
        int cell = -1;
        if (x >= g.xFace.front() && x <= g.xFace.back()) {
            auto it = std::upper_bound(g.xFace.begin(), g.xFace.end(), x);
            cell = std::max(0, std::min(solver.cells() - 1, static_cast<int>(it - g.xFace.begin()) - 1));
        }
        unsigned char r = 230, gg = 230, b = 235;
        if (cell >= 0) {
            double rho = solver.cellAverage(cell).rho;
            colour((rho - rhoMin) / std::max(rhoMax - rhoMin, 1e-12), r, gg, b);
        }
        for (int row = bandTop; row <= bandBot; ++row) {
            size_t k = (static_cast<size_t>(row) * width + col) * 3;
            img[k] = r; img[k + 1] = gg; img[k + 2] = b;
        }
    }
    auto drawV = [&](double x, unsigned char r, unsigned char gg, unsigned char b, int thick) {
        int c = static_cast<int>(std::lround(x / xmax * (width - 1)));
        for (int dc = -thick; dc <= thick; ++dc) {
            int col = c + dc;
            if (col < 0 || col >= width) continue;
            for (int row = bandTop - 8; row <= bandBot + 8; ++row) {
                if (row < 0 || row >= height) continue;
                size_t k = (static_cast<size_t>(row) * width + col) * 3;
                img[k] = r; img[k + 1] = gg; img[k + 2] = b;
            }
        }
    };
    drawV(g.xFace.front(), 30, 30, 30, 1);
    drawV(g.xFace.back(), 235, 235, 245, 4);
    drawV(g.xFace.back(), 20, 20, 25, 2);

    // A compact sparkline of piston displacement in the upper margin.
    int xpix = static_cast<int>(std::lround(g.xFace.back() / xmax * (width - 1)));
    int ypix = height / 10;
    for (int row = std::max(0, ypix - 3); row <= std::min(height - 1, ypix + 3); ++row)
        for (int col = std::max(0, xpix - 3); col <= std::min(width - 1, xpix + 3); ++col) {
            size_t k = (static_cast<size_t>(row) * width + col) * 3;
            img[k] = 220; img[k + 1] = 40; img[k + 2] = 40;
        }
    writePPM(path, width, height, img);
}

void writeDensityStillPPM(const std::string& path, const EulerALE1D& solver,
                          int width, int height, double rhoMin, double rhoMax) {
    SpringPiston p;
    p.x = solver.grid().xFace.back();
    p.maxX = solver.grid().xFace.back();
    writePistonFramePPM(path, solver, p, width, height, rhoMin, rhoMax);
}

} // namespace euler_fsi
