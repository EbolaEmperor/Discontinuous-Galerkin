#include "EulerALE2D.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace euler_fsi {

namespace {

constexpr double RHO_FLOOR = 1e-12;
constexpr double P_FLOOR = 1e-12;
constexpr double PI = 3.141592653589793238462643383279502884;

double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

State2D minmodState(const State2D& a, const State2D& b) {
    return {minmod(a.rho, b.rho), minmod(a.mx, b.mx),
            minmod(a.my, b.my), minmod(a.ene, b.ene)};
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
        std::cerr << "writePPM2D: cannot open " << path << "\n";
        return;
    }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
}

Grid2D withSweptVelocity(const Grid2D& g0, const Grid2D& g1, double dt) {
    Grid2D g = g0;
    if (dt <= 0.0) return g;
    if (g0.xFace.size() == g1.xFace.size()) {
        g.wXFace.resize(g.xFace.size());
        for (size_t i = 0; i < g.xFace.size(); ++i) g.wXFace[i] = (g1.xFace[i] - g0.xFace[i]) / dt;
    }
    if (g0.yFace.size() == g1.yFace.size()) {
        g.wYFace.resize(g.yFace.size());
        for (size_t j = 0; j < g.yFace.size(); ++j) g.wYFace[j] = (g1.yFace[j] - g0.yFace[j]) / dt;
    }
    return g;
}

} // namespace

State2D operator+(const State2D& a, const State2D& b) {
    return {a.rho + b.rho, a.mx + b.mx, a.my + b.my, a.ene + b.ene};
}
State2D operator-(const State2D& a, const State2D& b) {
    return {a.rho - b.rho, a.mx - b.mx, a.my - b.my, a.ene - b.ene};
}
State2D operator*(double s, const State2D& a) {
    return {s * a.rho, s * a.mx, s * a.my, s * a.ene};
}
State2D operator*(const State2D& a, double s) { return s * a; }
State2D operator/(const State2D& a, double s) {
    return {a.rho / s, a.mx / s, a.my / s, a.ene / s};
}

double pressure(const State2D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    double ke = 0.5 * (U.mx * U.mx + U.my * U.my) / rho;
    return (GAMMA - 1.0) * (U.ene - ke);
}

double soundSpeed(const State2D& U) {
    return std::sqrt(GAMMA * std::max(pressure(U), P_FLOOR) / std::max(U.rho, RHO_FLOOR));
}

Prim2D consToPrim(const State2D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    return {rho, U.mx / rho, U.my / rho, std::max(pressure(U), P_FLOOR)};
}

State2D primToCons(const Prim2D& q) {
    double rho = std::max(q.rho, RHO_FLOOR);
    double p = std::max(q.p, P_FLOOR);
    return {rho, rho * q.u, rho * q.v,
            p / (GAMMA - 1.0) + 0.5 * rho * (q.u * q.u + q.v * q.v)};
}

State2D fluxX(const State2D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    double u = U.mx / rho;
    double v = U.my / rho;
    double p = pressure(U);
    return {U.mx, U.mx * u + p, U.mx * v, (U.ene + p) * u};
}

State2D fluxY(const State2D& U) {
    double rho = std::max(U.rho, RHO_FLOOR);
    double u = U.mx / rho;
    double v = U.my / rho;
    double p = pressure(U);
    return {U.my, U.my * u, U.my * v + p, (U.ene + p) * v};
}

State2D aleRusanovFluxX(const State2D& UL, const State2D& UR, double wFace) {
    State2D FL = fluxX(UL) - wFace * UL;
    State2D FR = fluxX(UR) - wFace * UR;
    double uL = UL.mx / std::max(UL.rho, RHO_FLOOR);
    double uR = UR.mx / std::max(UR.rho, RHO_FLOOR);
    double lam = std::max(std::abs(uL - wFace) + soundSpeed(UL),
                          std::abs(uR - wFace) + soundSpeed(UR));
    return 0.5 * (FL + FR) - 0.5 * lam * (UR - UL);
}

State2D aleRusanovFluxY(const State2D& UB, const State2D& UT, double wFace) {
    State2D FB = fluxY(UB) - wFace * UB;
    State2D FT = fluxY(UT) - wFace * UT;
    double vB = UB.my / std::max(UB.rho, RHO_FLOOR);
    double vT = UT.my / std::max(UT.rho, RHO_FLOOR);
    double lam = std::max(std::abs(vB - wFace) + soundSpeed(UB),
                          std::abs(vT - wFace) + soundSpeed(UT));
    return 0.5 * (FB + FT) - 0.5 * lam * (UT - UB);
}

State2D movingWallFluxX(const State2D& Uinside, double wallVelocityX) {
    double p = std::max(pressure(Uinside), P_FLOOR);
    return {0.0, p, 0.0, p * wallVelocityX};
}

State2D movingWallFluxY(const State2D& Uinside, double wallVelocityY) {
    double p = std::max(pressure(Uinside), P_FLOOR);
    return {0.0, 0.0, p, p * wallVelocityY};
}

Grid2D Grid2D::affinePiston(int nxIn, int nyIn, double lengthX, double lengthY,
                            double velocityX) {
    Grid2D g;
    g.nx = nxIn;
    g.ny = nyIn;
    g.xFace.resize(nxIn + 1);
    g.yFace.resize(nyIn + 1);
    g.wXFace.resize(nxIn + 1);
    g.wYFace.assign(nyIn + 1, 0.0);
    for (int i = 0; i <= nxIn; ++i) {
        double xi = static_cast<double>(i) / nxIn;
        g.xFace[i] = lengthX * xi;
        g.wXFace[i] = velocityX * xi;
    }
    for (int j = 0; j <= nyIn; ++j) {
        double eta = static_cast<double>(j) / nyIn;
        g.yFace[j] = lengthY * eta;
    }
    return g;
}

Grid2D Grid2D::oscillatoryX(int nxIn, int nyIn, double lengthX, double lengthY,
                            double t, double amplitude, double omega) {
    Grid2D g;
    g.nx = nxIn;
    g.ny = nyIn;
    g.xFace.resize(nxIn + 1);
    g.yFace.resize(nyIn + 1);
    g.wXFace.resize(nxIn + 1);
    g.wYFace.assign(nyIn + 1, 0.0);
    for (int i = 0; i <= nxIn; ++i) {
        double eta = static_cast<double>(i) / nxIn;
        double shape = std::sin(2.0 * PI * eta);
        g.xFace[i] = lengthX * (eta + amplitude * shape * std::sin(omega * t));
        g.wXFace[i] = lengthX * amplitude * shape * omega * std::cos(omega * t);
    }
    for (int j = 0; j <= nyIn; ++j) {
        double eta = static_cast<double>(j) / nyIn;
        g.yFace[j] = lengthY * eta;
    }
    return g;
}

EulerALE2D::EulerALE2D(int nxIn, int nyIn)
    : nx_(nxIn), ny_(nyIn), q_(nxIn * nyIn) {}

void EulerALE2D::setGrid(const Grid2D& grid, double time) {
    grid_ = grid;
    time_ = time;
}

void EulerALE2D::setCellAverages(const std::function<Prim2D(double, double)>& prim,
                                 const Grid2D& grid) {
    grid_ = grid;
    q_.assign(cells(), {});
    const double gp[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) {
            double xa = grid.xFace[i], xb = grid.xFace[i + 1];
            double ya = grid.yFace[j], yb = grid.yFace[j + 1];
            double dx = xb - xa, dy = yb - ya;
            State2D integral{};
            for (double sx : gp)
                for (double sy : gp) {
                    double x = 0.5 * (xa + xb) + 0.5 * dx * sx;
                    double y = 0.5 * (ya + yb) + 0.5 * dy * sy;
                    integral = integral + (0.25 * dx * dy) * primToCons(prim(x, y));
                }
            q_[index(i, j)] = integral;
        }
}

void EulerALE2D::setUniform(const Prim2D& prim, const Grid2D& grid) {
    setCellAverages([&](double, double) { return prim; }, grid);
}

State2D EulerALE2D::cellAverage(int i, int j) const {
    return q_[index(i, j)] / grid_.area(i, j);
}

void EulerALE2D::reconstructSlopes(const std::vector<State2D>& U, const Grid2D& grid,
                                   const BoundaryCondition2D& bc, SlopeMode mode,
                                   std::vector<State2D>& sx, std::vector<State2D>& sy) const {
    sx.assign(cells(), {});
    sy.assign(cells(), {});
    if (mode == SlopeMode::FirstOrder) return;
    bool perX = (bc.left == BoundaryKind::Periodic && bc.right == BoundaryKind::Periodic);
    bool perY = (bc.bottom == BoundaryKind::Periodic && bc.top == BoundaryKind::Periodic);
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) {
            int id = index(i, j);
            if (perX || (i > 0 && i + 1 < nx_)) {
                int im = perX ? (i + nx_ - 1) % nx_ : i - 1;
                int ip = perX ? (i + 1) % nx_ : i + 1;
                double xm = grid.centerX(im), x = grid.centerX(i), xp = grid.centerX(ip);
                if (perX) {
                    double L = grid.lengthX();
                    if (im > i) xm -= L;
                    if (ip < i) xp += L;
                }
                State2D l = (U[id] - U[index(im, j)]) / std::max(x - xm, 1e-300);
                State2D r = (U[index(ip, j)] - U[id]) / std::max(xp - x, 1e-300);
                sx[id] = (mode == SlopeMode::Minmod) ? minmodState(l, r)
                                                     : (U[index(ip, j)] - U[index(im, j)]) / std::max(xp - xm, 1e-300);
            }
            if (perY || (j > 0 && j + 1 < ny_)) {
                int jm = perY ? (j + ny_ - 1) % ny_ : j - 1;
                int jp = perY ? (j + 1) % ny_ : j + 1;
                double ym = grid.centerY(jm), y = grid.centerY(j), yp = grid.centerY(jp);
                if (perY) {
                    double L = grid.lengthY();
                    if (jm > j) ym -= L;
                    if (jp < j) yp += L;
                }
                State2D l = (U[id] - U[index(i, jm)]) / std::max(y - ym, 1e-300);
                State2D r = (U[index(i, jp)] - U[id]) / std::max(yp - y, 1e-300);
                sy[id] = (mode == SlopeMode::Minmod) ? minmodState(l, r)
                                                     : (U[index(i, jp)] - U[index(i, jm)]) / std::max(yp - ym, 1e-300);
            }
        }
}

State2D EulerALE2D::boundaryFluxLeft(const State2D& U, const Grid2D& grid,
                                     const BoundaryCondition2D& bc) const {
    double w = grid.wXFace.front();
    switch (bc.left) {
        case BoundaryKind::FixedWall: return movingWallFluxX(U, 0.0);
        case BoundaryKind::MovingWall: return movingWallFluxX(U, w);
        case BoundaryKind::Dirichlet: return aleRusanovFluxX(primToCons(bc.leftPrimitive), U, w);
        case BoundaryKind::Transmissive:
        default: return aleRusanovFluxX(U, U, w);
    }
}

State2D EulerALE2D::boundaryFluxRight(const State2D& U, const Grid2D& grid,
                                      const BoundaryCondition2D& bc) const {
    double w = grid.wXFace.back();
    switch (bc.right) {
        case BoundaryKind::FixedWall: return movingWallFluxX(U, 0.0);
        case BoundaryKind::MovingWall: return movingWallFluxX(U, w);
        case BoundaryKind::Dirichlet: return aleRusanovFluxX(U, primToCons(bc.rightPrimitive), w);
        case BoundaryKind::Transmissive:
        default: return aleRusanovFluxX(U, U, w);
    }
}

State2D EulerALE2D::boundaryFluxBottom(const State2D& U, const Grid2D& grid,
                                       const BoundaryCondition2D& bc) const {
    double w = grid.wYFace.front();
    switch (bc.bottom) {
        case BoundaryKind::FixedWall: return movingWallFluxY(U, 0.0);
        case BoundaryKind::MovingWall: return movingWallFluxY(U, w);
        case BoundaryKind::Dirichlet: return aleRusanovFluxY(primToCons(bc.bottomPrimitive), U, w);
        case BoundaryKind::Transmissive:
        default: return aleRusanovFluxY(U, U, w);
    }
}

State2D EulerALE2D::boundaryFluxTop(const State2D& U, const Grid2D& grid,
                                    const BoundaryCondition2D& bc) const {
    double w = grid.wYFace.back();
    switch (bc.top) {
        case BoundaryKind::FixedWall: return movingWallFluxY(U, 0.0);
        case BoundaryKind::MovingWall: return movingWallFluxY(U, w);
        case BoundaryKind::Dirichlet: return aleRusanovFluxY(U, primToCons(bc.topPrimitive), w);
        case BoundaryKind::Transmissive:
        default: return aleRusanovFluxY(U, U, w);
    }
}

std::vector<State2D> EulerALE2D::rhs(const std::vector<State2D>& q, const Grid2D& grid,
                                    const BoundaryCondition2D& bc, SlopeMode slopes) const {
    std::vector<State2D> U(cells());
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) U[index(i, j)] = q[index(i, j)] / grid.area(i, j);

    std::vector<State2D> sx, sy;
    reconstructSlopes(U, grid, bc, slopes, sx, sy);

    auto fxId = [&](int f, int j) { return j * (nx_ + 1) + f; };
    auto fyId = [&](int i, int f) { return f * nx_ + i; };
    std::vector<State2D> fx((nx_ + 1) * ny_);
    std::vector<State2D> fy(nx_ * (ny_ + 1));
    bool perX = (bc.left == BoundaryKind::Periodic && bc.right == BoundaryKind::Periodic);
    bool perY = (bc.bottom == BoundaryKind::Periodic && bc.top == BoundaryKind::Periodic);

    for (int j = 0; j < ny_; ++j) {
        for (int f = 1; f < nx_; ++f) {
            int il = f - 1, ir = f;
            int l = index(il, j), r = index(ir, j);
            State2D UL = U[l] + sx[l] * (grid.xFace[f] - grid.centerX(il));
            State2D UR = U[r] + sx[r] * (grid.xFace[f] - grid.centerX(ir));
            fx[fxId(f, j)] = aleRusanovFluxX(UL, UR, grid.wXFace[f]);
        }
        if (perX) {
            int l = index(nx_ - 1, j), r = index(0, j);
            State2D UL = U[l] + sx[l] * (grid.xFace.back() - grid.centerX(nx_ - 1));
            State2D UR = U[r] + sx[r] * (grid.xFace.front() - grid.centerX(0));
            fx[fxId(0, j)] = aleRusanovFluxX(UL, UR, grid.wXFace.front());
            fx[fxId(nx_, j)] = fx[fxId(0, j)];
        } else {
            int l = index(0, j), r = index(nx_ - 1, j);
            State2D UL = U[l] + sx[l] * (grid.xFace.front() - grid.centerX(0));
            State2D UR = U[r] + sx[r] * (grid.xFace.back() - grid.centerX(nx_ - 1));
            fx[fxId(0, j)] = boundaryFluxLeft(UL, grid, bc);
            fx[fxId(nx_, j)] = boundaryFluxRight(UR, grid, bc);
        }
    }

    for (int i = 0; i < nx_; ++i) {
        for (int f = 1; f < ny_; ++f) {
            int jb = f - 1, jt = f;
            int b = index(i, jb), t = index(i, jt);
            State2D UB = U[b] + sy[b] * (grid.yFace[f] - grid.centerY(jb));
            State2D UT = U[t] + sy[t] * (grid.yFace[f] - grid.centerY(jt));
            fy[fyId(i, f)] = aleRusanovFluxY(UB, UT, grid.wYFace[f]);
        }
        if (perY) {
            int b = index(i, ny_ - 1), t = index(i, 0);
            State2D UB = U[b] + sy[b] * (grid.yFace.back() - grid.centerY(ny_ - 1));
            State2D UT = U[t] + sy[t] * (grid.yFace.front() - grid.centerY(0));
            fy[fyId(i, 0)] = aleRusanovFluxY(UB, UT, grid.wYFace.front());
            fy[fyId(i, ny_)] = fy[fyId(i, 0)];
        } else {
            int b = index(i, 0), t = index(i, ny_ - 1);
            State2D UB = U[b] + sy[b] * (grid.yFace.front() - grid.centerY(0));
            State2D UT = U[t] + sy[t] * (grid.yFace.back() - grid.centerY(ny_ - 1));
            fy[fyId(i, 0)] = boundaryFluxBottom(UB, grid, bc);
            fy[fyId(i, ny_)] = boundaryFluxTop(UT, grid, bc);
        }
    }

    std::vector<State2D> r(cells());
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) {
            r[index(i, j)] =
                grid.dy(j) * (fx[fxId(i, j)] - fx[fxId(i + 1, j)]) +
                grid.dx(i) * (fy[fyId(i, j)] - fy[fyId(i, j + 1)]);
        }
    return r;
}

bool EulerALE2D::step(double dt, const Grid2DFn& gridFn, const Boundary2DFn& bcFn,
                      SlopeMode slopes) {
    Grid2D g0 = gridFn(time_);
    Grid2D g1 = gridFn(time_ + dt);
    Grid2D g0Stage = withSweptVelocity(g0, g1, dt);
    Grid2D g1Stage = g1;
    g1Stage.wXFace = g0Stage.wXFace;
    g1Stage.wYFace = g0Stage.wYFace;

    auto r0 = rhs(q_, g0Stage, bcFn(time_), slopes);
    std::vector<State2D> qPred(cells());
    for (int k = 0; k < cells(); ++k) qPred[k] = q_[k] + dt * r0[k];

    auto r1 = rhs(qPred, g1Stage, bcFn(time_ + dt), slopes);
    std::vector<State2D> qNew(cells());
    for (int k = 0; k < cells(); ++k) qNew[k] = 0.5 * q_[k] + 0.5 * (qPred[k] + dt * r1[k]);

    q_ = qNew;
    grid_ = g1;
    time_ += dt;
    bool ok = minDensity() > 0.0 && minPressure() > 0.0;
    if (!ok && std::getenv("EULER_FSI_DEBUG")) {
        std::cerr << "EulerALE2D non-positive state: t=" << time_
                  << " dt=" << dt << " minRho=" << minDensity()
                  << " minP=" << minPressure() << "\n";
    }
    return ok;
}

std::array<double, 4> EulerALE2D::conservedTotals() const {
    std::array<double, 4> total{0.0, 0.0, 0.0, 0.0};
    for (const State2D& q : q_) {
        total[0] += q.rho;
        total[1] += q.mx;
        total[2] += q.my;
        total[3] += q.ene;
    }
    return total;
}

double EulerALE2D::maxWaveSpeed() const {
    double m = 0.0;
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) {
            State2D U = cellAverage(i, j);
            double rho = std::max(U.rho, RHO_FLOOR);
            double u = U.mx / rho, v = U.my / rho;
            m = std::max(m, std::sqrt(u * u + v * v) + soundSpeed(U));
        }
    return m;
}

double EulerALE2D::minDensity() const {
    double m = 1e300;
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) m = std::min(m, cellAverage(i, j).rho);
    return m;
}

double EulerALE2D::minPressure() const {
    double m = 1e300;
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) m = std::min(m, pressure(cellAverage(i, j)));
    return m;
}

double densityL2Error(const EulerALE2D& solver,
                      const std::function<double(double, double, double)>& exactRho,
                      double t) {
    const Grid2D& g = solver.grid();
    const double gp[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    double e = 0.0;
    for (int j = 0; j < solver.ny(); ++j)
        for (int i = 0; i < solver.nx(); ++i) {
            double xa = g.xFace[i], xb = g.xFace[i + 1], dx = xb - xa;
            double ya = g.yFace[j], yb = g.yFace[j + 1], dy = yb - ya;
            double avg = 0.0;
            for (double sx : gp)
                for (double sy : gp) {
                    double x = 0.5 * (xa + xb) + 0.5 * dx * sx;
                    double y = 0.5 * (ya + yb) + 0.5 * dy * sy;
                    avg += 0.25 * exactRho(x, y, t);
                }
            double d = solver.cellAverage(i, j).rho - avg;
            e += dx * dy * d * d;
        }
    return std::sqrt(e);
}

void writeDensityPPM2D(const std::string& path, const EulerALE2D& solver,
                       int width, int height, double rhoMin, double rhoMax) {
    std::vector<unsigned char> img(static_cast<size_t>(width) * height * 3, 245);
    const Grid2D& g = solver.grid();
    for (int row = 0; row < height; ++row) {
        double y = g.yFace.front() + (g.yFace.back() - g.yFace.front()) *
                                  (height - row - 0.5) / height;
        auto jyIt = std::upper_bound(g.yFace.begin(), g.yFace.end(), y);
        int j = std::max(0, std::min(solver.ny() - 1, static_cast<int>(jyIt - g.yFace.begin()) - 1));
        for (int col = 0; col < width; ++col) {
            double x = g.xFace.front() + (g.xFace.back() - g.xFace.front()) *
                                      (col + 0.5) / width;
            auto ixIt = std::upper_bound(g.xFace.begin(), g.xFace.end(), x);
            int i = std::max(0, std::min(solver.nx() - 1, static_cast<int>(ixIt - g.xFace.begin()) - 1));
            double rho = solver.cellAverage(i, j).rho;
            unsigned char r, gg, b;
            colour((rho - rhoMin) / std::max(rhoMax - rhoMin, 1e-12), r, gg, b);
            size_t k = (static_cast<size_t>(row) * width + col) * 3;
            img[k] = r; img[k + 1] = gg; img[k + 2] = b;
        }
    }
    writePPM(path, width, height, img);
}

void writePistonFramePPM2D(const std::string& path, const EulerALE2D& solver,
                           const SpringPiston& piston, int width, int height,
                           double rhoMin, double rhoMax) {
    std::vector<unsigned char> img(static_cast<size_t>(width) * height * 3, 238);
    const Grid2D& g = solver.grid();
    double xmin = g.xFace.front();
    double xmax = std::max(piston.maxX, g.xFace.back());
    double ymin = g.yFace.front();
    double ymax = g.yFace.back();
    double xspan = std::max(xmax - xmin, 1e-300);
    double yspan = std::max(ymax - ymin, 1e-300);

    auto setPixel = [&](int col, int row, unsigned char r, unsigned char gg, unsigned char b) {
        if (col < 0 || col >= width || row < 0 || row >= height) return;
        size_t k = (static_cast<size_t>(row) * width + col) * 3;
        img[k] = r;
        img[k + 1] = gg;
        img[k + 2] = b;
    };

    auto blendPixel = [&](int col, int row, unsigned char r, unsigned char gg,
                          unsigned char b, double alpha) {
        if (col < 0 || col >= width || row < 0 || row >= height) return;
        alpha = std::max(0.0, std::min(1.0, alpha));
        if (alpha <= 0.0) return;
        size_t k = (static_cast<size_t>(row) * width + col) * 3;
        double beta = 1.0 - alpha;
        img[k] = static_cast<unsigned char>(std::lround(beta * img[k] + alpha * r));
        img[k + 1] = static_cast<unsigned char>(std::lround(beta * img[k + 1] + alpha * gg));
        img[k + 2] = static_cast<unsigned char>(std::lround(beta * img[k + 2] + alpha * b));
    };

    auto coverage1D = [](double pixelCenter, double a, double b) {
        if (a > b) std::swap(a, b);
        double lo = std::max(pixelCenter - 0.5, a);
        double hi = std::min(pixelCenter + 0.5, b);
        return std::max(0.0, hi - lo);
    };

    auto xToPx = [&](double x) {
        return (x - xmin) / xspan * (width - 1);
    };
    auto yToPx = [&](double y) {
        return (ymax - y) / yspan * (height - 1);
    };

    auto drawAARect = [&](double x0, double y0, double x1, double y1,
                          unsigned char r, unsigned char gg, unsigned char b,
                          double opacity = 1.0) {
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);
        int c0 = std::max(0, static_cast<int>(std::floor(x0 - 1.0)));
        int c1 = std::min(width - 1, static_cast<int>(std::ceil(x1 + 1.0)));
        int r0 = std::max(0, static_cast<int>(std::floor(y0 - 1.0)));
        int r1 = std::min(height - 1, static_cast<int>(std::ceil(y1 + 1.0)));
        for (int row = r0; row <= r1; ++row) {
            double ay = coverage1D(row, y0, y1);
            if (ay <= 0.0) continue;
            for (int col = c0; col <= c1; ++col) {
                double ax = coverage1D(col, x0, x1);
                blendPixel(col, row, r, gg, b, opacity * ax * ay);
            }
        }
    };

    auto drawAALine = [&](double x0, double y0, double x1, double y1,
                          unsigned char r, unsigned char gg, unsigned char b,
                          double halfWidth = 0.9, double opacity = 1.0) {
        double xminPx = std::min(x0, x1) - halfWidth - 2.0;
        double xmaxPx = std::max(x0, x1) + halfWidth + 2.0;
        double yminPx = std::min(y0, y1) - halfWidth - 2.0;
        double ymaxPx = std::max(y0, y1) + halfWidth + 2.0;
        int c0 = std::max(0, static_cast<int>(std::floor(xminPx)));
        int c1 = std::min(width - 1, static_cast<int>(std::ceil(xmaxPx)));
        int r0 = std::max(0, static_cast<int>(std::floor(yminPx)));
        int r1 = std::min(height - 1, static_cast<int>(std::ceil(ymaxPx)));
        double dx = x1 - x0, dy = y1 - y0;
        double len2 = dx * dx + dy * dy;
        for (int row = r0; row <= r1; ++row)
            for (int col = c0; col <= c1; ++col) {
                double a = 0.0;
                if (len2 > 1e-300) {
                    a = ((col - x0) * dx + (row - y0) * dy) / len2;
                    a = std::max(0.0, std::min(1.0, a));
                }
                double px = x0 + a * dx;
                double py = y0 + a * dy;
                double dist = std::hypot(col - px, row - py);
                double alpha = std::max(0.0, std::min(1.0, halfWidth + 0.5 - dist));
                blendPixel(col, row, r, gg, b, opacity * alpha);
            }
    };

    auto drawAACircle = [&](double cx, double cy, double radius,
                            unsigned char r, unsigned char gg, unsigned char b,
                            double opacity = 1.0) {
        int c0 = std::max(0, static_cast<int>(std::floor(cx - radius - 1.0)));
        int c1 = std::min(width - 1, static_cast<int>(std::ceil(cx + radius + 1.0)));
        int r0 = std::max(0, static_cast<int>(std::floor(cy - radius - 1.0)));
        int r1 = std::min(height - 1, static_cast<int>(std::ceil(cy + radius + 1.0)));
        for (int row = r0; row <= r1; ++row)
            for (int col = c0; col <= c1; ++col) {
                double dist = std::hypot(col - cx, row - cy);
                double alpha = std::max(0.0, std::min(1.0, radius + 0.5 - dist));
                blendPixel(col, row, r, gg, b, opacity * alpha);
            }
    };

    for (int row = 0; row < height; ++row) {
        double y = ymax - yspan * (row + 0.5) / height;
        bool insideY = (y >= ymin && y <= ymax);
        int j = 0;
        if (insideY) {
            auto jyIt = std::upper_bound(g.yFace.begin(), g.yFace.end(), y);
            j = std::max(0, std::min(solver.ny() - 1, static_cast<int>(jyIt - g.yFace.begin()) - 1));
        }
        for (int col = 0; col < width; ++col) {
            double x = xmin + xspan * (col + 0.5) / width;
            unsigned char r = 232, gg = 232, b = 238;
            if (insideY && x >= g.xFace.front() && x <= g.xFace.back()) {
                auto ixIt = std::upper_bound(g.xFace.begin(), g.xFace.end(), x);
                int i = std::max(0, std::min(solver.nx() - 1, static_cast<int>(ixIt - g.xFace.begin()) - 1));
                double rho = solver.cellAverage(i, j).rho;
                colour((rho - rhoMin) / std::max(rhoMax - rhoMin, 1e-12), r, gg, b);
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

    drawAARect(left - 1.0, top, left + 1.0, bottom, 35, 35, 38);
    drawAARect(left, top - 1.0, wall, top + 1.0, 35, 35, 38);
    drawAARect(left, bottom - 1.0, wall, bottom + 1.0, 35, 35, 38);

    double slabHalf = std::max(4.0, width / 180.0);
    drawAARect(wall - slabHalf - 3.0, top - 4.0, wall + slabHalf + 3.0, bottom + 4.0,
               238, 238, 246);
    drawAARect(wall - slabHalf, top - 4.0, wall + slabHalf, bottom + 4.0,
               16, 16, 22);
    if (std::abs(rest - wall) > 2.0) {
        drawAARect(rest - 0.8, top, rest + 0.8, bottom, 150, 150, 158, 0.72);
    }

    double markerY = std::max(5.0, top - height / 12.0);
    drawAACircle(wall, markerY, 4.4, 220, 40, 40);

    double springY = std::min(height - 4.0, bottom + height / 10.0);
    double anchor = xToPx(xmax);
    double start = std::min(wall + slabHalf + 4.0, width - 1.0);
    double end = std::max(start, anchor - 6.0);
    double amp = std::max(3.0, height / 70.0);
    int turns = 8;
    double prevX = start, prevY = springY;
    for (int k = 1; k <= turns * 2; ++k) {
        double x = start + (end - start) * k / (turns * 2);
        double y = springY + ((k % 2) ? amp : -amp);
        drawAALine(prevX, prevY, x, y, 70, 70, 78, 0.8);
        prevX = x;
        prevY = y;
    }
    drawAARect(anchor - 2.0, springY - 18.0, anchor + 2.0, springY + 18.0, 70, 70, 78);

    writePPM(path, width, height, img);
}

} // namespace euler_fsi
