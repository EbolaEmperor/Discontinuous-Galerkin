#ifndef EULER_ALE_2D_H
#define EULER_ALE_2D_H

#include "EulerALE1D.h"

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace euler_fsi {

struct State2D {
    double rho = 0.0;
    double mx = 0.0;
    double my = 0.0;
    double ene = 0.0;
};

struct Prim2D {
    double rho = 1.0;
    double u = 0.0;
    double v = 0.0;
    double p = 1.0;
};

State2D operator+(const State2D& a, const State2D& b);
State2D operator-(const State2D& a, const State2D& b);
State2D operator*(double s, const State2D& a);
State2D operator*(const State2D& a, double s);
State2D operator/(const State2D& a, double s);

double pressure(const State2D& U);
double soundSpeed(const State2D& U);
Prim2D consToPrim(const State2D& U);
State2D primToCons(const Prim2D& q);
State2D fluxX(const State2D& U);
State2D fluxY(const State2D& U);
State2D aleRusanovFluxX(const State2D& UL, const State2D& UR, double wFace);
State2D aleRusanovFluxY(const State2D& UB, const State2D& UT, double wFace);
State2D movingWallFluxX(const State2D& Uinside, double wallVelocityX);
State2D movingWallFluxY(const State2D& Uinside, double wallVelocityY);

struct Grid2D {
    int nx = 0;
    int ny = 0;
    std::vector<double> xFace;
    std::vector<double> yFace;
    std::vector<double> wXFace;
    std::vector<double> wYFace;

    double dx(int i) const { return xFace[i + 1] - xFace[i]; }
    double dy(int j) const { return yFace[j + 1] - yFace[j]; }
    double area(int i, int j) const { return dx(i) * dy(j); }
    double centerX(int i) const { return 0.5 * (xFace[i] + xFace[i + 1]); }
    double centerY(int j) const { return 0.5 * (yFace[j] + yFace[j + 1]); }
    double lengthX() const { return xFace.empty() ? 0.0 : xFace.back() - xFace.front(); }
    double lengthY() const { return yFace.empty() ? 0.0 : yFace.back() - yFace.front(); }

    static Grid2D affinePiston(int nx, int ny, double lengthX, double lengthY,
                               double velocityX);
    static Grid2D oscillatoryX(int nx, int ny, double lengthX, double lengthY,
                               double t, double amplitude, double omega);
};

struct BoundaryCondition2D {
    BoundaryKind left = BoundaryKind::Transmissive;
    BoundaryKind right = BoundaryKind::Transmissive;
    BoundaryKind bottom = BoundaryKind::Transmissive;
    BoundaryKind top = BoundaryKind::Transmissive;
    Prim2D leftPrimitive;
    Prim2D rightPrimitive;
    Prim2D bottomPrimitive;
    Prim2D topPrimitive;
};

using Grid2DFn = std::function<Grid2D(double)>;
using Boundary2DFn = std::function<BoundaryCondition2D(double)>;

class EulerALE2D {
public:
    EulerALE2D(int nx, int ny);

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int cells() const { return nx_ * ny_; }
    double time() const { return time_; }
    const Grid2D& grid() const { return grid_; }
    int index(int i, int j) const { return j * nx_ + i; }

    void setGrid(const Grid2D& grid, double time);
    void setCellAverages(const std::function<Prim2D(double, double)>& prim,
                         const Grid2D& grid);
    void setUniform(const Prim2D& prim, const Grid2D& grid);

    State2D cellAverage(int i, int j) const;
    Prim2D cellPrimitive(int i, int j) const { return consToPrim(cellAverage(i, j)); }
    const std::vector<State2D>& conservedIntegrals() const { return q_; }

    bool step(double dt, const Grid2DFn& gridFn, const Boundary2DFn& bcFn,
              SlopeMode slopes = SlopeMode::Minmod);

    std::array<double, 4> conservedTotals() const;
    double maxWaveSpeed() const;
    double minDensity() const;
    double minPressure() const;

private:
    std::vector<State2D> rhs(const std::vector<State2D>& q, const Grid2D& grid,
                             const BoundaryCondition2D& bc, SlopeMode slopes) const;
    void reconstructSlopes(const std::vector<State2D>& U, const Grid2D& grid,
                           const BoundaryCondition2D& bc, SlopeMode mode,
                           std::vector<State2D>& sx, std::vector<State2D>& sy) const;
    State2D boundaryFluxLeft(const State2D& U, const Grid2D& grid,
                             const BoundaryCondition2D& bc) const;
    State2D boundaryFluxRight(const State2D& U, const Grid2D& grid,
                              const BoundaryCondition2D& bc) const;
    State2D boundaryFluxBottom(const State2D& U, const Grid2D& grid,
                               const BoundaryCondition2D& bc) const;
    State2D boundaryFluxTop(const State2D& U, const Grid2D& grid,
                            const BoundaryCondition2D& bc) const;

    int nx_ = 0;
    int ny_ = 0;
    double time_ = 0.0;
    Grid2D grid_;
    std::vector<State2D> q_;
};

double densityL2Error(const EulerALE2D& solver,
                      const std::function<double(double, double, double)>& exactRho,
                      double t);

void writeDensityPPM2D(const std::string& path, const EulerALE2D& solver,
                       int width, int height, double rhoMin, double rhoMax);

void writePistonFramePPM2D(const std::string& path, const EulerALE2D& solver,
                           const SpringPiston& piston, int width, int height,
                           double rhoMin, double rhoMax);

} // namespace euler_fsi

#endif
