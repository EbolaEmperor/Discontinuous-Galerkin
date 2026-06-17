#ifndef EULER_ALE_1D_H
#define EULER_ALE_1D_H

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace euler_fsi {

constexpr double GAMMA = 1.4;

struct State1D {
    double rho = 0.0;
    double mom = 0.0;
    double ene = 0.0;
};

struct Prim1D {
    double rho = 1.0;
    double u = 0.0;
    double p = 1.0;
};

State1D operator+(const State1D& a, const State1D& b);
State1D operator-(const State1D& a, const State1D& b);
State1D operator*(double s, const State1D& a);
State1D operator*(const State1D& a, double s);
State1D operator/(const State1D& a, double s);

double pressure(const State1D& U);
double soundSpeed(const State1D& U);
Prim1D consToPrim(const State1D& U);
State1D primToCons(const Prim1D& q);
State1D flux(const State1D& U);
State1D aleRusanovFlux(const State1D& UL, const State1D& UR, double wFace);
State1D movingWallFlux(const State1D& Uinside, double wallVelocity);

enum class BoundaryKind {
    Periodic,
    Transmissive,
    FixedWall,
    MovingWall,
    Dirichlet
};

struct BoundaryCondition1D {
    BoundaryKind left = BoundaryKind::Transmissive;
    BoundaryKind right = BoundaryKind::Transmissive;
    Prim1D leftPrimitive;
    Prim1D rightPrimitive;
    double leftWallVelocity = 0.0;
    double rightWallVelocity = 0.0;
};

enum class SlopeMode {
    Unlimited,
    Minmod,
    FirstOrder
};

struct Grid1D {
    std::vector<double> xFace;
    std::vector<double> wFace;

    int cells() const { return static_cast<int>(xFace.size()) - 1; }
    double dx(int i) const { return xFace[i + 1] - xFace[i]; }
    double center(int i) const { return 0.5 * (xFace[i] + xFace[i + 1]); }
    double length() const { return xFace.empty() ? 0.0 : xFace.back() - xFace.front(); }

    static Grid1D affinePiston(int n, double length, double velocity);
    static Grid1D oscillatory(int n, double length, double t, double amplitude, double omega);
};

using GridFn = std::function<Grid1D(double)>;
using BoundaryFn = std::function<BoundaryCondition1D(double)>;

class EulerALE1D {
public:
    explicit EulerALE1D(int cells);

    int cells() const { return n_; }
    double time() const { return time_; }
    const Grid1D& grid() const { return grid_; }

    void setGrid(const Grid1D& grid, double time);
    void setCellAverages(const std::function<Prim1D(double)>& prim, const Grid1D& grid);
    void setUniform(const Prim1D& prim, const Grid1D& grid);

    State1D cellAverage(int i) const;
    Prim1D cellPrimitive(int i) const { return consToPrim(cellAverage(i)); }
    const std::vector<State1D>& conservedIntegrals() const { return q_; }

    bool step(double dt, const GridFn& gridFn, const BoundaryFn& bcFn,
              SlopeMode slopes = SlopeMode::Minmod);

    std::array<double, 3> conservedTotals() const;
    double maxWaveSpeed() const;
    double minDensity() const;
    double minPressure() const;

private:
    std::vector<State1D> rhs(const std::vector<State1D>& q, const Grid1D& grid,
                             const BoundaryCondition1D& bc, SlopeMode slopes) const;
    std::vector<State1D> reconstructSlopes(const std::vector<State1D>& U,
                                           const Grid1D& grid,
                                           const BoundaryCondition1D& bc,
                                           SlopeMode mode) const;
    State1D boundaryFluxLeft(const State1D& U0, const Grid1D& grid,
                             const BoundaryCondition1D& bc) const;
    State1D boundaryFluxRight(const State1D& UN, const Grid1D& grid,
                              const BoundaryCondition1D& bc) const;

    int n_ = 0;
    double time_ = 0.0;
    Grid1D grid_;
    std::vector<State1D> q_;   // cell integrals over the current moving cells
};

struct SpringPiston {
    double x = 1.0;
    double v = 0.0;
    double mass = 1.0;
    double stiffness = 100.0;
    double damping = 1.0;
    double xRest = 1.0;
    double area = 1.0;
    double externalPressure = 1.0;
    double minX = 0.25;
    double maxX = 2.0;

    double acceleration(double gasPressure) const;
    double kineticEnergy() const;
    double springEnergy() const;
    double externalPressurePotential() const;
    void advanceExplicit(double dt, double gasPressure);
};

double densityL2Error(const EulerALE1D& solver,
                      const std::function<double(double, double)>& exactRho,
                      double t);

void writePistonFramePPM(const std::string& path, const EulerALE1D& solver,
                         const SpringPiston& piston, int width, int height,
                         double rhoMin, double rhoMax);

void writeDensityStillPPM(const std::string& path, const EulerALE1D& solver,
                          int width, int height, double rhoMin, double rhoMax);

} // namespace euler_fsi

#endif
