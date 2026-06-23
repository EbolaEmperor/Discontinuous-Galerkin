#ifndef EULER_DG_H
#define EULER_DG_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "FEM.h"
#include "Mesh.h"

// ===========================================================================
// Discontinuous-Galerkin (dP_k) solver for the 2-D compressible Euler equations
//
//        U_t + div F(U) = 0,     U = (rho, rho u, rho v, E),
//        F = ( F_x(U), F_y(U) ),  p = (gamma-1)(E - 1/2 rho|u|^2),
//
// with a 2nd-order IMEX (implicit-explicit) Runge-Kutta time discretisation:
//   * the inviscid flux divergence is treated EXPLICITLY (non-stiff hyperbolic
//     part) with a Rusanov / local Lax-Friedrichs numerical flux, and
//   * a Persson-Peraire sub-cell ARTIFICIAL-VISCOSITY diffusion -div(eps grad U)
//     -- the shock-capturing term -- is treated IMPLICITLY (stiff parabolic
//     part), so the stiff diffusion never restricts the time step.
// The two are combined by the L-stable, 2nd-order ARS(2,2,2) IMEX-RK scheme
// (Ascher-Ruuth-Spiteri).  When the artificial viscosity vanishes (smooth flow,
// the convergence test) the scheme reduces to its explicit 2-stage RK2 part, so
// the spatial order k+1 and the temporal order 2 are both clean to verify.
//
// Robustness for strong shocks (Mach-10 double Mach reflection) is provided by a
// Zhang-Shu positivity-preserving limiter applied to every stage, keeping the
// density and pressure positive at all quadrature points while preserving the
// cell average (and hence conservation) and high-order accuracy in smooth zones.
//
// Conservative variables live in the broken P_k space (elem2dof from
// FEM::getDOF, contiguous per-element blocks).  The whole state is one nDof x 4
// matrix; per-element coefficient blocks are U.middleRows(t*locDof, locDof).
// ===========================================================================

namespace euler {

using namespace Eigen;

constexpr double GAMMA = 1.4;

// --------------------------------------------------------------------------
// Pointwise gas dynamics on a conservative state U = (rho, rho u, rho v, E).
// --------------------------------------------------------------------------
inline double pressure(const Vector4d& U) {
    double rho = U(0);
    double ke  = 0.5 * (U(1) * U(1) + U(2) * U(2)) / rho;
    return (GAMMA - 1.0) * (U(3) - ke);
}
inline double pressure(const Vector4d& U, double rho_floor, double p_floor) {
    double rho = std::max(U(0), rho_floor);
    double ke  = 0.5 * (U(1) * U(1) + U(2) * U(2)) / rho;
    return std::max((GAMMA - 1.0) * (U(3) - ke), p_floor);
}
inline double soundSpeed(const Vector4d& U) {
    return std::sqrt(GAMMA * std::max(pressure(U), 1e-300) / std::max(U(0), 1e-300));
}

// Normal physical flux  F_n = F_x n_x + F_y n_y  (n need not be unit; pass unit).
inline Vector4d normalFlux(const Vector4d& U, double nx, double ny) {
    double rho = U(0), u = U(1) / rho, v = U(2) / rho, E = U(3);
    double p   = pressure(U);
    double un  = u * nx + v * ny;
    Vector4d Fn;
    Fn(0) = rho * un;
    Fn(1) = U(1) * un + p * nx;     // rho u * un + p nx
    Fn(2) = U(2) * un + p * ny;     // rho v * un + p ny
    Fn(3) = (E + p) * un;
    return Fn;
}

// x- and y- physical fluxes (used by the volume term  int F : grad phi).
inline void fluxes(const Vector4d& U, Vector4d& Fx, Vector4d& Fy) {
    double rho = U(0), u = U(1) / rho, v = U(2) / rho, E = U(3);
    double p   = pressure(U);
    Fx(0) = U(1);                 Fy(0) = U(2);
    Fx(1) = U(1) * u + p;         Fy(1) = U(1) * v;
    Fx(2) = U(2) * u;             Fy(2) = U(2) * v + p;
    Fx(3) = (E + p) * u;          Fy(3) = (E + p) * v;
}

// Maximal signal speed normal to (nx,ny):  |u.n| + c.
inline double maxWaveSpeed(const Vector4d& U, double nx, double ny) {
    double rho = std::max(U(0), 1e-300);
    double un  = (U(1) * nx + U(2) * ny) / rho;
    return std::abs(un) + soundSpeed(U);
}

// Rusanov (local Lax-Friedrichs) numerical normal flux between interior Um and
// exterior Up across a face with unit outward normal n.
inline Vector4d rusanov(const Vector4d& Um, const Vector4d& Up, double nx, double ny) {
    double lam = std::max(maxWaveSpeed(Um, nx, ny), maxWaveSpeed(Up, nx, ny));
    return 0.5 * (normalFlux(Um, nx, ny) + normalFlux(Up, nx, ny)) - 0.5 * lam * (Up - Um);
}

// HLLC (Harten-Lax-van Leer-Contact) numerical normal flux.  Unlike Rusanov,
// which damps EVERY wave by |u.n|+c, HLLC restores the contact/shear wave (speed
// S*), so it adds almost no dissipation to slip lines / shear layers -- essential
// for resolving the double-Mach-reflection slip-line Kelvin-Helmholtz vortex
// street.  Davis wave-speed estimates; reduces to the exact flux in smooth flow,
// so it preserves the design order.  (Toro, "Riemann Solvers", ch.10.)
inline Vector4d hllc(const Vector4d& UL, const Vector4d& UR, double nx, double ny) {
    double rL = UL(0), uL = UL(1) / rL, vL = UL(2) / rL, EL = UL(3), pL = pressure(UL);
    double rR = UR(0), uR = UR(1) / rR, vR = UR(2) / rR, ER = UR(3), pR = pressure(UR);
    double cL = std::sqrt(GAMMA * std::max(pL, 1e-300) / rL);
    double cR = std::sqrt(GAMMA * std::max(pR, 1e-300) / rR);
    double unL = uL * nx + vL * ny, unR = uR * nx + vR * ny;
    double SL = std::min(unL - cL, unR - cR);          // Davis wave-speed estimates
    double SR = std::max(unL + cL, unR + cR);
    Vector4d FnL = normalFlux(UL, nx, ny), FnR = normalFlux(UR, nx, ny);
    if (SL >= 0.0) return FnL;
    if (SR <= 0.0) return FnR;
    double den = rL * (SL - unL) - rR * (SR - unR);
    if (std::abs(den) < 1e-300) return rusanov(UL, UR, nx, ny);
    double Sstar = (pR - pL + rL * unL * (SL - unL) - rR * unR * (SR - unR)) / den;
    // star state on side K (normal velocity -> S*, tangential velocity unchanged)
    auto star = [&](double r, double u, double v, double E, double p, double un, double S) -> Vector4d {
        double rs = r * (S - un) / (S - Sstar);
        double us = u + (Sstar - un) * nx, vs = v + (Sstar - un) * ny;
        double Es = rs * (E / r + (Sstar - un) * (Sstar + p / (r * (S - un))));
        return Vector4d(rs, rs * us, rs * vs, Es);
    };
    if (Sstar >= 0.0) return FnL + SL * (star(rL, uL, vL, EL, pL, unL, SL) - UL);
    else              return FnR + SR * (star(rR, uR, vR, ER, pR, unR, SR) - UR);
}

// Convert primitive (rho,u,v,p) to conservative (rho,rho u,rho v,E).
inline Vector4d primToCons(double rho, double u, double v, double p) {
    return Vector4d(rho, rho * u, rho * v, p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v));
}

// Boundary ghost-state callback supplying the EXTERIOR Rusanov state on a
// boundary edge: (x, y, t, interior state, outward unit normal, edge tag).
using ExteriorStateFn =
    std::function<Vector4d(double, double, double, const Vector4d&, double, double, int)>;

// --------------------------------------------------------------------------
// Solver configuration (all overridable from JSON in the drivers).
// --------------------------------------------------------------------------
struct EulerConfig {
    int    time_order = 2;        // 1 = IMEX-Euler, 2 = ARS(2,2,2)
    bool   use_hllc = true;       // numerical flux: true = HLLC (resolves contacts/slip
                                  // lines -> crisp vortex street), false = Rusanov/LLF
    bool   use_av = true;         // Persson-Peraire artificial viscosity (implicit)
    bool   use_positivity = true; // Zhang-Shu positivity-preserving limiter
    bool   use_tvb = false;       // optional minmod/TVB slope limiter (P1 hard shocks)
    double tvb_M = 0.0;           // TVB constant (0 = pure minmod)
    // artificial viscosity
    double av_c     = 1.0;        // magnitude:  eps0_K = av_c * (h_K/k) * (|u|+c)_K
    int    av_indicator = 1;      // sensor variable: 0 = density, 1 = pressure (pressure
                                  // fires on shocks but NOT contacts/slip lines, so it
                                  // keeps the slip-line Kelvin-Helmholtz roll-ups sharp)
    double av_kappa = 1.0;        // sensor ramp half-width (decades)
    double av_s0    = -3.0;       // sensor threshold override (used iff av_s0_set)
    bool   av_s0_set = false;     // false -> auto s0 = -4 log10(k) (k>=2), -3 (k=1)
    int    av_refresh = 8;        // steps between AV-operator (eps) refreshes
    double sigma_ip = 8.0;        // SIPG penalty factor:  sigma = sigma_ip (k+1)^2
    // positivity floors
    double rho_floor = 1e-12;
    double p_floor   = 1e-12;
    // local time stepping (LTS): coarse cells advance with a larger step, fine
    // (shock + finest) cells subcycle; a flux register keeps coarse/fine
    // interfaces exactly conservative.  See stepLTS.
    bool   use_lts = false;
};

// --------------------------------------------------------------------------
// The DG-Euler integrator: holds the broken-P_k space data, the artificial-
// viscosity machinery, and advances the 4-field state with the IMEX-RK scheme.
// --------------------------------------------------------------------------
class EulerDG {
public:
    EulerDG(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
            const MatrixXi& edge, const MatrixXi& edge2side,
            const VectorXi& edgeTag, const EulerConfig& cfg);

    void setState(const MatrixXd& U0);              // nDof x 4
    const MatrixXd& state() const { return U_; }
    int nDof() const { return nDof_; }
    int stepsTaken() const { return nstep_; }

    // Advance one step of size dt to physical time tEnd; bc gives ghost states.
    // Returns false on a failed implicit solve.
    bool step(double dt, double tEnd, const ExteriorStateFn& bc);

    // LOCAL-TIME-STEPPING macro step (n-level recursive subcycling).  cellClass[t]
    // in {0..levels-1}: class 0 = finest (step dt0 = dtMacro/2^(levels-1)), class c
    // takes dt0*2^c, the coarsest (class levels-1) takes dtMacro in one step.  Each
    // coarse/fine class interface carries a flux register so mass/momentum/energy
    // are conserved to machine precision.  Caller must guarantee: (i) dt0 satisfies
    // each cell's CFL at its class; (ii) face-adjacent classes differ by <=1 (2:1
    // balanced); (iii) all AV cells (eps>0) and their neighbours are class 0 (the
    // implicit operator never crosses an interface).
    bool stepLTS(double dtMacro, double tEnd, const ExteriorStateFn& bc,
                 const std::vector<int>& cellClass, int levels);

    // Global maximal wave speed (for a CFL-based dt).
    double maxWaveSpeedGlobal() const;

    // DEBUG: min density & pressure over all flux quadrature points (volume +
    // both edge orientations) and the worst element index + its centroid.
    void quadMin(double& rmin, double& pmin, int& worstElem, double& cx, double& cy) const;
    // DEBUG: cell-average conservative state of element t.
    Vector4d cellMeanState(int t) const { return cellAverage(U_, t); }

    // Integral of each conserved variable over the domain (conservation check).
    Vector4d conservedTotals() const;

    // Nodal fields for visualisation / error measurement (length nDof).
    VectorXd densityField()  const { return U_.col(0); }
    VectorXd pressureField() const;
    VectorXd machField()     const;
    VectorXd velocityMagField() const;
    VectorXd schlierenField() const;   // exp(-k |grad rho| / max) shading, per element
    VectorXd avField()       const;    // last per-element eps broadcast to nodes

    // Inviscid spatial residual  R = int F:grad phi - oint Hhat phi  (nDof x 4),
    // i.e. dU/dt|_explicit = M^{-1} R.  Exposed for the convergence driver / tests.
    void inviscidResidual(const MatrixXd& U, double t, const ExteriorStateFn& bc,
                          MatrixXd& R) const;
    // Apply the (block-diagonal, affine) inverse mass matrix in place: R_e <- Minv_e R_e.
    void applyMassInverse(MatrixXd& R) const;

private:
    // ---- LTS helpers ----
    // COMPACT masked inviscid residual: Uc and Rc are packed (nc*locDof x 4), compact
    // row ci*locDof+i <-> class cell cells[ci], local dof i.  A SAME-class neighbour is
    // read from the compact working array Uc (via g2c_); another class reads the
    // macro-start snapshot ltsSnap_ (global); boundary faces use bc.
    void inviscidResidualMaskedCompact(const MatrixXd& Uc, const std::vector<int>& cells,
                                       double t, const ExteriorStateFn& bc, MatrixXd& Rc) const;
    // positivity limiter restricted to a cell list (the rest are untouched).
    void positivityLimitCells(MatrixXd& U, const std::vector<int>& cells) const;
    // Flux register: accumulate the time-integrated interface flux of the advancing
    // class `advClass` (against the COARSE cell's test functions) into reg[coarseClass]
    // for every class interface touching advClass.  weight = stage_b * dt.  The
    // advancing cell's state comes from the COMPACT stage array Uc (via g2c_); the
    // other side from ltsSnap_ (global).
    void accumulateReflux(const MatrixXd& Uc, int advClass, double weight,
                          std::vector<MatrixXd>& reg) const;
    // One masked ARS(2,2,2) substep of class `advClass` (implicit AV only for class 0);
    // accumulates its interface fluxes into the registers.
    bool advanceMaskedARS(const std::vector<int>& cells, int advClass, double dt,
                          double tStageBase, const ExteriorStateFn& bc,
                          std::vector<MatrixXd>& reg);
    // Recursive subcycling: advance the hierarchy one step of class `level`'s clock.
    void ltsIntegrate(int level, double dt, double tStageBase, const ExteriorStateFn& bc,
                      std::vector<MatrixXd>& reg, int levels);

    // ---- precompute (reference basis traces, mass, projection) ----
    void precompute();
    // ---- artificial viscosity ----
    void computeViscosity(const MatrixXd& U, VectorXd& epsK) const;   // sensor -> eps_K
    SparseMatrix<double> assembleAV(const VectorXd& epsK) const;      // SIPG variable-coeff -Lap (SPD)
    void refreshImplicit(const VectorXd& epsK, double a_dt);          // build A, K_aa, factor
    // solve (M + a_dt A) X = B (X,B are nDof x 4) via the partitioned scheme
    MatrixXd implicitSolve(const MatrixXd& B) const;
    void applyMass(const MatrixXd& U, MatrixXd& MU) const;            // MU = M U (block, parallel)
    // ---- limiters ----
    void limitCellCore(Matrix<double, Dynamic, 4>& Ue) const;   // Zhang-Shu squeeze of one locDof x 4 block
    void limitCell(MatrixXd& U, int t) const;   // gather cell t -> core -> scatter (global)
    void positivityLimit(MatrixXd& U) const;
    void tvbLimit(MatrixXd& U) const;
    // ---- helpers ----
    Vector4d cellAverage(const MatrixXd& U, int t) const;
    void elementBlock(const MatrixXd& U, int t, Matrix<double, Dynamic, 4>& Ue) const;

    FEM& fem_;
    Mesh& mesh_;
    const MatrixXi& elem2dof_;
    const MatrixXi& edge_;
    const MatrixXi& edge2side_;
    VectorXi edgeTag_;
    EulerConfig cfg_;

    int NT_, NE_, locDof_, nDof_, nstep_;
    double hmin_;

    // reference-element data
    MatrixXd MrefInv_;              // (locDof x locDof) inverse of unit-area mass
    VectorXd hK_;                   // per-element size = longest edge (P-P AV scale)
    SparseMatrix<double> Mblk_;     // block-diagonal DG mass (for implicit LHS)

    // implicit operator (K = M + gamma*dt*A) solved by a PARTITIONED scheme: the
    // artificial viscosity A is nonzero only on shock cells + their neighbours, so
    // K is EXACTLY block-diagonal between those "active" dofs and the rest.  The
    // decoupled rows are just the block mass M (trivial parallel inverse); only the
    // small active sub-block K_aa is factorised (Cholesky) and solved.  A_ is kept
    // for the explicit A*U2 term.  All refreshed once per av_refresh window.
    SparseMatrix<double> A_;                          // assembled artificial-viscosity operator
    SparseMatrix<double> Kaa_;                        // active sub-block of K (= M+gamma dt A)
    SparseMatrix<double> Aaa_;                        // A restricted to active x active (active indexing; for compact A*U2)
    SimplicialLDLT<SparseMatrix<double>> ldltA_;      // factor of Kaa_
    std::vector<int> activeDofs_;                     // global dofs touched by A
    bool implicitReady_ = false;
    double aDtCached_ = -1.0;       // the gamma*dt baked into Kaa_ / A_ scaling
    VectorXd epsFrozen_;            // viscosity used in the cached operator
    int refreshCountdown_ = 0;

    // opaque precomputed quadrature / trace tables (defined in the .cpp)
    struct Impl;
    std::shared_ptr<Impl> P_;

    MatrixXd U_;                    // nDof x 4 conservative state
    std::vector<int> cellClass_;    // LTS: per-element dt-class 0..levels-1 (set in stepLTS)
    std::vector<int> g2c_;          // LTS: global cell -> compact index for the advancing class (-1 otherwise)
    MatrixXd ltsSnap_;              // LTS: macro-start snapshot (cross-class ghost; no future-read)
    // LTS flux-register interface: each class-interface shared edge ("coarse" = the
    // higher-class / bigger-dt cell), with the local edge index (0..2) on each side.
    struct LtsEdge { int coarse, fine, coarseK, fineK; };
    std::vector<LtsEdge> ltsEdges_;
    std::vector<std::vector<int>> classCells_;   // compact per-class cell lists (set in stepLTS)
};

// --------------------------------------------------------------------------
// Structured right-triangle mesh of the rectangle [x0,x1] x [y0,y1] with nx x ny
// cells (each split into two CCW triangles).  Used for the vortex and DMR domains.
void makeRectMesh(Mesh& mesh, double x0, double x1, double y0, double y1, int nx, int ny);

// L2-project an initial primitive field (rho,u,v,p)(x,y) into the DG space,
// returning conservative coefficients (nDof x 4).
// --------------------------------------------------------------------------
MatrixXd projectInitial(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                        const std::function<Vector4d(double, double)>& primField);

// L2 error of a conserved component against an exact conservative field.
double l2Error(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
               const MatrixXd& U, int comp,
               const std::function<Vector4d(double, double)>& exactCons);

// --------------------------------------------------------------------------
// Visualisation: rasterise a scalar dP_k field onto an image and write a PPM,
// sampling the true high-order polynomial per pixel (smooth interior shading).
// --------------------------------------------------------------------------
enum Colormap { CM_VIRIDIS, CM_INFERNO, CM_GRAY, CM_GRAY_INV, CM_COOLWARM, CM_JET };
std::vector<unsigned char> renderScalarPPMImage(
                    FEM& fem, const Mesh& mesh,
                    const MatrixXi& elem2dof, const VectorXd& field, int W, int H,
                    double xmin, double xmax, double ymin, double ymax,
                    double vmin, double vmax, Colormap cm,
                    const std::function<bool(double, double)>& inDomain = {});
void writeScalarPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                    const MatrixXi& elem2dof, const VectorXd& field, int W, int H,
                    double xmin, double xmax, double ymin, double ymax,
                    double vmin, double vmax, Colormap cm,
                    const std::function<bool(double, double)>& inDomain = {});

// Numerical Schlieren: per-pixel |grad rho| from the TRUE dP_k polynomial gradient
// (not piecewise-constant), shaded  s = exp(-beta |grad rho| / max|grad rho|) in
// grayscale (white = smooth, dark = steep).  max|grad rho| is taken over the render
// window so the image auto-contrasts to its local features (e.g. a zoom on the
// slip line shows the Kelvin-Helmholtz roll-ups, not just the main shocks).
void writeSchlierenPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                       const MatrixXi& elem2dof, const VectorXd& rho, int W, int H,
                       double xmin, double xmax, double ymin, double ymax, double beta = 8.0);

} // namespace euler

#endif
