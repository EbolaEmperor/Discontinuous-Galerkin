#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include "FEM.h"
#include "Mesh.h"

using namespace Eigen;

// ===========================================================================
// Discontinuous-Galerkin (dP_k) solver for the 2-D incompressible Navier-Stokes
// equations
//        du/dt + (u.grad)u = -grad p + nu*Lap(u),   div u = 0,
// using a high-order stiffly-stable splitting (Karniadakis-Israeli-Orszag) in
// time: 2nd-order BDF2 for the time derivative, 2nd-order extrapolation (EX2)
// of the explicit convection, and a pressure-Poisson projection.  Each time
// step therefore costs only scalar elliptic solves whose system matrices are
// constant in time (factorised once):
//   - one pressure Poisson solve (SIPG -Lap, Dirichlet p where outflow),
//   - two velocity Helmholtz solves (gamma0/dt*M + nu*(SIPG -Lap)).
// The convection is treated explicitly with a local Lax-Friedrichs flux.
//
// All scalar fields (u, v, p) live in the same broken P_k space (elem2dof from
// FEM::getDOF).  The viscous/pressure operators reuse the existing SIPG volume
// + interior-penalty assembly (assembleK_Poi2D + assembleIP_Poi2D); boundary
// conditions are imposed weakly (Nitsche for Dirichlet, natural for Neumann).
// ===========================================================================

namespace ns {

// Per-edge boundary-condition codes (index aligned with the rows of `edge`).
//   velocity (bcU, bcV):  0 = interior, 1 = Dirichlet (Nitsche), 2 = Neumann (natural)
//   pressure (bcP):       0 = interior, 1 = Dirichlet (Nitsche, value presDir),
//                         2 = Neumann high-order (KIO consistent BC),
//                         3 = Neumann homogeneous (dp/dn = 0, slip-wall symmetry)
enum { BCN_INTERIOR = 0, BCN_DIRICHLET = 1, BCN_NEUMANN = 2, BCN_NEUMANN_HO = 2, BCN_NEUMANN_ZERO = 3 };
// Pressure equation modes:
//   projection:  Delta p = (1/dt) div(ha - dt N*)  (classic split projection PPE)
//   direct PPE:  Delta p = -div N* with KIO/Gresho-Sani pressure Neumann data.
enum { NSPRESSURE_PROJECTION = 0, NSPRESSURE_DIRECT_PPE = 1 };

struct BCData {
    VectorXi bcU, bcV, bcP;     // length = #edges
    // Boundary data as functions of (x, y, t):
    std::function<Vector2d(double, double, double)> velDir;  // Dirichlet velocity value
    std::function<Vector2d(double, double, double)> velAcc;  // d/dt of Dirichlet velocity
    std::function<double (double, double, double)>  presDir; // Dirichlet pressure value
    // DIAGNOSTIC ONLY: if set, the high-order pressure Neumann data on bcP==2 edges
    // uses the exact pressure gradient  dp/dn = n.presGrad  instead of the KIO
    // reconstruction n.(-N* + nu*Lap u* - a_b).  Used by the convergence study to
    // isolate the reconstruction error; leave null for real runs.
    std::function<Vector2d(double, double, double)> presGrad;
};

// Scalar DG mass matrix  M(i,j) = sum_K \int_K phi_i phi_j  (block-diagonal SPD).
SparseMatrix<double> assembleScalarMassDG(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof);

// Boundary Nitsche matrix for the SIPG -Laplacian on the edges flagged Dirichlet
// (isDir[e] != 0):   sum_e \int_e [ -(grad u.n) v - beta (grad v.n) u + (sigma/h) u v ].
SparseMatrix<double> assembleNitscheDirichlet(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                              const MatrixXi& edge, const MatrixXi& edge2side,
                                              const VectorXi& isDir, double sigma, double beta);

// DG weak first-derivative operators Gx, Gy (each nDof x nDof) with central
// interior fluxes and single-sided boundary traces, such that
//   (Gm p)_i = \int_Omega (d_m p) phi_i     (m = x, y),
// and the DG divergence of a vector field w is  Dvec(w) = Gx*w_x + Gy*w_y, while
// the DG gradient load of p is  (Gx p, Gy p).  These dual roles share the same
// matrices, so divergence and gradient are exact adjoints.
void assembleWeakGrad(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                      const MatrixXi& edge, const MatrixXi& edge2side,
                      SparseMatrix<double>& Gx, SparseMatrix<double>& Gy);

// ---------------------------------------------------------------------------
// Time integrator.  Holds the constant operators and their factorisations and
// advances (u,v,p) one step.  Boundary data is supplied through `bc`.
// ---------------------------------------------------------------------------
class NSIntegrator {
public:
    NSIntegrator(FEM& fem, Mesh& mesh, const MatrixXi& elem2dof,
                 const MatrixXi& edge, const MatrixXi& edge2side,
                 const BCData& bc, double nu, double dt, double sigma, double beta,
                 // gradDiv>0 enables optional coupled velocity grad-div; ppeDivDamping
                 // damps divergence through the Direct PPE RHS without coupling u/v.
                 double gradDiv = 0.0, int pressureMode = NSPRESSURE_PROJECTION,
                 double ppeDivDamping = 0.0);

    // Set the initial velocity (and optional pressure) coefficient vectors.
    void setInitial(const VectorXd& u0, const VectorXd& v0,
                    const VectorXd& p0 = VectorXd());

    // Advance one step, ending at physical time tEnd (= t^{n+1}); first step is
    // bootstrapped with first-order BDF1/EX1.  Returns false on solve failure.
    bool step(double tEnd);

    // Same as step(), but also adds an external nodal body-force load
    // (already integrated against the basis: each entry has units of  int f phi).
    // Used by the immersed-boundary FSI driver to inject the constraint forces
    // produced by the Cosserat filament.  Pass empty VectorXd to skip a side.
    bool stepWithBodyForce(double tEnd, const VectorXd& loadFx, const VectorXd& loadFy);

    // High-order pressure Neumann data uses the *rotational* form of the viscous
    // term,  n.(-nu*curl(omega)),  rather than the Laplacian form  n.(nu*Lap u).
    // The two agree for a divergence-free field, but the discrete velocity is not
    // exactly divergence-free, and the Laplacian form leaks the spurious term
    // n.grad(div u) into the boundary pressure -> an O(dt) splitting error that
    // caps the temporal order at 1.  The rotational form restores 2nd order
    // (E & Liu 1995; Karniadakis-Israeli-Orszag 1991).  Default true (= the fix);
    // expose it only so the order-reduction can be reproduced for verification.
    bool rotationalPressureBC = true;

    // DIAGNOSTIC ONLY: when false, drop the explicit convection entirely (so the
    // solver advances unsteady Stokes).  Lets the temporal-order study isolate the
    // convection treatment from the pressure/viscous splitting.  Default true.
    bool includeConvection = true;

    const VectorXd& u() const { return u_; }
    const VectorXd& v() const { return v_; }
    const VectorXd& p() const { return p_; }
    int stepsTaken() const { return n_; }

    // Vorticity field  omega = dv/dx - du/dy  (L2 projection onto the DG space).
    VectorXd vorticity() const;
    // Velocity magnitude field |u| at the DG nodes (for visualisation).
    VectorXd speed() const;

    // Net force on the cylinder boundary (BD_CYL edges, identified by bcU==Dirichlet
    // AND tagged as cylinder via `isCylEdge`).  Returns (Fx, Fy) = \oint (-p n + nu(grad u+grad u^T) n).
    void cylinderForce(const VectorXi& isCylEdge, double& Fx, double& Fy) const;

private:
    void buildOperators();
    // Explicit Lax-Friedrichs convection load c_m = \int N_m(u) phi  (m=x,y).
    void assembleConvection(const VectorXd& uu, const VectorXd& vv,
                            double t, VectorXd& cx, VectorXd& cy) const;
    // High-order pressure Neumann load  \int_{Gamma} g_N q  with
    //   g_N = n.[ -N* + nu*Lap(u*) - a_b ]  on bcP==2 edges.
    VectorXd assemblePressureNeumann(const VectorXd& us, const VectorXd& vs,
                                     double tEnd) const;
    // Dirichlet RHS lift for a velocity component (Nitsche), comp=0->u,1->v.
    VectorXd velDirichletLift(int comp, double tEnd) const;
    // Dirichlet RHS lift for pressure (Nitsche) on bcP==1 edges.
    VectorXd presDirichletLift(double tEnd) const;
    // Solve M x = b per component (block-diagonal SPD mass).
    VectorXd massSolve(const VectorXd& b) const;

    FEM& fem_;
    Mesh& mesh_;
    const MatrixXi& elem2dof_;
    const MatrixXi& edge_;
    const MatrixXi& edge2side_;
    BCData bc_;
    double nu_, dt_, sigma_, beta_;
    double gradDiv_, ppeDivDamping_;
    int pressureMode_;
    int nDof_, n_;

    SparseMatrix<double> M_, Gx_, Gy_;
    SparseMatrix<double> Au_, Av_, Ap_;       // SIPG -Lap with the three BC sets
    SparseMatrix<double> Hu_, Hv_;            // gamma0/dt*M + nu*A (BDF2 gamma0=3/2)
    SimplicialLDLT<SparseMatrix<double>> luM_, luHu_, luHv_, luAp_;
    // Helmholtz matrix uses gamma0=3/2; the very first (BDF1) step needs gamma0=1.
    SparseMatrix<double> Hu1_, Hv1_;
    SimplicialLDLT<SparseMatrix<double>> luHu1_, luHv1_;
    SparseMatrix<double> Huv_, Huv1_;         // coupled velocity matrices with grad-div
    SimplicialLDLT<SparseMatrix<double>> luHuv_, luHuv1_;

    VectorXd u_, v_, p_;          // current fields  (t^n)
    VectorXd uPrev_, vPrev_;      // previous fields (t^{n-1})
    VectorXd cxPrev_, cyPrev_;    // previous convection load (for EX2)
};

// ---------------------------------------------------------------------------
// Visualisation: rasterise a scalar DG field onto an Npix-wide image over a
// chosen window and write a binary PPM through a colormap.  `mask` (the geometry
// signed-distance) leaves the cylinder hole / exterior as background.
// ---------------------------------------------------------------------------
enum Colormap { CM_COOLWARM, CM_VIRIDIS };
// The field is sampled at each pixel through the true dP_k basis (smooth, high-order
// rendering), not just the vertex values.
void writeFieldPPM(const std::string& path, FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                   const VectorXd& field, int Wpix, int Hpix,
                   double xmin, double xmax, double ymin, double ymax,
                   double vmin, double vmax, Colormap cmap,
                   const std::function<bool(double, double)>& inDomain);

// ---------------------------------------------------------------------------
// Streakline visualisation: a swarm of massless tracer particles advected by
// the live (u,v) field.  Plotted on top of a dim |u| background, this produces
// a movie that "looks like fluid moving" -- particles stream past the cylinder
// and curl into the von Karman vortices, instead of just colour-coding vorticity.
// Implementation: mesh-locator-accelerated particle-in-element + RK2 advance,
// dead particles re-seeded at an inflow rake every step.  Frames are PPM, just
// like writeFieldPPM, so they fold into the same ffmpeg pipeline.
// ---------------------------------------------------------------------------

// Spatial index: O(1)-amortised point-in-element lookup on a static triangle
// mesh by bucketing each element's AABB into a uniform grid.
class MeshLocator {
public:
    // Build a grid covering the mesh extent; nx<=0 picks ~sqrt(NT) cells per side.
    void build(const Mesh& mesh, int nx = 0);
    // Locate the element containing (x,y).  `hint` (>=0) is checked first as
    // a warm cache; on success writes the barycentric coords (lam1,lam2,lam3)
    // matching mesh.elem(t,0..2).  Returns -1 if (x,y) is outside the mesh.
    int locate(const Mesh& mesh, double x, double y, int hint, double lam[3]) const;
private:
    double xmin_ = 0, ymin_ = 0, invDx_ = 0, invDy_ = 0;
    int nx_ = 0, ny_ = 0;
    std::vector<std::vector<int>> cells_;
};

// Population of tracer particles.  Coordinates are physical; rendering happens
// through writeParticlesPPM which uses (xa,xb,ya,yb) as the image window.
struct ParticleTracer {
    int N = 1500;
    double xa = 0, xb = 1, ya = 0, yb = 1;          // render / sim window
    double inflowX = 0;                             // re-injection rake (left edge)
    std::function<bool(double, double)> inDomain;   // alive iff inDomain(x,y)
    std::uint32_t rng = 12345u;
    // Per-particle ribbon trail (continuous streak instead of a flickering dot):
    //   trailLen    -- number of historical samples kept per particle
    //   trailStride -- only every k-th advance() call is recorded into the trail
    //                  (1 = every step; >1 stretches the visible streak without
    //                  growing memory).  Set trailLen=1 to recover plain dots.
    int trailLen    = 24;
    int trailStride = 1;
    // Per-particle state.
    std::vector<double> x, y;                          // current pos
    std::vector<int>    age;                           // frames since (re)spawn
    std::vector<char>   alive;
    std::vector<int>    hint;                          // last hosting element
    // Ring buffer of past positions per particle (size = N * trailLen).
    // trailHead[i] points at the most recently written slot for particle i.
    std::vector<double> trailX, trailY;
    std::vector<char>   trailValid;                    // 1 if slot was written since reset
    std::vector<int>    trailHead;
    int                 strideTick = 0;                // counter for trailStride

    // Allocate / scatter N particles uniformly inside the window (skipping cylinder).
    void reset();
    // Advance every particle by dt with RK2 against the supplied (u,v) DG fields,
    // reseeding dead/escaped particles at the inflow rake.
    void advance(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                 const VectorXd& uField, const VectorXd& vField,
                 const MeshLocator& loc, double dt);
};

// Render a streakline frame: a dimmed viridis(|u|) background + alpha-blended
// fading ribbon trails (no flickering dots).  `bgDim` in [0,1] mixes the |u|
// colour toward black -- 0 = pure black background, 1 = full viridis.
void writeParticlesPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                       const MatrixXi& elem2dof, const VectorXd& speedField,
                       int Wpix, int Hpix,
                       double xmin, double xmax, double ymin, double ymax,
                       double speedMin, double speedMax,
                       const ParticleTracer& tracer,
                       const std::function<bool(double, double)>& inDomain,
                       double bgDim = 0.12);

} // namespace ns

#endif
