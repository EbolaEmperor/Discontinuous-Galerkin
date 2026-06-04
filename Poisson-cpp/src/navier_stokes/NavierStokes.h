#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
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

struct BCData {
    VectorXi bcU, bcV, bcP;     // length = #edges
    // Boundary data as functions of (x, y, t):
    std::function<Vector2d(double, double, double)> velDir;  // Dirichlet velocity value
    std::function<Vector2d(double, double, double)> velAcc;  // d/dt of Dirichlet velocity
    std::function<double (double, double, double)>  presDir; // Dirichlet pressure value
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
                 const BCData& bc, double nu, double dt, double sigma, double beta);

    // Set the initial velocity (and optional pressure) coefficient vectors.
    void setInitial(const VectorXd& u0, const VectorXd& v0,
                    const VectorXd& p0 = VectorXd());

    // Advance one step, ending at physical time tEnd (= t^{n+1}); first step is
    // bootstrapped with first-order BDF1/EX1.  Returns false on solve failure.
    bool step(double tEnd);

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
    int nDof_, n_;

    SparseMatrix<double> M_, Gx_, Gy_;
    SparseMatrix<double> Au_, Av_, Ap_;       // SIPG -Lap with the three BC sets
    SparseMatrix<double> Hu_, Hv_;            // gamma0/dt*M + nu*A (BDF2 gamma0=3/2)
    SimplicialLDLT<SparseMatrix<double>> luM_, luHu_, luHv_, luAp_;
    // Helmholtz matrix uses gamma0=3/2; the very first (BDF1) step needs gamma0=1.
    SparseMatrix<double> Hu1_, Hv1_;
    SimplicialLDLT<SparseMatrix<double>> luHu1_, luHv1_;

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

} // namespace ns

#endif
