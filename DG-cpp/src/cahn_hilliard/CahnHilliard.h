#ifndef CAHN_HILLIARD_H
#define CAHN_HILLIARD_H

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <string>
#include <functional>
#include <map>
#include "FEM.h"
#include "Mesh.h"

using namespace Eigen;

// ---------------------------------------------------------------------------
// Cahn-Hilliard mixed-DG helpers.
//
// We solve, on Omega = [0,1]^2 with no-flux (homogeneous Neumann) BCs:
//      dc/dt = M_mob * Laplacian(mu),     mu = F'(c) - eps^2 * Laplacian(c),
// with the quartic double-well potential
//      F(c) = (1/4)(c^2 - 1)^2,   F'(c) = c^3 - c,   F''(c) = 3 c^2 - 1.
//
// The Laplacian is discretised with the existing symmetric interior-penalty DG
// (SIPG) bilinear form a(.,.) -> matrix A = assembleK_Poi2D + assembleIP_Poi2D
// (beta = 1).  Because the interior-penalty term is added on interior edges
// only, a(.,.) carries no boundary contribution, which is exactly the natural
// no-flux boundary condition; A is symmetric PSD with the constant vector in
// its nullspace (A * 1 = 0), so discrete mass is conserved.
// ---------------------------------------------------------------------------

// DG mass matrix  Mm(i,j) = sum_K \int_K phi_i phi_j   (block-diagonal per element,
// SPD).  Uses the element nodal Lagrange basis from FEM.
SparseMatrix<double> assembleMass_DG2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof);

// Nonlinear (explicit) load  b(i) = sum_K \int_K F'(c_h) phi_i  with F'(c) = c^3 - c,
// where c_h is the current DG field given by the coefficient vector c.
VectorXd assembleNonlinearCH(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                             const VectorXd& c);

// Discrete Cahn-Hilliard energy  E_h(c) = (eps^2/2) c^T A c + sum_K \int_K F(c_h).
double computeEnergyCH(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                       const VectorXd& c, const SparseMatrix<double>& A, double eps2);

// Total mass  = 1^T Mm c = \int_Omega c_h  (should stay constant in time).
double computeMassCH(const SparseMatrix<double>& Mm, const VectorXd& c);

// Random spinodal-decomposition initial condition: per-DOF
//      c(dof) = cbar + amp * U(-1, 1),
// reproducible through the supplied seed.  Returns a coefficient vector of length nDof.
VectorXd initSpinodal(int nDof, double cbar, double amp, unsigned seed);

// Rasterise the DG field c_h onto an Npix x Npix grid over the mesh bounding box
// and write it as a binary (P6) PPM image.  Values are mapped through a diverging
// "coolwarm" colormap clamped to [cmin, cmax].  Npix should be EVEN so the frames
// can be encoded with libx264 / yuv420p.  Image rows run top (max y) to bottom.
void writeFramePPM(const std::string& path, const Mesh& mesh,
                   const MatrixXi& elem2dof, const VectorXd& c,
                   int Npix, double cmin, double cmax);

// Generic load vector  L(i) = sum_K \int_K f(x,y) phi_i.  Used to assemble the
// manufactured-solution forcing term for the convergence test.
VectorXd assembleScalarLoad(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                            const std::function<double(double, double)>& f);

// ---------------------------------------------------------------------------
// Time integrator for the mixed-DG Cahn-Hilliard system.  Encapsulates the
// constant 2x2 block matrix and its one-time LU factorisation, and advances the
// field one step with a chosen time order:
//   timeOrder = 1 : first-order stabilised semi-implicit (backward-Euler type),
//   timeOrder = 2 : second-order SBDF2 (BDF2 + 2nd-order extrapolation of the
//                   explicit nonlinear term), bootstrapped by ONE first-order step.
// Both keep the system matrix constant in time, so it is factorised once (or
// twice for order 2: the bootstrap J1 and the SBDF2 J2).
// ---------------------------------------------------------------------------
class CHIntegrator {
public:
    CHIntegrator(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                 const SparseMatrix<double>& Mm, const SparseMatrix<double>& A,
                 double tau, double mob, double eps2, double S, int timeOrder);

    // Sets c^0 and resets the step counter.
    void setInitial(const VectorXd& c0);

    // Advances one step (t_n -> t_{n+1}).  `sourceLoad` is the optional forcing
    // vector \int g(.,t_{n+1}) phi (size nDof), added to the c-equation RHS; pass
    // an empty vector for the unforced problem.  Returns false on solve failure.
    bool step(const VectorXd& sourceLoad);
    bool step() { return step(VectorXd()); }

    const VectorXd& field() const { return c_; }
    int stepsTaken() const { return n_; }
    int timeOrder() const { return order_; }

private:
    void buildAndFactor(double a0, SparseMatrix<double>& J,
                        SparseLU<SparseMatrix<double>>& lu);

    FEM& fem_;
    const Mesh& mesh_;
    const MatrixXi& elem2dof_;
    const SparseMatrix<double>& Mm_;
    const SparseMatrix<double>& A_;
    double tau_, mob_, eps2_, S_;
    int order_;
    int nDof_;
    int n_;                 // number of steps taken
    VectorXd c_, cprev_;    // c^n and c^{n-1}
    SparseMatrix<double> J1_, J2_;
    SparseLU<SparseMatrix<double>> lu1_, lu2_;  // J1 (order 1 / bootstrap), J2 (SBDF2)
};

// ---------------------------------------------------------------------------
// Adaptive-step first-order (backward-Euler) integrator for the movie.  Each step
// uses a tau chosen from the relative solution change (small while the field changes
// fast, large once it coarsens slowly).  tau is snapped to a power-of-2 level grid in
// [tauMin, tauMax]; the tau-dependent block matrix is factorised at most ONCE per
// level and cached, so only a handful of factorisations happen over a whole run.
// Backward Euler is unconditionally energy-stable and conserves mass exactly at ANY
// tau, so adaptivity preserves both invariants.
// ---------------------------------------------------------------------------
class CHAdaptiveStepper {
public:
    CHAdaptiveStepper(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                      const SparseMatrix<double>& Mm, const SparseMatrix<double>& A,
                      double mob, double eps2, double S, double tauMin, double tauMax);

    void setInitial(const VectorXd& c0);

    // One backward-Euler step at a tau near `tauDesired` (snapped to the level grid;
    // factorisation built lazily and cached).  `src` is an optional forcing vector
    // added to the c-equation RHS.  Returns false on solve failure.
    bool step(double tauDesired, const VectorXd& src = VectorXd());

    double lastTau() const { return lastTau_; }            // tau actually used last step
    double lastRelChange() const { return lastRel_; }      // ||dc||_Mm / ||c||_Mm last step
    int numFactorizations() const { return static_cast<int>(cache_.size()); }
    const VectorXd& field() const { return c_; }

private:
    int levelFor(double tau) const;
    SparseLU<SparseMatrix<double>>* solverForLevel(int lvl); // lazy build + cache

    FEM& fem_;
    const Mesh& mesh_;
    const MatrixXi& elem2dof_;
    const SparseMatrix<double>& Mm_;
    const SparseMatrix<double>& A_;
    double mob_, eps2_, S_, tauMin_, tauMax_;
    int nDof_, maxLevel_;
    VectorXd c_;
    double lastTau_, lastRel_;
    std::map<int, SparseLU<SparseMatrix<double>>> cache_;  // level -> factorised J(tau)
};

#endif
