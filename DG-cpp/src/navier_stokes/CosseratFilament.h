#ifndef COSSERAT_FILAMENT_H
#define COSSERAT_FILAMENT_H

// ===========================================================================
// Discrete Cosserat / Kirchhoff rod (a.k.a. "discrete elastic rod") for the
// 2-D flexible-tail FSI demo.
//
//   - Geometry:  N+1 nodes X_i in R^2, N edges, edge lengths l_i, edge
//                tangents phi_i = atan2(dY, dX).
//   - Energies:  stretch  V_s = (1/2) EA  sum_i ((l_i - l0_i)^2 / l0_i)
//                bending  V_b = (1/2) EI  sum_i ((phi_i - phi_{i-1})^2 / lbar_i)
//                where lbar_i = (l0_{i-1} + l0_i)/2 is the rest length of
//                the dual cell at interior node i.
//   - Equation:  M * Xddot = -grad V(X) + F_ext - gamma * M * Xdot
//                with the clamp at node 0 (and optionally at node 1 to lock
//                the wall angle) imposed exactly by "freezing" those DOFs.
//   - Stepper:   Newmark-beta (beta = 1/4, gammaN = 1/2, unconditionally
//                stable for linear problems).  At each step solve the
//                non-linear residual  R(X^{n+1}) = 0  by 1-3 Newton iterates
//                with the analytic stiffness K = d/dX (-grad V).  Sparse
//                ~7-banded SPD system, factorised by SimplicialLDLT.
//
// Units / conventions: caller picks them.  In the FSI driver everything is
// non-dimensional with rho_f = U_inf = D = 1, so EI -> K_B and EA -> K_S.
// ===========================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>

namespace ns {

class CosseratFilament {
public:
    // Number of nodes  = N+1; degrees of freedom  = 2*(N+1).
    int N = 0;

    // Material / geometric parameters.
    double EA       = 1e4;   // tensile (axial) stiffness
    double EI       = 1e-2;  // bending stiffness
    double rhoLine  = 1.0;   // line density rho_s * h  (mass per unit length)
    double dampStr  = 0.0;   // structural Rayleigh-mass damping coeff (gamma)

    // Per-node external linear damping (size N+1, drag-type):
    //   F_drag_i = -dragCoef_i * (V_i - dragRef_i)
    // Applied implicitly inside step() so it can be arbitrarily stiff.
    // Both vectors default to size 0 (= no drag).  The driver re-fills them
    // every step from the sampled fluid velocity.
    Eigen::VectorXd dragCoef;     // (N+1) per-node coefficient (already x ds)
    Eigen::MatrixXd dragRef;      // (N+1) x 2 reference velocity per node

    // Reference (rest) configuration -- straight line by default.
    Eigen::MatrixXd X0;      // (N+1) x 2

    // Current state.
    Eigen::MatrixXd X;       // (N+1) x 2 positions
    Eigen::MatrixXd V;       // (N+1) x 2 velocities
    Eigen::MatrixXd A;       // (N+1) x 2 accelerations

    // Rest lengths of edges (size N) and dual-cell lengths at interior nodes (size N+1).
    Eigen::VectorXd l0;      // size N
    Eigen::VectorXd lbar;    // size N+1 (sentinels at 0 and N).

    // Lumped nodal masses (size N+1).  For uniform rest length, m_i = rhoLine * lbar_i.
    Eigen::VectorXd m;

    // Clamp DOFs.  "Frozen" indices (component-wise: 2*node + dim) are held at
    // the prescribed clamp value during the solve.  Callers set this through
    // setClampNodes(); the default clamps node 0 (positions).
    std::vector<int> frozen;
    Eigen::VectorXd  frozenVal;   // same length as `frozen`, target positions

    // ------------------------------------------------------------------ API
    // Initialise an N-segment straight rod from (xs, ys) to (xe, ye).
    void initStraight(int N, double xs, double ys, double xe, double ye,
                      double rhoLineIn, double EAin, double EIin);

    // Clamp node 0 fully (x and y).  Optionally clamp node 1 too -- this locks
    // the "tangent at the wall", giving a true cantilever instead of pin.
    void clampRoot(bool alsoClampNode1 = true);

    // External nodal load (size 2*(N+1)).  Caller fills this every step
    // (typically = IBM constraint force * arclength weight).  Reset internally.
    Eigen::VectorXd Fext;

    // Newmark-beta step from t -> t + dt, with EA, EI, rhoLine, dampStr,
    // Fext and frozen* fixed during the step.  Returns false on Newton failure.
    bool step(double dt);

    // Diagnostics.
    double tipDeflection() const;     // |X_N - X0_N|
    double potentialEnergy() const;
    double kineticEnergy()   const;
    double maxEdgeStrain()   const;   // max |l_i - l0_i| / l0_i

private:
    // Assemble residual r = M*A_pred + grad V(X) - Fext + dampStr * M * V_pred,
    // with A_pred, V_pred consistent with Newmark from (X^n, V^n, A^n).
    // Also accumulates the tangent K = d/dX (M/(beta dt^2) + grad V).
    void assembleResidualAndTangent(const Eigen::MatrixXd& Xtrial,
                                    double dt,
                                    const Eigen::MatrixXd& Xn,
                                    const Eigen::MatrixXd& Vn,
                                    const Eigen::MatrixXd& An,
                                    Eigen::VectorXd& res,
                                    Eigen::SparseMatrix<double>& K) const;
    // Eliminate frozen DOFs by zeroing rows/cols + diagonal placeholder.
    void applyClamps(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& res) const;
};

} // namespace ns

#endif
