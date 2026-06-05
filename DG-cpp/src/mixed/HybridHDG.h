#ifndef HYBRID_HDG_H
#define HYBRID_HDG_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "VecFEM.h"
#include "FEM.h"
#include "Mesh.h"
#include "ExactSolution.h"

// ---------------------------------------------------------------------------
// Hybridizable Discontinuous Galerkin (HDG) solver for the Poisson problem
//
//     sigma = grad(u),   -div(sigma) = f   in Omega,   u = g on dOmega.
//
// Internally we use the standard flux variable  q = -grad(u)  (so that
// div(q) = f) which yields a SYMMETRIC local solver and an SPD globally
// hybridized system.  The numerical trace  lambda = \hat{u}  lives on the
// mesh skeleton in P_k(e).
//
// Local solver on each element K (q in [P_k]^2, u in P_k):
//   ( q, v )_K - ( u, div v )_K + < lambda, v.n >_dK            = 0          (1)
//  -( q, grad w )_K + < q.n + tau (u - lambda), w >_dK          = (f, w)_K   (2)
// Conservativity / transmission (global equation for lambda, mu in P_k(e)):
//   sum_K < q.n + tau (u - lambda), mu >_dK = 0   on interior faces,
//   lambda = P_h g                                 on boundary faces.
//
// Static condensation of (q,u) per element gives the SPD Schur system
//   A * Lambda = b,     A = sum_K ( M_ll - S * Aloc^{-1} * R ),
// which is solved with CHOLMOD.  (q,u) are then recovered element-by-element.
//
// On return  sigmah = grad(u) = -q  is stored in the VecFEM dof layout and
// uh in the FEM dof layout, so the existing l2Error* / postprocessing
// routines can be reused verbatim.
// ---------------------------------------------------------------------------

struct HybridHDGResult {
    Eigen::VectorXd sigmah;   // sigma = grad u (= -q), VecFEM dof layout
    Eigen::VectorXd uh;       // u, FEM dof layout
    Eigen::VectorXd lambda;   // numerical trace, nLambda = NE*(k+1)
    int    nLambda    = 0;    // total trace dofs
    int    nLambdaFree= 0;    // interior (solved) trace dofs
    double tAssemble  = 0.0;   // local condensation -> global Schur assembly
    double tSolve     = 0.0;   // factor + solve of the SPD trace system
    double tRecover   = 0.0;   // local back-substitution of (q,u)
    std::string solverName;
};

HybridHDGResult solveHybridHDG(VecFEM& femQ, FEM& femU, Mesh& mesh,
                               const Eigen::MatrixXi& elem2dofSigma,
                               const Eigen::MatrixXi& elem2dofU,
                               const Eigen::MatrixXi& edge,
                               const Eigen::MatrixXi& edge2side,
                               double tau, const ExactSolution& sol,
                               int solver_type);

#endif
