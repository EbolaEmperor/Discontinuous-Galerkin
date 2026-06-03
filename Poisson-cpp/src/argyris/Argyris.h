#ifndef ARGYRIS_H
#define ARGYRIS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "Mesh.h"
#include "ExactSolution.h"

// ---------------------------------------------------------------------------
// Argyris element (C^1-conforming quintic triangle) for the biharmonic problem
//
//     Delta^2 u = f   in Omega,    u = g_1, d_n u = g_2 on dOmega  (clamped).
//
// The local space is P_5 (21-dimensional).  The 21 degrees of freedom are
//   * at each of the 3 vertices:  u, u_x, u_y, u_xx, u_xy, u_yy   (18 total)
//   * at each of the 3 edge midpoints:  the normal derivative d_n u   (3 total)
// which makes the global finite-element space a subspace of H^2 (C^1).
//
// We build the nodal basis per element by inverting the 21x21 matrix of DOF
// functionals applied to the raw degree-5 barycentric basis.  Vertex DOFs are
// expressed in the global Cartesian frame (so they are shared directly between
// elements) and each edge normal-derivative DOF uses a *globally fixed* edge
// normal, so the two elements sharing a face produce the identical functional
// and the DOF is shared without sign bookkeeping.  Together this yields a
// genuinely C^1 conforming discretization (no interior penalty needed).
//
// Weak form (conforming): find u_h with the clamped data imposed essentially,
//   sum_K \int_K D^2 u_h : D^2 v_h = \int_Omega f v_h     for all v_h.
// ---------------------------------------------------------------------------

struct ArgyrisResult {
    Eigen::VectorXd c;     // global dof vector
    double errL2 = 0.0;
    double errH1 = 0.0;    // full H1 norm
    double errH2 = 0.0;    // H2 seminorm |D^2(u-u_h)|
    int    nDof  = 0;
    int    nFree = 0;
    double tAssemble = 0.0;
    double tSolve    = 0.0;
    double tError    = 0.0;
    std::string solverName;
};

class ArgyrisFEM {
public:
    static constexpr int locDof = 21;

    Mesh& mesh;
    std::vector<Eigen::MatrixXd> Dlam;   // 2x3 per element
    Eigen::VectorXd area;                // NT
    std::vector<Eigen::MatrixXd> Rplain; // 3x6 per element: [Hxx;Hyy;Hxy] <- barycentric Hess
    std::vector<Eigen::MatrixXd> coef;   // 21x21 per element (columns = nodal basis)

    Eigen::MatrixXi edge, edge2side;     // skeleton (sorted endpoints, adjacency)
    Eigen::MatrixXi elem2dof;            // NT x 21 global dof indices
    int nDof = 0;
    int nNode = 0, nEdge = 0;

    explicit ArgyrisFEM(Mesh& m);

    // Assemble the conforming biharmonic stiffness  \int_K D^2 phi_i : D^2 phi_j.
    Eigen::SparseMatrix<double> assembleStiffness() const;

    // Load vector  \int_K f phi_i  with f = Delta^2 u_exact.
    Eigen::VectorXd assembleLoad(const ExactSolution& sol) const;

    // Clamped boundary data imposed strongly from the exact solution:
    // all 6 DOFs at boundary vertices + the normal-derivative DOF on boundary
    // edges are pinned to the nodal interpolant of u_exact.  Returns the lifted
    // vector c (boundary entries set) and the list of free dofs.
    void strongBC(const ExactSolution& sol, Eigen::VectorXd& c,
                  std::vector<int>& freeDof) const;

    // L2 error, full H1 error, and H2 seminorm error.
    void errors(const Eigen::VectorXd& c, const ExactSolution& sol,
                double& errL2, double& errH1, double& errH2) const;

private:
    void buildLocalBases();
    void buildDOFMap();
};

#endif
