#ifndef IB_COUPLER_H
#define IB_COUPLER_H

// ===========================================================================
// Immersed-boundary transfer operators between the DG (dP_k) flow space and
// a 1-D Lagrangian filament (CosseratFilament).
//
//   * meshToRod : sample (u, v) at every rod node X_k via the DG basis
//                 (the same routine used by the streakline tracer).
//   * rodToMesh : scatter a per-node force F_k onto the DG nodal-load vector
//                 (one Eigen vector per component).  Each F_k is multiplied
//                 by the dual-cell arclength weight  Δs_k  before scattering,
//                 so the discrete pairing is
//                     <F, I(u)>_{rod}  =  sum_k F_k · u_h(X_k) · Δs_k
//                 and the scatter operator S = I^T  satisfies the adjoint
//                 identity  <F, I u>_rod = <S F, u>_DG  exactly.
//
// This is the consistent / FE-immersed-boundary form: the DG basis itself is
// the transfer kernel (Heltai & Costanzo 2012; Boffi & Gastaldi 2003).  No
// regularised δ function -> no smearing knob, exact discrete adjointness.
// ===========================================================================

#include <Eigen/Dense>
#include "FEM.h"
#include "Mesh.h"
#include "NavierStokes.h"          // MeshLocator
#include "CosseratFilament.h"

namespace ns {

// Compute the dual-cell arclength weights for the rod nodes (size N+1).
// Trapezoidal: ds_0 = l_0/2,  ds_N = l_{N-1}/2,  ds_i = (l_{i-1}+l_i)/2.
// Uses the *current* edge lengths from rod.X (not the rest lengths) so the
// integrals stay accurate when the rod stretches.
Eigen::VectorXd rodArclengthWeights(const CosseratFilament& rod);

// Sample (u, v) at every rod node X_k.  `hint` (size N+1) caches the last
// containing element per node and is updated in place.  On nodes that fall
// outside the mesh the returned velocity is zero (alive[k] = false).
struct RodVelocitySample {
    Eigen::MatrixXd uv;       // (N+1) x 2
    std::vector<char> alive;  // size N+1
};
RodVelocitySample meshToRod(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
                            const Eigen::VectorXd& uField, const Eigen::VectorXd& vField,
                            const MeshLocator& loc,
                            const CosseratFilament& rod,
                            std::vector<int>& hint);

// Scatter per-node forces (size N+1, two components) onto the DG nodal-load
// vector (size nDof, two components).  Each F_k is multiplied by Δs_k before
// scatter.  loadX, loadY are *added to*, so the caller can pre-zero or
// accumulate as wanted.  Nodes outside the mesh are silently skipped.
void rodToMesh(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
               const MeshLocator& loc,
               const CosseratFilament& rod,
               const Eigen::MatrixXd& F,             // (N+1) x 2
               const Eigen::VectorXd& dsWeights,     // (N+1)
               std::vector<int>& hint,
               Eigen::VectorXd& loadX, Eigen::VectorXd& loadY);

// Overlay a polyline (white, anti-aliased) onto a binary PPM file in place.
// The polyline is given in physical coords; the file is treated as the image
// of the world rectangle [xmin, xmax] x [ymin, ymax] (note: pixel row 0 = top).
// Used by the FSI driver to draw the filament on every rendered frame.
void overlayPolylineOnPPM(const std::string& path,
                          const Eigen::MatrixXd& X,        // (N+1) x 2 polyline nodes
                          double xmin, double xmax, double ymin, double ymax,
                          int width = 2, int rgb = 0xFFFFFF);

} // namespace ns

#endif
