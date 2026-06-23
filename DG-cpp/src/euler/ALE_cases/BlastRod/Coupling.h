#ifndef BLAST_ROD_COUPLING_H
#define BLAST_ROD_COUPLING_H

#include "Geometry.h"
#include "ElasticSolid.h"

namespace euler_ale {

ModalState advanceBlastRodSymplectic(ModalState state, double fluidForce,
                                     const BlastRodParams& params, double dt);
double blastRodGeneralizedForce(const Space& space, const MatrixXd& U, const BlastRodMap& map,
                                double time, double pExt, double* meanPressure = nullptr,
                                double* drag = nullptr);
double blastRodLoadSolidFromFluid(const Space& space, const MatrixXd& U, ElasticSolid2D& solid,
                                  double pExt, double* meanPressure = nullptr,
                                  double* drag = nullptr);

} // namespace euler_ale

#endif
