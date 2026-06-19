#ifndef ALE_DG_STATE_H
#define ALE_DG_STATE_H

#include "Core.h"

namespace euler_ale {

void applyPrimitiveBounds(MatrixXd& U, double rhoFloor, double pFloor, double speedMax);
MatrixXd interpolateDGToSpace(Space& oldSpace, const MatrixXd& oldU, Space& newSpace);

} // namespace euler_ale

#endif
