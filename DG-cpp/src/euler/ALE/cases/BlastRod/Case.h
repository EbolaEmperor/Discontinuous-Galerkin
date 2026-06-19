#ifndef BLAST_ROD_CASE_H
#define BLAST_ROD_CASE_H

#include "Checkpoint.h"
#include "Coupling.h"
#include "DGState.h"
#include "Domain.h"
#include "Geometry.h"
#include "Output.h"

namespace euler_ale {

int runBlastRod(bool quick, bool freshStart = false);

} // namespace euler_ale

#endif
