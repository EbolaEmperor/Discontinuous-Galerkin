#ifndef BLAST_ROD_CASE_H
#define BLAST_ROD_CASE_H

#include "Checkpoint.h"
#include "Coupling.h"
#include "DGState.h"
#include "Domain.h"
#include "Geometry.h"
#include "Output.h"

#include <string>

namespace euler_ale {

int runBlastRod(bool quick, bool freshStart = false,
                const std::string& configPath = std::string(),
                int frameOverride = 0, double tEndOverride = 0.0);

} // namespace euler_ale

#endif
