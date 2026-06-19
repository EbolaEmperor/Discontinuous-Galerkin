#ifndef BLAST_ROD_DOMAIN_H
#define BLAST_ROD_DOMAIN_H

#include "Geometry.h"
#include "ElasticSolid.h"

namespace euler_ale {

double blastRodFluidSdf(const BlastRodGeom& geom, double x, double y);
Mesh makeBlastRodMesh(const BlastRodGeom& geom, double h, int maxIter, bool verbose);
Mesh makeCurrentSolidBlastRodMesh(const BlastRodGeom& geom, const ElasticSolid2D& solid,
                                  double h, int maxIter, bool verbose);
int blastRodBoundaryTag(double x, double y, double t, const BlastRodMap& map);
int blastRodBoundaryTag(double x, double y, double t, const SolidALEMap& map,
                        const BlastRodGeom& geom);

} // namespace euler_ale

#endif
