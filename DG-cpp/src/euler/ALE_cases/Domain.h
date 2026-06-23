#ifndef BLAST_ROD_DOMAIN_H
#define BLAST_ROD_DOMAIN_H

#include "Geometry.h"
#include "ElasticSolid.h"

namespace euler_ale {

struct BlastRodMeshOptions {
    double hNearFactor = 0.25;
    double gradeRadius = 0.36;
    double seedOffsetX = 0.0;
    double seedOffsetY = 0.0;
    unsigned randomSeed = 12345u;
};

double blastRodFluidSdf(const BlastRodGeom& geom, double x, double y);
Mesh makeBlastRodMesh(const BlastRodGeom& geom, double h, int maxIter, bool verbose,
                      const BlastRodMeshOptions& options = {});
Mesh makeBlastRodSolidReferenceMesh(const BlastRodGeom& geom, double h,
                                    int maxIter, bool verbose);
Mesh makeCurrentBlastRodSolidInteriorMesh(const ElasticSolid2D& solid, double h,
                                          int maxIter, bool verbose);
Mesh makeCurrentSolidBlastRodMesh(const BlastRodGeom& geom, const ElasticSolid2D& solid,
                                  double h, int maxIter, bool verbose,
                                  const BlastRodMeshOptions& options = {});
int blastRodBoundaryTag(double x, double y, double t, const BlastRodMap& map);
int blastRodBoundaryTag(double x, double y, double t, const SolidALEMap& map,
                        const BlastRodGeom& geom);

} // namespace euler_ale

#endif
