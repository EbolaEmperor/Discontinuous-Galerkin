#ifndef ALE_BODY_FITTED_MESH_H
#define ALE_BODY_FITTED_MESH_H

#include "ElasticSolid.h"
#include "MeshGen.h"

#include <functional>
#include <utility>
#include <vector>

namespace euler_ale {

struct SolidBodyMeshSpec {
    double xa = 0.0;
    double xb = 1.0;
    double ya = 0.0;
    double yb = 1.0;
    double hFar = 0.02;
    double hNearFactor = 0.25;
    double gradeRadius = 0.36;
    int maxIter = 160;
    bool verbose = false;
    std::function<double(double, double)> solidDistance;
    std::vector<std::pair<Vector2d, Vector2d>> fixedSegments;
    const ElasticSolid2D* solid = nullptr;
    bool addMovingSolidBoundary = false;
};

double signedBox(double x, double y, double x0, double x1, double y0, double y1);
double segmentProjection(const Vector2d& p, const Vector2d& a, const Vector2d& b,
                         double& dist2);
double signedDistancePolygon(const std::vector<Vector2d>& poly, double x, double y);
std::vector<Vector2d> solidBoundaryLoop(const ElasticSolid2D& solid);
void addGradedFixedSegment(DistanceMeshSpec& spec, const Vector2d& a, const Vector2d& b,
                           const std::function<double(double, double)>& sizeFn,
                           double hMin);
Mesh makeSolidBodyFittedMesh(const SolidBodyMeshSpec& spec);

} // namespace euler_ale

#endif
