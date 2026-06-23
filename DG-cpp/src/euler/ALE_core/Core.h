#ifndef ALE_CORE_H
#define ALE_CORE_H

#include <Eigen/Dense>

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "DG.h"
#include "FEM.h"
#include "Mesh.h"

namespace euler_ale {

using namespace Eigen;

extern const int TAG_INTERIOR;
extern const int TAG_MOVING_WALL;
extern const int TAG_OUTFLOW;
extern const int TAG_SLIP_WALL;
extern const int TAG_EXACT;

using RefMapFn = std::function<Vector2d(const Vector2d&, double)>;
using MeshVelocityFn = std::function<Vector2d(double, double, double)>;
using MaxMeshSpeedFn = std::function<double(double)>;
using ALEBCFn = std::function<Vector4d(double, double, double, const Vector4d&, double, double, int, double)>;
using Tagger = std::function<int(double, double, double)>;

struct EdgeOnElem;
struct Space;
class SolidALEMap;

class ALEAdaptiveForest {
public:
    ALEAdaptiveForest(const Mesh& baseMesh, int order, int nComp);
    ~ALEAdaptiveForest();

    ALEAdaptiveForest(const ALEAdaptiveForest&) = delete;
    ALEAdaptiveForest& operator=(const ALEAdaptiveForest&) = delete;
    ALEAdaptiveForest(ALEAdaptiveForest&&) noexcept;
    ALEAdaptiveForest& operator=(ALEAdaptiveForest&&) noexcept;

    int numLeaves() const;
    int gen(int flatElem) const;
    void syncFromState(const MatrixXd& U, const MatrixXi& e2d, int locDof);
    MatrixXd gatherState(const MatrixXi& e2d, int locDof, int nDof) const;
    void buildMesh(Mesh& mesh, const RefMapFn& map, double time);
    std::pair<int, int> adapt(const std::vector<int>& flag, int maxGen);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

struct EdgeOnElem {
    int et = 0;
    int dir = 0;
    Vector2d nout = Vector2d::Zero();
    double he = 0.0;
};

struct Space {
    Mesh mesh;
    std::unique_ptr<FEM> fem;
    MatrixXi e2d;
    MatrixXi edge;
    MatrixXi e2s;
    VectorXi tag;
    int nDof = 0;
    double minH = 0.0;
    std::vector<int> volumeVelocityHints;
    std::vector<int> edgeVelocityHints;
};

struct StaticSpace {
    Mesh referenceMesh;
    Space space;
};

struct ModalState {
    double q = 0.0;
    double qd = 0.0;
};

double triArea(const Mesh& mesh, int elem);
double longestEdge(const Mesh& mesh, int elem);
double hCFL(const Mesh& mesh, int elem);
EdgeOnElem edgeOnElem(const Mesh& mesh, int elem, int n1, int n2);

Vector4d aleRusanov(const Vector4d& UL, const Vector4d& UR, double nx, double ny, double wn);
Vector4d movingWallGhost(const Vector4d& Um, double nx, double ny, double wn);
Vector4d characteristicPressureOutletGhost(const Vector4d& Um, double nx, double ny,
                                           double wn, const Vector4d& farPrimitive);

void rebuildSpace(ALEAdaptiveForest& forest, int order, const RefMapFn& refMap,
                  double time, const Tagger& tagger, Space& space);
void initializeStaticSpace(ALEAdaptiveForest& forest, int order, const RefMapFn& refMap,
                           double time, const Tagger& tagger, StaticSpace& staticSpace);
void updateStaticSpace(StaticSpace& staticSpace, const RefMapFn& refMap, double time);
MatrixXd advanceOne(ALEAdaptiveForest& forest, int order, const RefMapFn& refMap,
                    double time, double dt, const Tagger& tagger,
                    const MeshVelocityFn& meshVelocity, const ALEBCFn& bc,
                    const MatrixXd& U);
MatrixXd advanceOneStatic(StaticSpace& spaceNow, StaticSpace& spaceNext,
                          const RefMapFn& refMap, double time, double dt,
                          const MeshVelocityFn& meshVelocity, const ALEBCFn& bc,
                          const MatrixXd& U);
MatrixXd advanceOneStatic(StaticSpace& spaceNow, StaticSpace& spaceNext,
                          const RefMapFn& refMap, double time, double dt,
                          const SolidALEMap& meshMap, const ALEBCFn& bc,
                          const MatrixXd& U);
double estimateDt(const Space& space, const MatrixXd& U, const MaxMeshSpeedFn& maxMeshSpeed,
                  double time, int order, double cfl);

std::vector<int> computeAMRFlags(const Space& space, const MatrixXd& U,
                                 const ALEAdaptiveForest& forest, int maxGen,
                                 double refineThreshold, double coarsenThreshold,
                                 int bufferLayers, bool allowCoarsen);

bool readPPM(const std::string& path, int& width, int& height,
             std::vector<unsigned char>& image);
void writePPM(const std::string& path, int width, int height,
              const std::vector<unsigned char>& image);
bool downsamplePPM(const std::string& srcPath, const std::string& dstPath,
                   int factor);
std::vector<unsigned char> downsampleImage(const std::vector<unsigned char>& src,
                                           int srcWidth, int srcHeight,
                                           int factor,
                                           int& dstWidth, int& dstHeight);
void drawLine(std::vector<unsigned char>& image, int width, int height,
              int x0, int y0, int x1, int y1);
void overlayMesh(std::vector<unsigned char>& image, int width, int height,
                 const Mesh& mesh, double xa, double xb, double ya, double yb);
void overlayMesh(const std::string& path, const Mesh& mesh,
                 double xa, double xb, double ya, double yb);

} // namespace euler_ale

#endif
