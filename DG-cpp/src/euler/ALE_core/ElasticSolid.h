#ifndef ELASTIC_SOLID_H
#define ELASTIC_SOLID_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Mesh.h"

#include <memory>
#include <string>
#include <vector>

namespace euler_ale {

using namespace Eigen;

struct Space;

enum SolidBoundarySide {
    SOLID_LEFT = 1,
    SOLID_RIGHT = 2,
    SOLID_TOP = 3,
    SOLID_BOTTOM = 4
};

struct SolidBoundarySegment {
    int a;
    int b;
    int side;
};

struct SolidMaterial {
    double density;
    double thickness;
    double young;
    double poisson;
    double damping;

    SolidMaterial();
};

struct SolidMeshQuality {
    double minAngleDeg = 0.0;
    double meanAngleDeg = 0.0;
    double minArea = 0.0;
    double maxArea = 0.0;
    double minEdge = 0.0;
    int invertedElements = 0;
};

class ElasticSolid2D {
public:
    ElasticSolid2D();

    void buildRectangularBeam(double x0, double x1, double y0, double y1,
                              int nx, int ny, const SolidMaterial& material);
    void buildRoundedRootBeam(double x0, double x1, double y0, double y1,
                              double rootRadius, int nx, int ny,
                              const SolidMaterial& material);
    void resetReferenceMesh(const Mesh& referenceMesh, const SolidMaterial& material);
    void remeshToReferenceMesh(const Mesh& referenceMesh);
    void remeshToCurrentMesh(const Mesh& currentMesh);
    void setFixedNodesInDisk(const Vector2d& center, double radius,
                             bool clearExisting = true);
    void setAllBoundarySegmentsMoving();
    void clearExternalForces();
    void addBoundaryTractionAt(const Vector2d& point, const Vector2d& traction, double measure);
    void smoothMovingBoundaryForces(int passes, double blend);
    void advanceExplicit(double dt);
    void setState(const MatrixXd& nodes, const MatrixXd& velocities);

    int numNodes() const;
    int numElements() const;
    double totalMass() const;
    double stableTimeStep(double cfl) const;
    double strainEnergy() const;
    double kineticEnergy() const;
    double maxNodeSpeed() const;
    double tipDisplacementX() const;
    double tipVelocityX() const;
    SolidMeshQuality meshQuality() const;
    SolidMeshQuality currentMeshQuality() const;

    const MatrixXd& referenceNodes() const;
    const MatrixXd& currentNodes() const;
    const MatrixXd& velocities() const;
    const MatrixXd& externalForces() const;
    const VectorXd& lumpedMasses() const;
    const VectorXi& fixedMask() const;
    const MatrixXi& elements() const;
    const std::vector<SolidBoundarySegment>& boundarySegments() const;
    const std::vector<SolidBoundarySegment>& movingBoundarySegments() const;

private:
    int nodeIndex(int i, int j) const;
    void assembleStiffness();
    void appendBoundarySegments();
    int closestMovingSegment(const Vector2d& point, double& s, double& dist2) const;
    void appendBoundarySegmentsFromMesh();

    int nx_;
    int ny_;
    double x0_;
    double x1_;
    double y0_;
    double y1_;
    SolidMaterial material_;
    MatrixXd X_;
    MatrixXd x_;
    MatrixXd v_;
    MatrixXd f_;
    MatrixXi elem_;
    VectorXd mass_;
    VectorXi fixed_;
    std::vector<int> tipNodes_;
    std::vector<SolidBoundarySegment> boundarySegments_;
    std::vector<SolidBoundarySegment> movingBoundarySegments_;
    SparseMatrix<double> stiffness_;
};

double loadSolidFromFluidPressure(const Space& space, const MatrixXd& U,
                                  ElasticSolid2D& solid, double pExt,
                                  double* meanPressure = nullptr,
                                  double* drag = nullptr,
                                  double* lift = nullptr);

} // namespace euler_ale

#endif
