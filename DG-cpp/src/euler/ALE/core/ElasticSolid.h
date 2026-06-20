#ifndef ELASTIC_SOLID_H
#define ELASTIC_SOLID_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <memory>
#include <string>
#include <vector>

namespace euler_ale {

using namespace Eigen;

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
    double velocityLimit;
    double displacementLimit;

    SolidMaterial();
};

struct SolidMeshQuality {
    double minAngleDeg = 0.0;
    double meanAngleDeg = 0.0;
    double minArea = 0.0;
    double maxArea = 0.0;
    double minEdge = 0.0;
};

class ElasticSolid2D {
public:
    ElasticSolid2D();

    void buildRectangularBeam(double x0, double x1, double y0, double y1,
                              int nx, int ny, const SolidMaterial& material);
    void buildRoundedRootBeam(double x0, double x1, double y0, double y1,
                              double rootRadius, int nx, int ny,
                              const SolidMaterial& material);
    void clearExternalForces();
    void addBoundaryTractionAt(const Vector2d& point, const Vector2d& traction, double measure);
    void advanceExplicit(double dt);
    void setState(const MatrixXd& nodes, const MatrixXd& velocities);

    int numNodes() const;
    int numElements() const;
    double totalMass() const;
    double stableTimeStep(double cfl) const;
    double maxNodeSpeed() const;
    double tipDisplacementX() const;
    double tipVelocityX() const;
    SolidMeshQuality meshQuality() const;

    const MatrixXd& referenceNodes() const;
    const MatrixXd& currentNodes() const;
    const MatrixXd& velocities() const;
    const MatrixXi& elements() const;
    const std::vector<SolidBoundarySegment>& boundarySegments() const;
    const std::vector<SolidBoundarySegment>& movingBoundarySegments() const;

private:
    int nodeIndex(int i, int j) const;
    void assembleStiffness();
    void appendBoundarySegments();
    int closestMovingSegment(const Vector2d& point, double& s, double& dist2) const;

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

class SolidALEMap {
public:
    SolidALEMap();
    ~SolidALEMap();

    SolidALEMap(const SolidALEMap&) = delete;
    SolidALEMap& operator=(const SolidALEMap&) = delete;
    SolidALEMap(SolidALEMap&&) noexcept;
    SolidALEMap& operator=(SolidALEMap&&) noexcept;

    void setSolid(const ElasticSolid2D* solid);
    void setDomain(double xa, double xb, double ya, double yb);
    void setInfluence(double normalDecay, double wallMargin);
    void setReferenceNodes(const MatrixXd& referenceNodes);
    void setCurrent(double time, const MatrixXd& nodes, const MatrixXd& velocities);
    void setMotion(double t0, double t1, const MatrixXd& nodes0, const MatrixXd& velocities0,
                   const MatrixXd& nodes1, const MatrixXd& velocities1);

    Vector2d refToPhys(const Vector2d& X, double time) const;
    Vector2d velocityAt(double x, double y, double time) const;
    Vector2d velocityAtCached(double x, double y, double time, int& segmentHint) const;
    double maxMeshSpeed(double time) const;
    double distanceToBoundary(double x, double y, double time) const;

private:
    struct SegmentTree;

    double alpha(double time) const;
    double influence(double x, double y, double distance) const;
    Vector2d nodeAt(int i, double a) const;
    Vector2d nodeStepVelocity(int i) const;
    void rebuildReferenceTree();
    void rebuildMotionTree();
    int closestReferenceSegment(const Vector2d& point, double& s, double& dist2) const;
    int closestCurrentSegment(const Vector2d& point, double a, double& s, double& dist2) const;
    int closestCurrentSegment(const Vector2d& point, double a, int segmentHint,
                              double& s, double& dist2) const;

    const ElasticSolid2D* solid_;
    double xa_;
    double xb_;
    double ya_;
    double yb_;
    double normalDecay_;
    double wallMargin_;
    double t0_;
    double t1_;
    MatrixXd nodes0_;
    MatrixXd nodes1_;
    MatrixXd velocities0_;
    MatrixXd velocities1_;
    MatrixXd referenceNodes_;
    std::unique_ptr<SegmentTree> referenceTree_;
    std::unique_ptr<SegmentTree> currentTree_;
};

void overlaySolidMesh(const std::string& path, const ElasticSolid2D& solid,
                      double xa, double xb, double ya, double yb, bool drawInterior);
void overlaySolidMesh(std::vector<unsigned char>& image, int width, int height,
                      const ElasticSolid2D& solid,
                      double xa, double xb, double ya, double yb,
                      bool drawInterior);

} // namespace euler_ale

#endif
