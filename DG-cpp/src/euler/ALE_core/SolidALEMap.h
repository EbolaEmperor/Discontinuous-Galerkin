#ifndef ALE_SOLID_ALE_MAP_H
#define ALE_SOLID_ALE_MAP_H

#include "ElasticSolid.h"

#include <memory>
#include <string>
#include <vector>

namespace euler_ale {

using namespace Eigen;

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
    const MatrixXd& referenceNodes() const;

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
