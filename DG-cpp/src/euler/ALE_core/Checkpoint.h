#ifndef ALE_CHECKPOINT_H
#define ALE_CHECKPOINT_H

#include "Core.h"

#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace euler_ale {

class ElasticSolid2D;
class SolidALEMap;
struct SolidMaterial;

struct CheckpointMilestone {
    double fraction = 0.0;
    std::string label;
};

struct RunCheckpoint {
    bool quick = false;
    int ord = 0;
    int nFrames = 0;
    double tEnd = 0.0;
    double h = 0.0;
    double time = 0.0;
    double nextFrame = 0.0;
    int step = 0;
    int frame = 0;
    int remeshCount = 0;
    std::vector<int> milestoneDone;
    Mesh referenceMesh;
    Mesh currentMesh;
    Mesh solidReferenceMesh;
    MatrixXd U;
    MatrixXd solidNodes;
    MatrixXd solidVelocities;
    MatrixXd aleReferenceNodes;
};

using CheckpointSolidConfigureFn = std::function<void(ElasticSolid2D&)>;

struct CheckpointFlowRestoreResult {
    bool compatible = false;
    bool usedSavedPhysicalMesh = false;
    bool interpolated = false;
    bool resizedTopology = false;
    bool sameElementTopology = true;
    double meshNodeDiff = 0.0;
    int checkpointDof = 0;
    int targetDof = 0;
    int savedMeshDof = 0;
};

std::vector<CheckpointMilestone> checkpointSchedule(bool quick);
std::filesystem::path checkpointPath(const std::string& casePrefix, bool quick,
                                     const std::string& label);
bool writeCheckpointAtomic(const std::filesystem::path& path,
                           const RunCheckpoint& checkpoint);
void pruneOldCheckpoints(const std::string& casePrefix, bool quick, int keepCount);
std::optional<RunCheckpoint> loadLatestCheckpoint(const std::string& casePrefix, bool quick,
                                                  int ord, int nFrames,
                                                  double tEnd, double h,
                                                  int solidNodes,
                                                  bool allowExtension = false);
std::optional<RunCheckpoint> loadCheckpointByLabel(const std::string& casePrefix, bool quick,
                                                   const std::string& label,
                                                   int ord, int nFrames,
                                                   double tEnd, double h,
                                                   int solidNodes,
                                                   bool allowExtension = false);

double maxAbsMatrixDiff(const MatrixXd& a, const MatrixXd& b);
bool sameElementMatrix(const MatrixXi& a, const MatrixXi& b);
void buildSpaceOnMesh(const Mesh& mesh, int ord, Space& space);
bool restoreCheckpointSolidState(const RunCheckpoint& checkpoint,
                                 const SolidMaterial& material,
                                 ElasticSolid2D& solid,
                                 SolidALEMap& map,
                                 const CheckpointSolidConfigureFn& configureSolid = {});
bool restoreCheckpointAleMap(const RunCheckpoint& checkpoint,
                             const MatrixXd& aleReferenceNodes,
                             const ElasticSolid2D& solid,
                             SolidALEMap& map);
CheckpointFlowRestoreResult restoreCheckpointFlowState(const RunCheckpoint& checkpoint,
                                                       int ord,
                                                       Space& targetSpace,
                                                       MatrixXd& U);
RunCheckpoint makeRunCheckpoint(bool quick, int ord, int nFrames,
                                double tEnd, double h, double time,
                                double nextFrame, int step, int frame,
                                int remeshCount,
                                const std::vector<int>& milestoneDone,
                                const Mesh& referenceMesh,
                                const Space& currentSpace,
                                const ElasticSolid2D& solid,
                                const SolidALEMap& map,
                                const MatrixXd& U);

} // namespace euler_ale

#endif
