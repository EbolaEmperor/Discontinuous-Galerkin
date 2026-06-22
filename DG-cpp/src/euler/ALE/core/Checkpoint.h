#ifndef ALE_CHECKPOINT_H
#define ALE_CHECKPOINT_H

#include "Core.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace euler_ale {

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
    Mesh solidReferenceMesh;
    MatrixXd U;
    MatrixXd solidNodes;
    MatrixXd solidVelocities;
    MatrixXd aleReferenceNodes;
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

} // namespace euler_ale

#endif
