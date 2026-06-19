#ifndef ALE_CHECKPOINT_H
#define ALE_CHECKPOINT_H

#include "Core.h"

#include <array>
#include <filesystem>
#include <optional>
#include <string>

namespace euler_ale {

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
    std::array<int, 3> milestoneDone{{0, 0, 0}};
    Mesh referenceMesh;
    MatrixXd U;
    MatrixXd solidNodes;
    MatrixXd solidVelocities;
    MatrixXd aleReferenceNodes;
};

extern const std::array<double, 3> kCheckpointFractions;
extern const std::array<const char*, 3> kCheckpointLabels;

std::filesystem::path checkpointPath(const std::string& casePrefix, bool quick,
                                     const std::string& label);
bool writeCheckpointAtomic(const std::filesystem::path& path,
                           const RunCheckpoint& checkpoint);
std::optional<RunCheckpoint> loadLatestCheckpoint(const std::string& casePrefix, bool quick,
                                                  int ord, int nFrames,
                                                  double tEnd, double h,
                                                  int solidNodes);

} // namespace euler_ale

#endif
