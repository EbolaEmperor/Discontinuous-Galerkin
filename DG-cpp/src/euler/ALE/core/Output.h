#ifndef ALE_OUTPUT_H
#define ALE_OUTPUT_H

#include <string>

namespace euler_ale {

int runOutputCommand(const std::string& command);
void clearFrameDirectory(const std::string& dir);
void pruneFramesFrom(const std::string& dir, int firstInvalidFrame);
void trimDiagnosticsToTime(const std::string& csvPath, double time);

} // namespace euler_ale

#endif
