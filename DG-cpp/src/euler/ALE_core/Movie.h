#ifndef EULER_ALE_MOVIE_H
#define EULER_ALE_MOVIE_H

#include <string>
#include <vector>

namespace euler_ale {

struct MovieChannel {
    std::string dir;
    std::string video;
    std::string still;
    std::string stillPng;
    int crf = 16;
};

struct MovieWriter {
    std::string prefix;
    int fps = 30;
    std::vector<MovieChannel> channels;

    void setup();
    void pruneFrom(int frame);
    std::string framePath(int channel, int frame) const;
    void finalize(int lastFrame);
};

int runOutputCommand(const std::string& command);
void clearFrameDirectory(const std::string& dir);
void pruneFramesFrom(const std::string& dir, int firstInvalidFrame);
void trimDiagnosticsToTime(const std::string& csvPath, double time);

} // namespace euler_ale

#endif
