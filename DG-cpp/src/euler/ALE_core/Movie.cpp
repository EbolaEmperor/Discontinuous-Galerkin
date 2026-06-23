#include "Movie.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace euler_ale {
namespace {

int frameIndexFromPath(const std::filesystem::path& path) {
    std::string stem = path.stem().string();
    const std::string prefix = "frame_";
    if (stem.rfind(prefix, 0) != 0) return -1;
    try {
        return std::stoi(stem.substr(prefix.size()));
    } catch (const std::exception&) {
        return -1;
    }
}

} // namespace

int runOutputCommand(const std::string& command) {
    std::cout << command << "\n";
    int rc = std::system(command.c_str());
    if (rc != 0) {
        std::cerr << "Warning: output command failed with code " << rc << "\n";
    }
    return rc;
}

void clearFrameDirectory(const std::string& dir) {
    namespace fs = std::filesystem;
    if (!fs::exists(dir)) return;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.path().extension() == ".ppm") fs::remove(e.path());
    }
}

void pruneFramesFrom(const std::string& dir, int firstInvalidFrame) {
    namespace fs = std::filesystem;
    if (!fs::exists(dir)) return;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.path().extension() != ".ppm") continue;
        int idx = frameIndexFromPath(e.path());
        if (idx >= firstInvalidFrame) fs::remove(e.path());
    }
}

void trimDiagnosticsToTime(const std::string& csvPath, double time) {
    namespace fs = std::filesystem;
    if (!fs::exists(csvPath)) return;
    std::ifstream in(csvPath);
    if (!in) return;
    fs::path tmp = fs::path(csvPath);
    tmp += ".tmp";
    std::ofstream out(tmp, std::ios::trunc);
    if (!out) return;

    std::string line;
    if (std::getline(in, line)) out << line << "\n";
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        size_t comma = line.find(',');
        std::string first = (comma == std::string::npos) ? line : line.substr(0, comma);
        try {
            double rowTime = std::stod(first);
            if (rowTime <= time + 1e-10) out << line << "\n";
        } catch (const std::exception&) {
        }
    }
    out.close();
    if (!out) return;
    std::error_code ec;
    fs::rename(tmp, csvPath, ec);
    if (ec) {
        fs::remove(csvPath, ec);
        ec.clear();
        fs::rename(tmp, csvPath, ec);
    }
}

void MovieWriter::setup() {
    namespace fs = std::filesystem;
    fs::create_directories("out");
    for (auto& ch : channels) {
        fs::create_directories(ch.dir);
        clearFrameDirectory(ch.dir);
    }
}

void MovieWriter::pruneFrom(int frame) {
    for (auto& ch : channels) {
        pruneFramesFrom(ch.dir, frame);
    }
}

std::string MovieWriter::framePath(int channel, int frame) const {
    char fn[512];
    std::snprintf(fn, sizeof(fn), "%s/frame_%05d.ppm",
                  channels[channel].dir.c_str(), frame);
    return std::string(fn);
}

void MovieWriter::finalize(int lastFrame) {
    namespace fs = std::filesystem;
    for (auto& ch : channels) {
        std::string lastFramePath = framePath(
            static_cast<int>(&ch - channels.data()), lastFrame);
        if (fs::exists(lastFramePath)) {
            std::error_code ec;
            fs::copy_file(lastFramePath, ch.still,
                          fs::copy_options::overwrite_existing, ec);
        }
        if (!ch.stillPng.empty()) {
            runOutputCommand("ffmpeg -y -i " + ch.still +
                             " -frames:v 1 -update 1 " + ch.stillPng);
        }
        runOutputCommand("ffmpeg -y -framerate " + std::to_string(fps) +
                         " -i " + ch.dir + "/frame_%05d.ppm"
                         " -c:v libx264 -pix_fmt yuv420p -crf " +
                         std::to_string(ch.crf) + " " + ch.video);
    }
}

} // namespace euler_ale
