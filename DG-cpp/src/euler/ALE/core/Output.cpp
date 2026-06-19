#include "Output.h"

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

} // namespace euler_ale
