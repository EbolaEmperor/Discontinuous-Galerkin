#include "Checkpoint.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace euler_ale {
namespace {

void writeMatrixXd(std::ostream& out, const std::string& name, const MatrixXd& m) {
    out << name << " " << m.rows() << " " << m.cols() << "\n";
    out << std::setprecision(17);
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            if (j) out << " ";
            out << m(i, j);
        }
        out << "\n";
    }
}

void writeMatrixXi(std::ostream& out, const std::string& name, const MatrixXi& m) {
    out << name << " " << m.rows() << " " << m.cols() << "\n";
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            if (j) out << " ";
            out << m(i, j);
        }
        out << "\n";
    }
}

bool readMatrixXd(std::istream& in, const std::string& name, MatrixXd& m) {
    std::string got;
    int rows = 0;
    int cols = 0;
    if (!(in >> got >> rows >> cols) || got != name || rows < 0 || cols < 0) return false;
    m.resize(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (!(in >> m(i, j))) return false;
    return true;
}

bool readMatrixXi(std::istream& in, const std::string& name, MatrixXi& m) {
    std::string got;
    int rows = 0;
    int cols = 0;
    if (!(in >> got >> rows >> cols) || got != name || rows < 0 || cols < 0) return false;
    m.resize(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (!(in >> m(i, j))) return false;
    return true;
}

bool loadCheckpoint(const std::filesystem::path& path, RunCheckpoint& cp) {
    std::ifstream in(path);
    if (!in) return false;
    std::string magic;
    int version = 0;
    if (!(in >> magic >> version) || magic != "ALE_RUN_CHECKPOINT" || version != 1)
        return false;

    auto readScalar = [&](const std::string& name, auto& value) {
        std::string got;
        if (!(in >> got >> value) || got != name) return false;
        return true;
    };

    int quickFlag = 0;
    if (!readScalar("quick", quickFlag)) return false;
    cp.quick = quickFlag != 0;
    if (!readScalar("ord", cp.ord)) return false;
    if (!readScalar("n_frames", cp.nFrames)) return false;
    if (!readScalar("t_end", cp.tEnd)) return false;
    if (!readScalar("h", cp.h)) return false;
    if (!readScalar("time", cp.time)) return false;
    if (!readScalar("next_frame", cp.nextFrame)) return false;
    if (!readScalar("step", cp.step)) return false;
    if (!readScalar("frame", cp.frame)) return false;
    if (!readScalar("remesh_count", cp.remeshCount)) return false;

    std::string got;
    if (!(in >> got >> cp.milestoneDone[0] >> cp.milestoneDone[1] >>
          cp.milestoneDone[2]) || got != "milestone_done") {
        return false;
    }

    if (!readMatrixXd(in, "mesh_node", cp.referenceMesh.node)) return false;
    if (!readMatrixXi(in, "mesh_elem", cp.referenceMesh.elem)) return false;
    if (!readMatrixXd(in, "solid_nodes", cp.solidNodes)) return false;
    if (!readMatrixXd(in, "solid_velocities", cp.solidVelocities)) return false;
    if (!readMatrixXd(in, "ale_reference_nodes", cp.aleReferenceNodes)) return false;
    if (!readMatrixXd(in, "dg_state", cp.U)) return false;
    if (!(in >> got) || got != "END") return false;
    return true;
}

bool compatibleCheckpoint(const RunCheckpoint& cp, bool quick, int ord,
                          int nFrames, double tEnd, double h, int solidNodes) {
    if (cp.quick != quick || cp.ord != ord || cp.nFrames != nFrames) return false;
    if (std::abs(cp.tEnd - tEnd) > 1e-12 || std::abs(cp.h - h) > 1e-12) return false;
    if (cp.referenceMesh.node.cols() != 2 || cp.referenceMesh.elem.cols() != 3) return false;
    if (cp.solidNodes.rows() != solidNodes || cp.solidNodes.cols() != 2) return false;
    if (cp.solidVelocities.rows() != solidNodes || cp.solidVelocities.cols() != 2) return false;
    if (cp.aleReferenceNodes.rows() != solidNodes || cp.aleReferenceNodes.cols() != 2) return false;
    if (cp.U.cols() != 4 || cp.time < -1e-14 || cp.time > tEnd + 1e-12) return false;
    return true;
}

} // namespace

const std::array<double, 3> kCheckpointFractions{{0.25, 0.50, 0.75}};
const std::array<const char*, 3> kCheckpointLabels{{"25", "50", "75"}};

std::filesystem::path checkpointPath(const std::string& casePrefix, bool quick,
                                     const std::string& label) {
    std::string prefix = quick ? casePrefix + "_quick_checkpoint_" :
                                 casePrefix + "_checkpoint_";
    return std::filesystem::path("out") / (prefix + label + ".chk");
}

bool writeCheckpointAtomic(const std::filesystem::path& path,
                           const RunCheckpoint& cp) {
    namespace fs = std::filesystem;
    fs::create_directories(path.parent_path());
    fs::path tmp = path;
    tmp += ".tmp";
    {
        std::ofstream out(tmp, std::ios::trunc);
        if (!out) {
            std::cerr << "Warning: cannot write checkpoint " << tmp << "\n";
            return false;
        }
        out << "ALE_RUN_CHECKPOINT 1\n";
        out << "quick " << (cp.quick ? 1 : 0) << "\n";
        out << "ord " << cp.ord << "\n";
        out << "n_frames " << cp.nFrames << "\n";
        out << std::setprecision(17);
        out << "t_end " << cp.tEnd << "\n";
        out << "h " << cp.h << "\n";
        out << "time " << cp.time << "\n";
        out << "next_frame " << cp.nextFrame << "\n";
        out << "step " << cp.step << "\n";
        out << "frame " << cp.frame << "\n";
        out << "remesh_count " << cp.remeshCount << "\n";
        out << "milestone_done " << cp.milestoneDone[0] << " "
            << cp.milestoneDone[1] << " " << cp.milestoneDone[2] << "\n";
        writeMatrixXd(out, "mesh_node", cp.referenceMesh.node);
        writeMatrixXi(out, "mesh_elem", cp.referenceMesh.elem);
        writeMatrixXd(out, "solid_nodes", cp.solidNodes);
        writeMatrixXd(out, "solid_velocities", cp.solidVelocities);
        writeMatrixXd(out, "ale_reference_nodes", cp.aleReferenceNodes);
        writeMatrixXd(out, "dg_state", cp.U);
        out << "END\n";
        if (!out) {
            std::cerr << "Warning: failed while writing checkpoint " << tmp << "\n";
            return false;
        }
    }

    std::error_code ec;
    fs::rename(tmp, path, ec);
    if (ec) {
        fs::remove(path, ec);
        ec.clear();
        fs::rename(tmp, path, ec);
    }
    if (ec) {
        std::cerr << "Warning: cannot finalize checkpoint " << path
                  << ": " << ec.message() << "\n";
        return false;
    }
    return true;
}

std::optional<RunCheckpoint> loadLatestCheckpoint(const std::string& casePrefix, bool quick,
                                                  int ord, int nFrames,
                                                  double tEnd, double h,
                                                  int solidNodes) {
    namespace fs = std::filesystem;
    for (int i = static_cast<int>(kCheckpointLabels.size()) - 1; i >= 0; --i) {
        fs::path path = checkpointPath(casePrefix, quick, kCheckpointLabels[i]);
        if (!fs::exists(path)) continue;
        RunCheckpoint cp;
        if (!loadCheckpoint(path, cp)) {
            std::cerr << "Warning: checkpoint unreadable, skipping " << path << "\n";
            continue;
        }
        if (!compatibleCheckpoint(cp, quick, ord, nFrames, tEnd, h, solidNodes)) {
            std::cerr << "Warning: checkpoint incompatible, skipping " << path << "\n";
            continue;
        }
        return cp;
    }
    return std::nullopt;
}

} // namespace euler_ale
