#include "Checkpoint.h"

#include "DGState.h"
#include "ElasticSolid.h"
#include "SolidALEMap.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

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
    if (!(in >> magic >> version) || magic != "ALE_RUN_CHECKPOINT" ||
        (version != 1 && version != 2 && version != 3 && version != 4))
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
    if (!(in >> got) || got != "milestone_done") return false;
    if (version == 1) {
        cp.milestoneDone.assign(3, 0);
        if (!(in >> cp.milestoneDone[0] >> cp.milestoneDone[1] >>
              cp.milestoneDone[2])) {
            return false;
        }
    } else {
        int n = 0;
        if (!(in >> n) || n < 0 || n > 1000) return false;
        cp.milestoneDone.assign(n, 0);
        for (int i = 0; i < n; ++i)
            if (!(in >> cp.milestoneDone[i])) return false;
    }

    if (!readMatrixXd(in, "mesh_node", cp.referenceMesh.node)) return false;
    if (!readMatrixXi(in, "mesh_elem", cp.referenceMesh.elem)) return false;
    if (version >= 4) {
        if (!readMatrixXd(in, "current_mesh_node", cp.currentMesh.node)) return false;
        if (!readMatrixXi(in, "current_mesh_elem", cp.currentMesh.elem)) return false;
    } else {
        cp.currentMesh.node.resize(0, 2);
        cp.currentMesh.elem.resize(0, 3);
    }
    if (version >= 3) {
        if (!readMatrixXd(in, "solid_ref_node", cp.solidReferenceMesh.node)) return false;
        if (!readMatrixXi(in, "solid_ref_elem", cp.solidReferenceMesh.elem)) return false;
    } else {
        cp.solidReferenceMesh.node.resize(0, 2);
        cp.solidReferenceMesh.elem.resize(0, 3);
    }
    if (!readMatrixXd(in, "solid_nodes", cp.solidNodes)) return false;
    if (!readMatrixXd(in, "solid_velocities", cp.solidVelocities)) return false;
    if (!readMatrixXd(in, "ale_reference_nodes", cp.aleReferenceNodes)) return false;
    if (!readMatrixXd(in, "dg_state", cp.U)) return false;
    if (!(in >> got) || got != "END") return false;
    return true;
}

bool compatibleCheckpoint(const RunCheckpoint& cp, bool quick, int ord,
                          int nFrames, double tEnd, double h, int solidNodes,
                          bool allowExtension) {
    if (cp.quick != quick || cp.ord != ord) return false;
    if (std::abs(cp.h - h) > 1e-12) return false;
    if (allowExtension) {
        if (cp.nFrames > nFrames || cp.tEnd > tEnd + 1e-12) return false;
        if (!(cp.nFrames > 0) || !(nFrames > 0) || !(cp.tEnd > 0.0) ||
            !(tEnd > 0.0)) return false;
        double oldFrameDt = cp.tEnd / static_cast<double>(cp.nFrames);
        double newFrameDt = tEnd / static_cast<double>(nFrames);
        if (std::abs(oldFrameDt - newFrameDt) > 1e-10 * std::max(1.0, oldFrameDt))
            return false;
    } else {
        if (cp.nFrames != nFrames) return false;
        if (std::abs(cp.tEnd - tEnd) > 1e-12) return false;
    }
    if (cp.referenceMesh.node.cols() != 2 || cp.referenceMesh.elem.cols() != 3) return false;
    if (cp.currentMesh.node.rows() > 0) {
        if (cp.currentMesh.node.cols() != 2 || cp.currentMesh.elem.cols() != 3 ||
            cp.currentMesh.elem.rows() <= 0) {
            return false;
        }
    }
    int expectedSolidNodes = solidNodes;
    if (cp.solidReferenceMesh.node.rows() > 0) {
        if (cp.solidReferenceMesh.node.cols() != 2 ||
            cp.solidReferenceMesh.elem.cols() != 3 ||
            cp.solidReferenceMesh.elem.rows() <= 0) {
            return false;
        }
        expectedSolidNodes = static_cast<int>(cp.solidReferenceMesh.node.rows());
    }
    if (cp.solidNodes.rows() != expectedSolidNodes || cp.solidNodes.cols() != 2) return false;
    if (cp.solidVelocities.rows() != expectedSolidNodes ||
        cp.solidVelocities.cols() != 2) return false;
    if (cp.aleReferenceNodes.rows() != expectedSolidNodes ||
        cp.aleReferenceNodes.cols() != 2) return false;
    if (cp.U.cols() != 4 || cp.time < -1e-14 || cp.time > tEnd + 1e-12) return false;
    return true;
}

} // namespace

std::vector<CheckpointMilestone> checkpointSchedule(bool quick) {
    if (quick) return {{0.25, "25"}, {0.50, "50"}, {0.75, "75"}};
    std::vector<CheckpointMilestone> schedule;
    schedule.reserve(20);
    for (int pct = 5; pct <= 100; pct += 5) {
        schedule.push_back({0.01 * static_cast<double>(pct), std::to_string(pct)});
    }
    return schedule;
}

std::filesystem::path checkpointPath(const std::string& casePrefix, bool quick,
                                     const std::string& label) {
    std::string prefix = quick ? casePrefix + "_quick_checkpoint_" :
                                 casePrefix + "_checkpoint_";
    return std::filesystem::path("out") / (prefix + label + ".chk");
}

double maxAbsMatrixDiff(const MatrixXd& a, const MatrixXd& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return 1e300;
    double m = 0.0;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            m = std::max(m, std::abs(a(i, j) - b(i, j)));
    return m;
}

bool sameElementMatrix(const MatrixXi& a, const MatrixXi& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            if (a(i, j) != b(i, j)) return false;
    return true;
}

void buildSpaceOnMesh(const Mesh& mesh, int ord, Space& sp) {
    sp.mesh = mesh;
    sp.fem = std::make_unique<FEM>(ord, sp.mesh, false);
    sp.fem->getDOF(sp.mesh, sp.e2d, sp.nDof);
    sp.mesh.getEdge2Side(sp.edge, sp.e2s);
    sp.minH = 1e300;
    for (int k = 0; k < sp.mesh.elem.rows(); ++k) {
        sp.minH = std::min(sp.minH, hCFL(sp.mesh, k));
    }
    sp.tag = VectorXi::Zero(sp.edge.rows());
    sp.volumeVelocityHints.clear();
    sp.edgeVelocityHints.clear();
}

bool restoreCheckpointSolidState(const RunCheckpoint& cp,
                                 const SolidMaterial& material,
                                 ElasticSolid2D& solid,
                                 SolidALEMap& map,
                                 const CheckpointSolidConfigureFn& configureSolid) {
    if (cp.solidReferenceMesh.node.rows() > 0) {
        solid.resetReferenceMesh(cp.solidReferenceMesh, material);
    }
    if (configureSolid) configureSolid(solid);
    map.setSolid(&solid);
    if (cp.solidNodes.rows() != solid.numNodes() || cp.solidNodes.cols() != 2 ||
        cp.solidVelocities.rows() != solid.numNodes() || cp.solidVelocities.cols() != 2) {
        return false;
    }
    solid.setState(cp.solidNodes, cp.solidVelocities);
    return true;
}

bool restoreCheckpointAleMap(const RunCheckpoint& cp,
                             const MatrixXd& aleReferenceNodes,
                             const ElasticSolid2D& solid,
                             SolidALEMap& map) {
    if (aleReferenceNodes.rows() != solid.numNodes() ||
        aleReferenceNodes.cols() != 2) {
        return false;
    }
    map.setReferenceNodes(aleReferenceNodes);
    map.setCurrent(cp.time, solid.currentNodes(), solid.velocities());
    return true;
}

CheckpointFlowRestoreResult restoreCheckpointFlowState(const RunCheckpoint& cp,
                                                       int ord,
                                                       Space& targetSp,
                                                       MatrixXd& U) {
    CheckpointFlowRestoreResult result;
    result.checkpointDof = static_cast<int>(cp.U.rows());
    result.targetDof = targetSp.nDof;
    if (cp.U.cols() != 4) return result;

    const bool hasSavedMesh =
        cp.currentMesh.node.rows() > 0 && cp.currentMesh.elem.rows() > 0;
    if (hasSavedMesh) {
        Space checkpointSp;
        buildSpaceOnMesh(cp.currentMesh, ord, checkpointSp);
        result.usedSavedPhysicalMesh = true;
        result.savedMeshDof = checkpointSp.nDof;
        if (cp.U.rows() == checkpointSp.nDof) {
            const bool sameNodeShape =
                cp.currentMesh.node.rows() == targetSp.mesh.node.rows() &&
                cp.currentMesh.node.cols() == targetSp.mesh.node.cols();
            result.meshNodeDiff = sameNodeShape ?
                maxAbsMatrixDiff(cp.currentMesh.node, targetSp.mesh.node) : 1e300;
            result.sameElementTopology =
                sameElementMatrix(cp.currentMesh.elem, targetSp.mesh.elem);
            result.resizedTopology = !sameNodeShape || !result.sameElementTopology;
            const bool samePhysicalMesh =
                sameNodeShape && result.meshNodeDiff <= 1e-11 &&
                result.sameElementTopology;
            if (samePhysicalMesh && cp.U.rows() == targetSp.nDof) {
                U = cp.U;
            } else {
                U = interpolateDGToSpace(checkpointSp, cp.U, targetSp);
                result.interpolated = true;
            }
            result.compatible = true;
            return result;
        }
    }

    if (cp.U.rows() == targetSp.nDof) {
        U = cp.U;
        result.compatible = true;
    }
    return result;
}

RunCheckpoint makeRunCheckpoint(bool quick, int ord, int nFrames,
                                double tEnd, double h, double time,
                                double nextFrame, int step, int frame,
                                int remeshCount,
                                const std::vector<int>& milestoneDone,
                                const Mesh& referenceMesh,
                                const Space& currentSpace,
                                const ElasticSolid2D& solid,
                                const SolidALEMap& map,
                                const MatrixXd& U) {
    RunCheckpoint cp;
    cp.quick = quick;
    cp.ord = ord;
    cp.nFrames = nFrames;
    cp.tEnd = tEnd;
    cp.h = h;
    cp.time = time;
    cp.nextFrame = nextFrame;
    cp.step = step;
    cp.frame = frame;
    cp.remeshCount = remeshCount;
    cp.milestoneDone = milestoneDone;
    cp.referenceMesh = referenceMesh;
    cp.currentMesh = currentSpace.mesh;
    cp.solidReferenceMesh.node = solid.referenceNodes();
    cp.solidReferenceMesh.elem = solid.elements();
    cp.U = U;
    cp.solidNodes = solid.currentNodes();
    cp.solidVelocities = solid.velocities();
    cp.aleReferenceNodes = map.referenceNodes();
    return cp;
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
        out << "ALE_RUN_CHECKPOINT 4\n";
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
        out << "milestone_done " << cp.milestoneDone.size();
        for (int done : cp.milestoneDone) out << " " << done;
        out << "\n";
        writeMatrixXd(out, "mesh_node", cp.referenceMesh.node);
        writeMatrixXi(out, "mesh_elem", cp.referenceMesh.elem);
        writeMatrixXd(out, "current_mesh_node", cp.currentMesh.node);
        writeMatrixXi(out, "current_mesh_elem", cp.currentMesh.elem);
        writeMatrixXd(out, "solid_ref_node", cp.solidReferenceMesh.node);
        writeMatrixXi(out, "solid_ref_elem", cp.solidReferenceMesh.elem);
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

void pruneOldCheckpoints(const std::string& casePrefix, bool quick, int keepCount) {
    namespace fs = std::filesystem;
    if (keepCount < 0) keepCount = 0;
    std::vector<fs::path> existing;
    for (const auto& m : checkpointSchedule(quick)) {
        fs::path path = checkpointPath(casePrefix, quick, m.label);
        if (fs::exists(path)) existing.push_back(path);
    }
    int removeCount = static_cast<int>(existing.size()) - keepCount;
    for (int i = 0; i < removeCount; ++i) {
        std::error_code ec;
        fs::remove(existing[i], ec);
        if (ec) {
            std::cerr << "Warning: cannot remove old checkpoint "
                      << existing[i] << ": " << ec.message() << "\n";
        }
    }
}

std::optional<RunCheckpoint> loadLatestCheckpoint(const std::string& casePrefix, bool quick,
                                                  int ord, int nFrames,
                                                  double tEnd, double h,
                                                  int solidNodes,
                                                  bool allowExtension) {
    namespace fs = std::filesystem;
    std::vector<CheckpointMilestone> schedule = checkpointSchedule(quick);
    for (int i = static_cast<int>(schedule.size()) - 1; i >= 0; --i) {
        fs::path path = checkpointPath(casePrefix, quick, schedule[i].label);
        if (!fs::exists(path)) continue;
        RunCheckpoint cp;
        if (!loadCheckpoint(path, cp)) {
            std::cerr << "Warning: checkpoint unreadable, skipping " << path << "\n";
            continue;
        }
        if (!compatibleCheckpoint(cp, quick, ord, nFrames, tEnd, h, solidNodes,
                                  allowExtension)) {
            std::cerr << "Warning: checkpoint incompatible, skipping " << path << "\n";
            continue;
        }
        return cp;
    }
    return std::nullopt;
}

std::optional<RunCheckpoint> loadCheckpointByLabel(const std::string& casePrefix, bool quick,
                                                   const std::string& label,
                                                   int ord, int nFrames,
                                                   double tEnd, double h,
                                                   int solidNodes,
                                                   bool allowExtension) {
    namespace fs = std::filesystem;
    fs::path path = checkpointPath(casePrefix, quick, label);
    if (!fs::exists(path)) {
        std::cerr << "Warning: requested checkpoint does not exist: "
                  << path << "\n";
        return std::nullopt;
    }
    RunCheckpoint cp;
    if (!loadCheckpoint(path, cp)) {
        std::cerr << "Warning: requested checkpoint unreadable: "
                  << path << "\n";
        return std::nullopt;
    }
    if (!compatibleCheckpoint(cp, quick, ord, nFrames, tEnd, h, solidNodes,
                              allowExtension)) {
        std::cerr << "Warning: requested checkpoint incompatible: "
                  << path << "\n";
        return std::nullopt;
    }
    return cp;
}

} // namespace euler_ale
