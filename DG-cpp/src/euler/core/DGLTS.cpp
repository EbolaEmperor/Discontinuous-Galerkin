#include "DGImpl.h"

#include <chrono>
#include <cmath>
#include <iostream>

namespace euler {

using namespace Eigen;

void EulerDG::inviscidResidualMaskedCompact(const MatrixXd& Uc, const std::vector<int>& cells,
                                            double t, const ExteriorStateFn& bc, MatrixXd& Rc) const {
    const Impl& P = *P_;
    const int nc = static_cast<int>(cells.size());
    Rc.resize((size_t)nc * locDof_, 4);
    const int nqv = static_cast<int>(P.wv.size());
    const bool hllcFlux = cfg_.use_hllc;
    parallel_ranges(nc, [&](int lo, int hi) {
        Vector4d Fx, Fy;
        MatrixXd Ue(locDof_, 4), Unb(locDof_, 4), G(2, locDof_);
        Matrix<double, Dynamic, 4> Re(locDof_, 4);
        for (int ci = lo; ci < hi; ++ci) {
            int tt = cells[ci];
            int cls = cellClass_[tt];
            Ue = Uc.middleRows((size_t)ci * locDof_, locDof_);
            Re.setZero();
            double area = fem_.area(tt);
            for (int q = 0; q < nqv; ++q) {
                Vector4d Uq = (P.phiV[q] * Ue).transpose();
                fluxes(Uq, Fx, Fy);
                G.noalias() = fem_.Dlam[tt] * P.dphiV[q];
                double wa = P.wv(q) * area;
                Re.noalias() += (wa * G.row(0).transpose()) * Fx.transpose();
                Re.noalias() += (wa * G.row(1).transpose()) * Fy.transpose();
            }
            for (int k = 0; k < 3; ++k) {
                const Impl::EE& r = P.ee[tt][k];
                bool interior = (r.nb != -1);
                if (interior) {
                    if (cellClass_[r.nb] == cls)
                        Unb = Uc.middleRows((size_t)g2c_[r.nb] * locDof_, locDof_);
                    else
                        for (int i = 0; i < locDof_; ++i) Unb.row(i) = ltsSnap_.row(elem2dof_(r.nb, i));
                }
                int tag = edgeTag_.size() ? edgeTag_(r.ei) : 0;
                for (int q = 0; q < P.nqe; ++q) {
                    const RowVectorXd& pa = P.ephi[r.et][r.dir][q];
                    Vector4d Um = (pa * Ue).transpose();
                    Vector4d Up;
                    if (interior) {
                        Up = (P.ephi[r.et_nb][r.dir_nb][q] * Unb).transpose();
                    } else {
                        double l1 = P.quad1d(q, 0), l2 = P.quad1d(q, 1);
                        Vector2d Pp = l1 * mesh_.node.row(r.n1).transpose() + l2 * mesh_.node.row(r.n2).transpose();
                        Up = bc ? bc(Pp.x(), Pp.y(), t, Um, r.nx, r.ny, tag) : Um;
                    }
                    Vector4d Hn = hllcFlux ? hllc(Um, Up, r.nx, r.ny) : rusanov(Um, Up, r.nx, r.ny);
                    double whe = P.w1d(q) * r.he;
                    Re.noalias() -= (whe * pa.transpose()) * Hn.transpose();
                }
            }
            Rc.middleRows((size_t)ci * locDof_, locDof_) = Re;
        }
    });
}

void EulerDG::accumulateReflux(const MatrixXd& Uc, int advClass, double weight,
                               std::vector<MatrixXd>& reg) const {
    const Impl& P = *P_;
    const bool hllcFlux = cfg_.use_hllc;
    MatrixXd Uadv(locDof_, 4), Uoth(locDof_, 4);
    for (const LtsEdge& le : ltsEdges_) {
        int cc = cellClass_[le.coarse];
        bool advIsCoarse;
        if (advClass == cc)          advIsCoarse = true;
        else if (advClass == cc - 1) advIsCoarse = false;
        else continue;
        const Impl::EE& eeC = P.ee[le.coarse][le.coarseK];
        int advCell = advIsCoarse ? le.coarse : le.fine;
        const Impl::EE& eeAdv = advIsCoarse ? eeC : P.ee[le.fine][le.fineK];
        const Impl::EE& eeOth = advIsCoarse ? P.ee[le.fine][le.fineK] : eeC;
        int othCell = advIsCoarse ? le.fine : le.coarse;
        int adv0 = g2c_[advCell] * locDof_;
        for (int i = 0; i < locDof_; ++i) {
            Uadv.row(i) = Uc.row(adv0 + i);
            Uoth.row(i) = ltsSnap_.row(elem2dof_(othCell, i));
        }
        Matrix<double, Dynamic, 4> contrib = Matrix<double, Dynamic, 4>::Zero(locDof_, 4);
        for (int q = 0; q < P.nqe; ++q) {
            Vector4d Um = (P.ephi[eeAdv.et][eeAdv.dir][q] * Uadv).transpose();
            Vector4d Up = (P.ephi[eeOth.et][eeOth.dir][q] * Uoth).transpose();
            Vector4d Hn = hllcFlux ? hllc(Um, Up, eeAdv.nx, eeAdv.ny) : rusanov(Um, Up, eeAdv.nx, eeAdv.ny);
            const RowVectorXd& phiC = P.ephi[eeC.et][eeC.dir][q];
            double whe = P.w1d(q) * eeAdv.he;
            contrib.noalias() += (weight * whe * phiC.transpose()) * Hn.transpose();
        }
        for (int i = 0; i < locDof_; ++i) reg[cc].row(elem2dof_(le.coarse, i)) += contrib.row(i);
    }
}

bool EulerDG::advanceMaskedARS(const std::vector<int>& cells, int advClass, double dt, double tN,
                               const ExteriorStateFn& bc, std::vector<MatrixXd>& reg) {
    const double g   = 1.0 - std::sqrt(2.0) / 2.0;
    const double del = -1.0 / std::sqrt(2.0);
    const bool avOn = cfg_.use_av && (advClass == 0);
    const int nc = static_cast<int>(cells.size());
    if (nc == 0) return true;

    for (int ci = 0; ci < nc; ++ci) g2c_[cells[ci]] = ci;

    MatrixXd Uc((size_t)nc * locDof_, 4);
    parallel_ranges(nc, [&](int lo, int hi) {
        for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
            for (int i = 0; i < locDof_; ++i) Uc.row((size_t)ci * locDof_ + i) = U_.row(elem2dof_(t, i)); }
    });

    std::vector<int> actRow;
    if (avOn && activeDofs_.size()) {
        int na = (int)activeDofs_.size(); actRow.resize(na);
        for (int a = 0; a < na; ++a) { int gd = activeDofs_[a]; int ci = g2c_[gd / locDof_];
            actRow[a] = (ci >= 0) ? ci * locDof_ + gd % locDof_ : -1; }
    }

    auto massC = [&](const MatrixXd& X, MatrixXd& MX) {
        MX.resize((size_t)nc * locDof_, 4);
        parallel_ranges(nc, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Ue(locDof_, 4);
            for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
                Ue.noalias() = (fem_.area(t) * P_->Mref) * X.middleRows((size_t)ci * locDof_, locDof_);
                MX.middleRows((size_t)ci * locDof_, locDof_) = Ue; }
        });
    };
    auto massInvC = [&](MatrixXd& R) {
        parallel_ranges(nc, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Re(locDof_, 4);
            for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
                Re.noalias() = (MrefInv_ / fem_.area(t)) * R.middleRows((size_t)ci * locDof_, locDof_);
                R.middleRows((size_t)ci * locDof_, locDof_) = Re; }
        });
    };
    auto limitC = [&](MatrixXd& X) {
        parallel_ranges(nc, [&](int lo, int hi) {
            Matrix<double, Dynamic, 4> Ue(locDof_, 4);
            for (int ci = lo; ci < hi; ++ci) {
                Ue = X.middleRows((size_t)ci * locDof_, locDof_);
                limitCellCore(Ue);
                X.middleRows((size_t)ci * locDof_, locDof_) = Ue; }
        });
    };
    const int na = avOn ? (int)activeDofs_.size() : 0;
    auto activeSolve = [&](const MatrixXd& Baa) -> MatrixXd {
        MatrixXd Xaa(na, 4);
        std::array<std::thread, 4> th;
        for (int c = 0; c < 4; ++c) th[c] = std::thread([&, c] { Xaa.col(c) = ldltA_.solve(Baa.col(c)); });
        for (auto& tt : th) tt.join();
        return Xaa;
    };
    auto Fload = [&](const MatrixXd& X, double tt) { MatrixXd R; inviscidResidualMaskedCompact(X, cells, tt, bc, R); return R; };
    auto gatherActive = [&](const MatrixXd& Bc, MatrixXd& Baa) {
        for (int a = 0; a < na; ++a) { if (actRow[a] >= 0) Baa.row(a) = Bc.row(actRow[a]); else Baa.row(a).setZero(); }
    };

    MatrixXd MUn; massC(Uc, MUn);
    MatrixXd f1 = Fload(Uc, tN);
    accumulateReflux(Uc, advClass, del * dt, reg);

    MatrixXd B1 = MUn + (g * dt) * f1;
    MatrixXd U2 = B1; massInvC(U2);
    MatrixXd Xaa1;
    if (na) {
        MatrixXd Baa(na, 4); gatherActive(B1, Baa);
        Xaa1 = activeSolve(Baa);
        for (int a = 0; a < na; ++a) if (actRow[a] >= 0) U2.row(actRow[a]) = Xaa1.row(a);
    }
    if (cfg_.use_positivity) limitC(U2);

    MatrixXd f2 = Fload(U2, tN + g * dt);
    accumulateReflux(U2, advClass, (1.0 - del) * dt, reg);

    MatrixXd rhs3 = MUn + dt * (del * f1 + (1.0 - del) * f2);
    MatrixXd U3 = rhs3; massInvC(U3);
    if (na) {
        MatrixXd AU = Aaa_ * Xaa1;
        MatrixXd Baa(na, 4);
        for (int a = 0; a < na; ++a) {
            Eigen::RowVector4d b; if (actRow[a] >= 0) b = rhs3.row(actRow[a]); else b.setZero();
            Baa.row(a) = b - ((1.0 - g) * dt) * AU.row(a);
        }
        MatrixXd Xaa3 = activeSolve(Baa);
        for (int a = 0; a < na; ++a) if (actRow[a] >= 0) U3.row(actRow[a]) = Xaa3.row(a);
    }
    if (cfg_.use_positivity) limitC(U3);

    parallel_ranges(nc, [&](int lo, int hi) {
        for (int ci = lo; ci < hi; ++ci) { int t = cells[ci];
            for (int i = 0; i < locDof_; ++i) U_.row(elem2dof_(t, i)) = U3.row((size_t)ci * locDof_ + i); }
    });
    for (int ci = 0; ci < nc; ++ci) g2c_[cells[ci]] = -1;
    return true;
}

void EulerDG::ltsIntegrate(int level, double dt, double tN, const ExteriorStateFn& bc,
                           std::vector<MatrixXd>& reg, int levels) {
    if (level >= 1) reg[level].setZero(nDof_, 4);
    advanceMaskedARS(classCells_[level], level, dt, tN, bc, reg);
    if (level == 0) return;
    ltsIntegrate(level - 1, 0.5 * dt, tN,            bc, reg, levels);
    ltsIntegrate(level - 1, 0.5 * dt, tN + 0.5 * dt, bc, reg, levels);
    const std::vector<int>& cl = classCells_[level];
    parallel_ranges(static_cast<int>(cl.size()), [&](int lo, int hi) {
        Matrix<double, Dynamic, 4> Re(locDof_, 4);
        for (int ci = lo; ci < hi; ++ci) { int t = cl[ci];
            for (int i = 0; i < locDof_; ++i) Re.row(i) = reg[level].row(elem2dof_(t, i));
            Re = (MrefInv_ / fem_.area(t)) * Re;
            for (int i = 0; i < locDof_; ++i) U_.row(elem2dof_(t, i)) += Re.row(i); }
    });
    if (cfg_.use_positivity) positivityLimitCells(U_, cl);
}

bool EulerDG::stepLTS(double dtMacro, double tEnd, const ExteriorStateFn& bc,
                      const std::vector<int>& cellClass, int levels) {
    cellClass_ = cellClass;
    if ((int)g2c_.size() != NT_) g2c_.assign(NT_, -1);
    ltsSnap_ = U_;
    const double tN  = tEnd - dtMacro;
    const double g   = 1.0 - std::sqrt(2.0) / 2.0;
    const double dt0 = dtMacro / static_cast<double>(1 << (levels - 1));

    classCells_.assign(levels, {});
    for (int t = 0; t < NT_; ++t) classCells_[cellClass_[t]].push_back(t);

    ltsEdges_.clear();
    for (int t = 0; t < NT_; ++t)
        for (int k = 0; k < 3; ++k) {
            int nb = P_->ee[t][k].nb;
            if (nb < 0 || t > nb) continue;
            if (cellClass_[t] == cellClass_[nb]) continue;
            int coarse = cellClass_[t] > cellClass_[nb] ? t : nb;
            int fine   = cellClass_[t] > cellClass_[nb] ? nb : t;
            int ck = -1, fk = -1, ei = P_->ee[t][k].ei;
            for (int kk = 0; kk < 3; ++kk) { if (P_->ee[coarse][kk].ei == ei) ck = kk; if (P_->ee[fine][kk].ei == ei) fk = kk; }
            ltsEdges_.push_back({coarse, fine, ck, fk});
        }

    static const bool LP = std::getenv("LTSPROF") != nullptr;
    static double tAV = 0, tInt = 0; static int np = 0;
    auto cA = std::chrono::high_resolution_clock::now();
    if (cfg_.use_av) {
        computeViscosity(U_, epsFrozen_);
        for (int t = 0; t < NT_; ++t) if (cellClass_[t] != 0) epsFrozen_(t) = 0.0;
        refreshImplicit(epsFrozen_, g * dt0);
    }
    if (LP) tAV += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - cA).count();

    std::vector<MatrixXd> reg(levels);
    auto c0 = std::chrono::high_resolution_clock::now();
    ltsIntegrate(levels - 1, dtMacro, tN, bc, reg, levels);
    if (LP) {
        auto c1 = std::chrono::high_resolution_clock::now();
        tInt += std::chrono::duration<double>(c1 - c0).count(); ++np;
        if (np % 50 == 0) std::cerr << "[ltsprof] AV-refresh=" << (1000*tAV/np) << " integrate="
                                    << (1000*tInt/np) << " ms/macro (" << np << " macros)\n";
    }
    ++nstep_;
    return true;
}

} // namespace euler
