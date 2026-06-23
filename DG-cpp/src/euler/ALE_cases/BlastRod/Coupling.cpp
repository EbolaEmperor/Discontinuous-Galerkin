#include "Coupling.h"

#include <algorithm>

namespace euler_ale {

ModalState advanceBlastRodSymplectic(ModalState state, double fluidForce,
                                     const BlastRodParams& params, double dt) {
    double a = (fluidForce - params.stiffness * state.q - params.damping * state.qd) /
               params.mass;
    state.qd += dt * a;
    state.q += dt * state.qd;
    if (state.q < params.qMin) {
        state.q = params.qMin;
        state.qd = std::max(0.0, state.qd);
    }
    if (state.q > params.qMax) {
        state.q = params.qMax;
        state.qd = std::min(0.0, state.qd);
    }
    return state;
}

double blastRodGeneralizedForce(const Space& sp, const MatrixXd& U, const BlastRodMap& map,
                                double time, double pExt, double* meanPressure, double* drag) {
    (void)time;
    MatrixXd q1d;
    VectorXd w1d;
    sp.fem->quad1d(q1d, w1d);
    int nqe = static_cast<int>(w1d.size());
    int locDof = sp.fem->locDof;
    double force = 0.0, fx = 0.0, area = 0.0, pInt = 0.0;
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) != TAG_MOVING_WALL) continue;
        int tt = (sp.e2s(e, 0) >= 0) ? sp.e2s(e, 0) : sp.e2s(e, 1);
        if (tt < 0) continue;
        int n1 = sp.edge(e, 0), n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, tt, n1, n2);
        Matrix<double, Dynamic, 4> Ue(locDof, 4);
        for (int i = 0; i < locDof; ++i) Ue.row(i) = U.row(sp.e2d(tt, i));
        for (int q = 0; q < nqe; ++q) {
            double l1 = q1d(q, 0), l2 = q1d(q, 1);
            Vector3d lam = Vector3d::Zero();
            if (eo.et == 0) {
                lam(0) = (eo.dir == 0) ? l1 : l2;
                lam(1) = (eo.dir == 0) ? l2 : l1;
            } else if (eo.et == 1) {
                lam(1) = (eo.dir == 0) ? l1 : l2;
                lam(2) = (eo.dir == 0) ? l2 : l1;
            } else {
                lam(2) = (eo.dir == 0) ? l1 : l2;
                lam(0) = (eo.dir == 0) ? l2 : l1;
            }
            RowVectorXd phi = sp.fem->computeBasisValue_all(lam.transpose()).row(0);
            Vector4d Uq = (phi * Ue).transpose();
            double p = euler::pressure(Uq);
            Vector2d Pp = l1 * sp.mesh.node.row(n1).transpose() +
                          l2 * sp.mesh.node.row(n2).transpose();
            double ds = w1d(q) * eo.he;
            double dp = p - pExt;
            Vector2d dqShape = map.beamDqAtPhysical(Pp, time);
            force += dp * eo.nout.dot(dqShape) * ds;
            fx += dp * eo.nout.x() * ds;
            pInt += p * ds;
            area += ds;
        }
    }
    if (meanPressure) *meanPressure = (area > 0.0) ? pInt / area : pExt;
    if (drag) *drag = fx;
    return force;
}

double blastRodLoadSolidFromFluid(const Space& sp, const MatrixXd& U, ElasticSolid2D& solid,
                                  double pExt, double* meanPressure, double* drag) {
    return loadSolidFromFluidPressure(sp, U, solid, pExt, meanPressure, drag);
}

} // namespace euler_ale
