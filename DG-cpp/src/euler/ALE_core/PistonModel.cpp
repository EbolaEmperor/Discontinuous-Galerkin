#include "PistonModel.h"

#include <algorithm>
#include <cmath>

namespace euler_ale {

using namespace Eigen;

double PistonMap::s(double t) const {
    if (!prescribed) {
        double den = std::max(1e-14, t1 - t0);
        double a = std::clamp((t - t0) / den, 0.0, 1.0);
        return (1.0 - a) * s0 + a * s1;
    }
    return amp * std::sin(2.0 * M_PI * t / period);
}

double PistonMap::sd(double t) const {
    if (!prescribed) {
        (void)t;
        return (s1 - s0) / std::max(1e-14, t1 - t0);
    }
    return amp * (2.0 * M_PI / period) * std::cos(2.0 * M_PI * t / period);
}

Vector2d PistonMap::refToPhys(const Vector2d& X, double t) const {
    double st = s(t);
    return Vector2d(st + (1.0 - st) * X.x(), X.y());
}

Vector2d PistonMap::velocityAt(double x, double, double t) const {
    double st = s(t);
    double X = (x - st) / std::max(1e-12, 1.0 - st);
    return Vector2d(sd(t) * (1.0 - X), 0.0);
}

double PistonMap::xLeft(double t) const { return s(t); }

void PistonMap::setLinearMotion(double a, double b, double xa, double xb) {
    prescribed = false;
    t0 = a;
    t1 = b;
    s0 = xa;
    s1 = xb;
}

double PistonMap::maxMeshSpeed(double t) const { return std::abs(sd(t)); }

int pistonBoundaryTag(double x, double y, double t, const PistonMap& map) {
    double xl = map.xLeft(t);
    double tol = 1e-7;
    if (std::abs(x - xl) < tol) return TAG_MOVING_WALL;
    if (std::abs(x - 1.0) < tol) return TAG_OUTFLOW;
    if (std::abs(y) < tol || std::abs(y - 1.0) < tol) return TAG_SLIP_WALL;
    return TAG_INTERIOR;
}

double pistonPressureForce(const Space& sp, const MatrixXd& U, double pExt,
                           double* meanPressure) {
    MatrixXd q1d;
    VectorXd w1d;
    sp.fem->quad1d(q1d, w1d);
    int nqe = static_cast<int>(w1d.size());
    int locDof = sp.fem->locDof;
    double force = 0.0, area = 0.0, pInt = 0.0;
    for (int e = 0; e < sp.edge.rows(); ++e) {
        if (sp.tag(e) != TAG_MOVING_WALL) continue;
        int t = (sp.e2s(e, 0) >= 0) ? sp.e2s(e, 0) : sp.e2s(e, 1);
        if (t < 0) continue;
        int n1 = sp.edge(e, 0), n2 = sp.edge(e, 1);
        EdgeOnElem eo = edgeOnElem(sp.mesh, t, n1, n2);
        Matrix<double, Dynamic, 4> Ue(locDof, 4);
        for (int i = 0; i < locDof; ++i) Ue.row(i) = U.row(sp.e2d(t, i));
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
            double ds = w1d(q) * eo.he;
            force += (p - pExt) * eo.nout.x() * ds;
            pInt += p * ds;
            area += ds;
        }
    }
    if (meanPressure) *meanPressure = (area > 0.0) ? pInt / area : pExt;
    return force;
}

double structureAcceleration(const PistonState& state, double fluidForce,
                             const FSIParams& params) {
    return (fluidForce - params.stiffness * state.x - params.damping * state.v) / params.mass;
}

PistonState advanceStructureSymplectic(PistonState state, double fluidForce,
                                       const FSIParams& params, double dt) {
    double a = structureAcceleration(state, fluidForce, params);
    state.v += dt * a;
    state.x += dt * state.v;
    if (state.x < params.xMin) {
        state.x = params.xMin;
        state.v = std::max(0.0, state.v);
    }
    if (state.x > params.xMax) {
        state.x = params.xMax;
        state.v = std::min(0.0, state.v);
    }
    return state;
}

} // namespace euler_ale
