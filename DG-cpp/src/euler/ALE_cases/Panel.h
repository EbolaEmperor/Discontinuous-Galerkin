#ifndef PANEL_H
#define PANEL_H

#include "Core.h"

namespace euler_ale {

struct PanelMap {
    double t0 = 0.0;
    double t1 = 1.0;
    double q0 = 0.0;
    double q1 = 0.0;

    double mode(double x) const;
    double q(double t) const;
    double qd(double t) const;
    double topY(double x, double t) const;
    Vector2d refToPhys(const Vector2d& X, double t) const;
    Vector2d velocityAt(double x, double y, double t) const;
    void setLinearMotion(double a, double b, double qa, double qb);
    double maxMeshSpeed(double t) const;
};

struct PanelParams {
    double mass = 0.35;
    double stiffness = 85.0;
    double damping = 0.35;
    double pExt = 1.0;
    double qMin = -0.07;
    double qMax = 0.07;
};

double panelAcceleration(const ModalState& state, double fluidForce, const PanelParams& params);
ModalState advancePanelSymplectic(ModalState state, double fluidForce,
                                  const PanelParams& params, double dt);
int panelBoundaryTag(double x, double y, double t, const PanelMap& map);
double panelGeneralizedForce(const Space& sp, const MatrixXd& U, const PanelMap& map,
                             double time, double pExt, double* meanPressure = nullptr);

int runPanel(bool quick);

} // namespace euler_ale

#endif
