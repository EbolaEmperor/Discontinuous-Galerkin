#ifndef EXACT_SOLUTION_CH_H
#define EXACT_SOLUTION_CH_H

#include <cmath>

// ---------------------------------------------------------------------------
// Manufactured solution for the Cahn-Hilliard convergence test on [0,1]^2 with
// no-flux (homogeneous Neumann) boundary conditions.
//
//   c_e(x,y,t) = cos(pi x) cos(pi y) * exp(-t)
//   mu_e       = F'(c_e) - eps^2 * Lap(c_e),     F'(c) = c^3 - c
//
// Because cos(pi x) (resp. cos(pi y)) has zero normal derivative at x=0,1
// (resp. y=0,1), BOTH d_n c_e = 0 and d_n mu_e = 0 on the boundary -- so the
// interior-only SIPG discretisation is consistent and the convergence rates are
// clean.  A forcing term is added ONLY to the c-equation so that c_e solves it
// exactly:
//
//   dc/dt = M Lap(mu) + g,   g(x,y,t) = d_t c_e - M Lap(mu_e).
//
// Using Lap(phi) = -2 pi^2 phi (phi = cos(pi x)cos(pi y)), a'(t) = -a(t), and
// Lap(u^3) = 6 u |grad u|^2 + 3 u^2 Lap(u):
//   d_t c_e   = -c_e
//   Lap(c_e)  = -2 pi^2 c_e
//   mu_e      = c_e^3 + (2 pi^2 eps^2 - 1) c_e
//   Lap(mu_e) = 6 c_e a^2 |grad phi|^2 - 6 pi^2 c_e^3 - 2 pi^2 (2 pi^2 eps^2 - 1) c_e
//   |grad phi|^2 = pi^2 [ sin^2(pi x) cos^2(pi y) + cos^2(pi x) sin^2(pi y) ]
// ---------------------------------------------------------------------------
struct ExactSolutionCH {
    double mob;   // mobility M
    double eps2;  // eps^2

    explicit ExactSolutionCH(double mob_ = 1.0, double eps2_ = 0.01)
        : mob(mob_), eps2(eps2_) {}

    double a(double t) const { return std::exp(-t); }

    double c(double x, double y, double t) const {
        return std::cos(M_PI * x) * std::cos(M_PI * y) * a(t);
    }

    void grad(double x, double y, double t, double& cx, double& cy) const {
        double at = a(t);
        cx = -M_PI * std::sin(M_PI * x) * std::cos(M_PI * y) * at;
        cy = -M_PI * std::cos(M_PI * x) * std::sin(M_PI * y) * at;
    }

    // Forcing term for the c-equation:  g = d_t c_e - M * Lap(mu_e).
    double source(double x, double y, double t) const {
        const double pi2 = M_PI * M_PI;
        double at = a(t);
        double cx = std::cos(M_PI * x), sx = std::sin(M_PI * x);
        double cy = std::cos(M_PI * y), sy = std::sin(M_PI * y);
        double ce = cx * cy * at;                                  // c_e
        double gradphi2 = pi2 * (sx * sx * cy * cy + cx * cx * sy * sy); // |grad phi|^2
        double lap_mu = 6.0 * ce * (at * at) * gradphi2
                      - 6.0 * pi2 * ce * ce * ce
                      - 2.0 * pi2 * (2.0 * pi2 * eps2 - 1.0) * ce;
        double dct = -ce;                                          // d_t c_e
        return dct - mob * lap_mu;
    }

    // The source is separable in time:  g(x,y,t) = a(t) * g1(x,y) + a(t)^3 * g3(x,y).
    // The convergence driver precomputes the two space loads \int g1 phi, \int g3 phi
    // ONCE and rescales them per step, instead of re-assembling g every step.
    double g1(double x, double y) const {
        const double pi2 = M_PI * M_PI;
        double phi = std::cos(M_PI * x) * std::cos(M_PI * y);
        return phi * (-1.0 + 2.0 * pi2 * mob * (2.0 * pi2 * eps2 - 1.0));
    }
    double g3(double x, double y) const {
        const double pi2 = M_PI * M_PI;
        double cx = std::cos(M_PI * x), sx = std::sin(M_PI * x);
        double cy = std::cos(M_PI * y), sy = std::sin(M_PI * y);
        double phi = cx * cy;
        double gradphi2 = pi2 * (sx * sx * cy * cy + cx * cx * sy * sy);
        return 6.0 * mob * (pi2 * phi * phi * phi - phi * gradphi2);
    }
};

#endif
