#include "Spoon.h"

#include "FEM.h"
#include "Mesh.h"
#include "NavierStokes.h"      // MeshLocator + basis-evaluation helpers

#include <algorithm>
#include <cmath>

namespace ns {

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::RowVectorXd;

namespace {
// Cubic smoothstep s(0)=0, s(1)=1, s'(0)=s'(1)=0.
inline double smoothstep(double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    return x * x * (3.0 - 2.0 * x);
}
} // namespace

Spoon::CrescentGeom Spoon::crescentGeom() const
{
    double r2f = std::min(std::max(crescentR2frac, 0.30), 0.97);
    double o   = std::min(std::max(crescentOffset, 0.15), 1.5);
    double etaH = (1.0 - r2f * r2f + o * o) / (2.0 * o);          // eta_horn / R1
    double erH  = std::sqrt(std::max(1e-4, 1.0 - etaH * etaH));   // er_horn  / R1
    double R1 = aRad / erH;                                       // horns land at +-aRad
    CrescentGeom g;
    g.R1 = R1; g.R2 = r2f * R1; g.eta1 = 0.0; g.eta2 = o * R1;
    return g;
}

void Spoon::build()
{
    wMark = dmark * dmark;
    double cb = std::cos(bearing0), sb = std::sin(bearing0);
    Vector2d er(cb, sb);              // radial (long-axis) unit at phi = 0
    Vector2d et(-sb, cb);             // tangential (thin-axis) unit at phi = 0

    std::vector<Vector2d> pts;
    if (!crescent) {
        // Ellipse: march a grid over the bounding box; keep interior markers.
        int nx = std::max(1, (int)std::ceil(aRad / dmark));
        int ny = std::max(1, (int)std::ceil(bThk / dmark));
        for (int i = -nx; i <= nx; ++i) {
            double xi = i * dmark;                    // along radial axis
            for (int j = -ny; j <= ny; ++j) {
                double eta = j * dmark;               // along thin axis
                double q = (xi * xi) / (aRad * aRad) + (eta * eta) / (bThk * bThk);
                if (q > 1.0) continue;
                pts.push_back((rPivot + xi) * er + eta * et);
            }
        }
    } else {
        // Crescent MOON (lune): inside outer circle C1(centre on the et axis, radius
        // R1) AND outside smaller circle C2 (radius R2=r2frac*R1, shifted +offset*R1
        // along +et).  The notch (opening) faces +et = the motion direction; the two
        // tips taper to sharp horns at +-aRad radially.  Local coords: xi along er
        // (radial), eta along et (tangential).
        CrescentGeom g = crescentGeom();   // R1,R2,eta1,eta2, and the centring shift
        double bb = g.R1 + 0.06;
        int nxi  = std::max(2, (int)std::ceil(2.0 * bb / dmark));
        int neta = std::max(2, (int)std::ceil((g.eta2 + g.R2 - (g.eta1 - g.R1)) / dmark));
        std::vector<Vector2d> loc;         // (xi, eta) interior points
        double etaLo = g.eta1 - g.R1, etaHi = g.eta2 + g.R2;
        for (int i = 0; i <= nxi; ++i) {
            double xi = -bb + 2.0 * bb * i / nxi;
            for (int j = 0; j <= neta; ++j) {
                double eta = etaLo + (etaHi - etaLo) * j / neta;
                double d1 = xi * xi + (eta - g.eta1) * (eta - g.eta1);
                double d2 = xi * xi + (eta - g.eta2) * (eta - g.eta2);
                if (d1 < g.R1 * g.R1 && d2 > g.R2 * g.R2) loc.emplace_back(xi, eta);
            }
        }
        // Centre the crescent on the blade centre (subtract the eta-centroid; xi is
        // ~0 by symmetry) so it sits at radius ~rPivot from the pivot.
        double mEta = 0.0;
        for (auto& p : loc) mEta += p.y();
        if (!loc.empty()) mEta /= (double)loc.size();
        for (auto& p : loc)
            pts.push_back((rPivot + p.x()) * er + (p.y() - mEta) * et);
        cresEtaShift_ = mEta;              // remember for outline()
    }
    r0.resize((int)pts.size(), 2);
    for (int k = 0; k < (int)pts.size(); ++k) r0.row(k) = pts[k].transpose();
    X.setZero(r0.rows(), 2);
    V.setZero(r0.rows(), 2);
    hint.assign(r0.rows(), -1);
}

void Spoon::update(double t)
{
    // ---- forcing envelope gain(t) ----
    if (t < tEnter) {
        gain = (tEnter > 0.0) ? smoothstep(t / tEnter) : 1.0;
    } else if (t <= tEndStroke()) {
        gain = 1.0;
    } else if (t < tDone()) {
        gain = (tLift > 0.0) ? 1.0 - smoothstep((t - tEndStroke()) / tLift) : 0.0;
    } else {
        gain = 0.0;
    }

    // ---- swept angle phi(t) and angular velocity omega(t): QUADRATIC profile ----
    // omega(tau) = omax * 4 tau (1 - tau), tau = u/T in [0,1]:  starts and ends at
    // zero, peaks omax = 1.5*sweep/T at mid-stroke.  Its integral
    //   phi(tau) = omax * T * (2 tau^2 - (4/3) tau^3)   sweeps 0 -> sweep exactly.
    double T = tStroke;
    double omax = peakOmega();
    double u = t - tEnter;                       // time into the stroke
    if (u <= 0.0) {
        phi = 0.0; omega = 0.0;
    } else if (u >= T) {
        phi = sweep; omega = 0.0;
    } else {
        double tau = u / T;
        omega = omax * 4.0 * tau * (1.0 - tau);
        phi   = omax * T * (2.0 * tau * tau - (4.0 / 3.0) * tau * tau * tau);
    }

    // ---- rigid rotation about the pivot ----
    double c = std::cos(phi), s = std::sin(phi);
    Vector2d piv(px, py);
    for (int k = 0; k < r0.rows(); ++k) {
        double rx = r0(k, 0), ry = r0(k, 1);
        double dx = c * rx - s * ry;             // Rot(phi) * r0
        double dy = s * rx + c * ry;
        X(k, 0) = px + dx;
        X(k, 1) = py + dy;
        // V = omega x (X - pivot) = omega * (-(X-piv).y, (X-piv).x)
        V(k, 0) = -omega * dy;
        V(k, 1) =  omega * dx;
    }
}

Vector2d Spoon::applyForcing(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                             const VectorXd& uField, const VectorXd& vField,
                             const MeshLocator& loc, double dt,
                             VectorXd& loadFx, VectorXd& loadFy)
{
    Vector2d net(0.0, 0.0);
    int Mm = M();
    if (Mm == 0 || gain <= 1e-6 || ibAlpha <= 0.0) return net;
    int locDof = fem.locDof;
    if ((int)hint.size() != Mm) hint.assign(Mm, -1);
    const double alpha = ibAlpha / dt;
    const double cap = std::max(0.0, ibForceCap);

    for (int k = 0; k < Mm; ++k) {
        double xk = X(k, 0), yk = X(k, 1);
        double lam[3];
        int e = loc.locate(mesh, xk, yk, hint[k], lam);
        if (e < 0) continue;                     // marker outside the mesh
        hint[k] = e;
        Vector3d L(lam[0], lam[1], lam[2]);
        RowVectorXd phi = fem.computeBasisValue_all(L.transpose()).row(0);
        double uh = 0.0, vh = 0.0;
        for (int i = 0; i < locDof; ++i) {
            int g = elem2dof(e, i);
            uh += phi(i) * uField(g);
            vh += phi(i) * vField(g);
        }
        double fx = alpha * gain * (V(k, 0) - uh) * wMark;
        double fy = alpha * gain * (V(k, 1) - vh) * wMark;
        if (cap > 0.0) {
            double nrm = std::hypot(fx, fy);
            if (nrm > cap) { double sc = cap / nrm; fx *= sc; fy *= sc; }
        }
        for (int i = 0; i < locDof; ++i) {
            int g = elem2dof(e, i);
            loadFx(g) += phi(i) * fx;
            loadFy(g) += phi(i) * fy;
        }
        net(0) += fx; net(1) += fy;
    }
    return net;
}

bool Spoon::contains(double x, double y) const
{
    // Transform (x,y) back to the blade's local frame (un-rotate by phi around pivot).
    double dx = x - px, dy = y - py;
    double c = std::cos(phi), s = std::sin(phi);
    double lx = c * dx + s * dy;   // local "rotated-back" x
    double ly = -s * dx + c * dy;  // local "rotated-back" y
    // Now express in the (er, et) frame at bearing0:
    double cb = std::cos(bearing0), sb = std::sin(bearing0);
    double xi  = cb * lx + sb * ly - rPivot;   // along radial (long axis), centred on blade
    double eta = -sb * lx + cb * ly;            // along tangential (thin axis)
    if (!crescent) {
        // Ellipse test
        double q = (xi * xi) / (aRad * aRad) + (eta * eta) / (bThk * bThk);
        return q <= 1.0;
    } else {
        // Lune test: inside outer circle C1, outside inner circle C2
        CrescentGeom g = crescentGeom();
        double etaL = eta + cresEtaShift_;   // undo the centring shift
        double d1 = xi * xi + (etaL - g.eta1) * (etaL - g.eta1);
        double d2 = xi * xi + (etaL - g.eta2) * (etaL - g.eta2);
        return (d1 < g.R1 * g.R1) && (d2 > g.R2 * g.R2);
    }
}

MatrixXd Spoon::outline(int nSeg) const
{
    if (nSeg < 4) nSeg = 4;
    double cb = std::cos(bearing0), sb = std::sin(bearing0);
    Vector2d er(cb, sb), et(-sb, cb);
    std::vector<Vector2d> off;                       // offsets from pivot at phi = 0
    if (!crescent) {
        for (int i = 0; i <= nSeg; ++i) {
            double th = 2.0 * M_PI * (double)i / nSeg;
            off.push_back((rPivot + aRad * std::cos(th)) * er + (bThk * std::sin(th)) * et);
        }
    } else {
        // Crescent moon outline: outer arc (circle C1, the back, between the horns)
        // then inner arc (circle C2, the concave notch), reversed, closing at the
        // horns.  theta measured from +et (top) clockwise; horns at +-th1 on C1.
        CrescentGeom g = crescentGeom();
        double s = cresEtaShift_;                    // same centring as build()
        double etaH = (g.R1 * g.R1 - g.R2 * g.R2 + g.eta2 * g.eta2) / (2.0 * g.eta2);
        double erH  = std::sqrt(std::max(1e-6, g.R1 * g.R1 - etaH * etaH));
        double th1 = std::atan2(erH, etaH - g.eta1);            // horn angle on C1
        double th2 = std::atan2(erH, etaH - g.eta2);            // horn angle on C2
        int na = std::max(3, nSeg / 2);
        auto add = [&](double cx_eta, double R, double th) {   // point on a circle (xi,eta)
            double xi = R * std::sin(th), eta = cx_eta + R * std::cos(th);
            off.push_back((rPivot + xi) * er + (eta - s) * et);
        };
        for (int i = 0; i <= na; ++i)                 // outer arc C1: th1 -> 2pi-th1 (through pi, the back)
            add(g.eta1, g.R1, th1 + (2.0 * M_PI - 2.0 * th1) * i / na);
        for (int i = 0; i <= na; ++i)                 // inner arc C2: 2pi-th2 -> th2 (reversed, the notch)
            add(g.eta2, g.R2, (2.0 * M_PI - th2) - (2.0 * M_PI - 2.0 * th2) * i / na);
        off.push_back(off.front());                   // close
    }
    double c = std::cos(phi), s = std::sin(phi);
    MatrixXd poly((int)off.size(), 2);
    for (int i = 0; i < (int)off.size(); ++i) {
        double dx = c * off[i].x() - s * off[i].y();
        double dy = s * off[i].x() + c * off[i].y();
        poly(i, 0) = px + dx;
        poly(i, 1) = py + dy;
    }
    return poly;
}

Vector2d Spoon::centre() const
{
    double cb = std::cos(bearing0), sb = std::sin(bearing0);
    Vector2d off = rPivot * Vector2d(cb, sb);
    double c = std::cos(phi), s = std::sin(phi);
    return Vector2d(px + c * off.x() - s * off.y(),
                    py + s * off.x() + c * off.y());
}

} // namespace ns
