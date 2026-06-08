#include "IBCoupler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ns {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::RowVectorXd;
using Eigen::Vector3d;

VectorXd rodArclengthWeights(const CosseratFilament& rod)
{
    int Nv = rod.N + 1;
    VectorXd ds(Nv); ds.setZero();
    // Edge lengths from current positions (consistent with the DG basis the
    // forces will pair against; rest-length weights would systematically bias
    // when the rod stretches even slightly).
    std::vector<double> lcur(rod.N);
    for (int i = 0; i < rod.N; ++i) {
        Vector2d e = rod.X.row(i + 1) - rod.X.row(i);
        lcur[i] = e.norm();
    }
    if (rod.N == 0) return ds;
    ds(0) = 0.5 * lcur[0];
    for (int i = 1; i < rod.N; ++i) ds(i) = 0.5 * (lcur[i - 1] + lcur[i]);
    ds(rod.N) = 0.5 * lcur[rod.N - 1];
    return ds;
}

RodVelocitySample meshToRod(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
                            const VectorXd& uField, const VectorXd& vField,
                            const MeshLocator& loc,
                            const CosseratFilament& rod,
                            std::vector<int>& hint)
{
    int Nv = rod.N + 1;
    int locDof = fem.locDof;
    if ((int)hint.size() != Nv) hint.assign(Nv, -1);

    RodVelocitySample s;
    s.uv.setZero(Nv, 2);
    s.alive.assign(Nv, 0);

    for (int k = 0; k < Nv; ++k) {
        double xk = rod.X(k, 0), yk = rod.X(k, 1);
        double lam[3];
        int t = loc.locate(mesh, xk, yk, hint[k], lam);
        if (t < 0) { s.alive[k] = 0; continue; }
        hint[k] = t;
        // Evaluate the dP_k basis at lam, accumulate u, v.
        Vector3d L(lam[0], lam[1], lam[2]);
        RowVectorXd phi = fem.computeBasisValue_all(L.transpose()).row(0);
        double u = 0.0, v = 0.0;
        for (int i = 0; i < locDof; ++i) {
            int g = elem2dof(t, i);
            u += phi(i) * uField(g);
            v += phi(i) * vField(g);
        }
        s.uv(k, 0) = u; s.uv(k, 1) = v;
        s.alive[k] = 1;
    }
    return s;
}

void rodToMesh(FEM& fem, const Mesh& mesh, const Eigen::MatrixXi& elem2dof,
               const MeshLocator& loc,
               const CosseratFilament& rod,
               const MatrixXd& F,
               const VectorXd& dsWeights,
               std::vector<int>& hint,
               VectorXd& loadX, VectorXd& loadY)
{
    int Nv = rod.N + 1;
    int locDof = fem.locDof;
    if ((int)hint.size() != Nv) hint.assign(Nv, -1);

    for (int k = 0; k < Nv; ++k) {
        double xk = rod.X(k, 0), yk = rod.X(k, 1);
        double lam[3];
        int t = loc.locate(mesh, xk, yk, hint[k], lam);
        if (t < 0) continue;     // node outside the mesh -> drop its load
        hint[k] = t;
        Vector3d L(lam[0], lam[1], lam[2]);
        RowVectorXd phi = fem.computeBasisValue_all(L.transpose()).row(0);
        double w = dsWeights(k);
        double Fx = F(k, 0) * w, Fy = F(k, 1) * w;
        for (int i = 0; i < locDof; ++i) {
            int g = elem2dof(t, i);
            loadX(g) += phi(i) * Fx;
            loadY(g) += phi(i) * Fy;
        }
    }
}

// ---------------------------------------------------------------------------
// PPM polyline overlay (read / paint / write back)
// ---------------------------------------------------------------------------

namespace {
// Skip whitespace and comment lines in PPM headers.
void skipWsComments(std::istream& in)
{
    int c = in.peek();
    while (c != EOF) {
        if (std::isspace(c)) { in.get(); c = in.peek(); }
        else if (c == '#') {
            std::string dummy; std::getline(in, dummy);
            c = in.peek();
        } else break;
    }
}
} // anon

void overlayPolylineOnPPM(const std::string& path,
                          const MatrixXd& X,
                          double xmin, double xmax, double ymin, double ymax,
                          int width, int rgb)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) { std::cerr << "overlayPolylineOnPPM: cannot open " << path << "\n"; return; }
    std::string magic; in >> magic;
    if (magic != "P6") {
        std::cerr << "overlayPolylineOnPPM: not a P6 PPM (" << path << ")\n";
        return;
    }
    skipWsComments(in);
    int W = 0, H = 0, maxv = 0;
    in >> W; skipWsComments(in);
    in >> H; skipWsComments(in);
    in >> maxv;
    in.get();   // single whitespace after maxv
    std::vector<unsigned char> img((size_t)W * H * 3);
    in.read(reinterpret_cast<char*>(img.data()), (std::streamsize)img.size());
    in.close();

    auto plot = [&](int col, int row, double a, unsigned char r, unsigned char g, unsigned char b) {
        if (col < 0 || col >= W || row < 0 || row >= H) return;
        if (a < 0) a = 0; if (a > 1) a = 1;
        size_t idx = ((size_t)row * W + col) * 3;
        img[idx]     = (unsigned char)((1 - a) * img[idx]     + a * r);
        img[idx + 1] = (unsigned char)((1 - a) * img[idx + 1] + a * g);
        img[idx + 2] = (unsigned char)((1 - a) * img[idx + 2] + a * b);
    };
    auto plotDisc = [&](double xp, double yp, int radius,
                        unsigned char r, unsigned char g, unsigned char b) {
        double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
        double col = (xp - xmin) / dx - 0.5;
        double row = (ymax - yp) / dy - 0.5;
        int c0 = (int)std::lround(col), r0 = (int)std::lround(row);
        for (int dr = -radius; dr <= radius; ++dr)
            for (int dc = -radius; dc <= radius; ++dc) {
                double dist = std::hypot((double)dc, (double)dr);
                if (dist <= radius + 0.5) {
                    double a = std::min(1.0, std::max(0.0, radius + 0.5 - dist));
                    plot(c0 + dc, r0 + dr, a, r, g, b);
                }
            }
    };
    auto plotLineThick = [&](double x0w, double y0w, double x1w, double y1w, int radius,
                             unsigned char r, unsigned char g, unsigned char b) {
        // Sample along the segment in pixel-space and stamp a disc at each step.
        double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
        double c0 = (x0w - xmin) / dx - 0.5;
        double rr0 = (ymax - y0w) / dy - 0.5;
        double c1 = (x1w - xmin) / dx - 0.5;
        double rr1 = (ymax - y1w) / dy - 0.5;
        double L = std::hypot(c1 - c0, rr1 - rr0);
        int steps = std::max(1, (int)std::ceil(L * 1.2));
        for (int s = 0; s <= steps; ++s) {
            double t = (double)s / steps;
            double cx = c0 + t * (c1 - c0);
            double cr = rr0 + t * (rr1 - rr0);
            // direct disc plotting in pixel space (no extra scaling)
            int ic = (int)std::lround(cx), ir = (int)std::lround(cr);
            for (int drow = -radius; drow <= radius; ++drow)
                for (int dcol = -radius; dcol <= radius; ++dcol) {
                    double dist = std::hypot((double)dcol, (double)drow);
                    if (dist <= radius + 0.5) {
                        double a = std::min(1.0, std::max(0.0, radius + 0.5 - dist));
                        plot(ic + dcol, ir + drow, a, r, g, b);
                    }
                }
        }
    };

    unsigned char R = (unsigned char)((rgb >> 16) & 0xFF);
    unsigned char G = (unsigned char)((rgb >> 8) & 0xFF);
    unsigned char B = (unsigned char)( rgb        & 0xFF);
    int N1 = (int)X.rows();
    int radius = std::max(1, width);
    for (int i = 0; i + 1 < N1; ++i) {
        plotLineThick(X(i, 0), X(i, 1), X(i + 1, 0), X(i + 1, 1), radius, R, G, B);
    }
    // Slightly larger marker at the clamped root and the free tip.
    if (N1 > 0) plotDisc(X(0,     0), X(0,     1), radius + 2, R, G, B);
    if (N1 > 1) plotDisc(X(N1 - 1, 0), X(N1 - 1, 1), radius + 1, R, G, B);

    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "overlayPolylineOnPPM: cannot rewrite " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

} // namespace ns
