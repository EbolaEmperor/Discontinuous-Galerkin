#include "DG.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace euler {

using namespace Eigen;

namespace {
void cmap(double t, Colormap cm, unsigned char& R, unsigned char& G, unsigned char& B) {
    t = std::min(1.0, std::max(0.0, t));
    double col[3];
    auto lerp = [](const double* a, const double* b, double s, double* o){ for (int k=0;k<3;++k) o[k]=a[k]+s*(b[k]-a[k]); };
    if (cm == CM_GRAY)      { col[0]=col[1]=col[2]=255.0*t; }
    else if (cm == CM_GRAY_INV) { col[0]=col[1]=col[2]=255.0*(1.0-t); }
    else if (cm == CM_COOLWARM) {
        const double lo[3]={59,76,192}, mid[3]={242,242,242}, hi[3]={180,4,38};
        if (t<0.5) lerp(lo,mid,t/0.5,col); else lerp(mid,hi,(t-0.5)/0.5,col);
    } else if (cm == CM_JET) {
        const double c0[3]={0,0,131},c1[3]={0,60,170},c2[3]={5,255,255},c3[3]={255,255,0},c4[3]={250,0,0};
        const double* C[5]={c0,c1,c2,c3,c4}; double s=t*4.0; int k=std::min(3,(int)s); lerp(C[k],C[k+1],s-k,col);
    } else if (cm == CM_INFERNO) {
        const double c0[3]={0,0,4},c1[3]={87,16,110},c2[3]={188,55,84},c3[3]={249,142,9},c4[3]={252,255,164};
        const double* C[5]={c0,c1,c2,c3,c4}; double s=t*4.0; int k=std::min(3,(int)s); lerp(C[k],C[k+1],s-k,col);
    } else {
        const double c0[3]={68,1,84},c1[3]={59,82,139},c2[3]={33,145,140},c3[3]={94,201,98},c4[3]={253,231,37};
        const double* C[5]={c0,c1,c2,c3,c4}; double s=t*4.0; int k=std::min(3,(int)s); lerp(C[k],C[k+1],s-k,col);
    }
    R=(unsigned char)std::lround(col[0]); G=(unsigned char)std::lround(col[1]); B=(unsigned char)std::lround(col[2]);
}
} // namespace

std::vector<unsigned char> renderScalarPPMImage(
                    FEM& fem, const Mesh& mesh,
                    const MatrixXi& elem2dof, const VectorXd& field, int W, int H,
                    double xmin, double xmax, double ymin, double ymax,
                    double vmin, double vmax, Colormap cm,
                    const std::function<bool(double, double)>& inDomain) {
    int locDof = fem.locDof;
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    std::vector<unsigned char> img((size_t)W * H * 3, 255);
    int NT = static_cast<int>(mesh.elem.rows());
    MatrixXi midx = numSplit3(fem.ord);
    VectorXd fe(locDof), mco(locDof);
    for (int t = 0; t < NT; ++t) {
        Vector2d P1 = mesh.node.row(mesh.elem(t,0)), P2 = mesh.node.row(mesh.elem(t,1)), P3 = mesh.node.row(mesh.elem(t,2));
        for (int i = 0; i < locDof; ++i) fe(i) = field(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        Vector2d v0 = P2 - P1, v1 = P3 - P1;
        double det = v0.x()*v1.y() - v1.x()*v0.y();
        if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;
        double txmin=std::min({P1.x(),P2.x(),P3.x()}), txmax=std::max({P1.x(),P2.x(),P3.x()});
        double tymin=std::min({P1.y(),P2.y(),P3.y()}), tymax=std::max({P1.y(),P2.y(),P3.y()});
        int col0=std::max(0,(int)std::floor((txmin-xmin)/dx)-1), col1=std::min(W-1,(int)std::ceil((txmax-xmin)/dx)+1);
        int row0=std::max(0,(int)std::floor((ymax-tymax)/dy)-1), row1=std::min(H-1,(int)std::ceil((ymax-tymin)/dy)+1);
        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double x = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double y = ymax - (row + 0.5) * dy;
                double px = x - P1.x(), py = y - P1.y();
                double l2 = (px*v1.y() - v1.x()*py)*invd, l3 = (v0.x()*py - px*v0.y())*invd, l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                if (inDomain && !inDomain(x, y)) continue;
                double val = 0.0;
                for (int j = 0; j < locDof; ++j) {
                    double s = mco(j);
                    for (int e = 0; e < midx(0,j); ++e) s *= l1;
                    for (int e = 0; e < midx(1,j); ++e) s *= l2;
                    for (int e = 0; e < midx(2,j); ++e) s *= l3;
                    val += s;
                }
                double s = (val - vmin) / (vmax - vmin);
                unsigned char R,G,B; cmap(s, cm, R, G, B);
                size_t idx = ((size_t)row*W + col)*3;
                img[idx]=R; img[idx+1]=G; img[idx+2]=B;
            }
        }
    }
    return img;
}

void writeScalarPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                    const MatrixXi& elem2dof, const VectorXd& field, int W, int H,
                    double xmin, double xmax, double ymin, double ymax,
                    double vmin, double vmax, Colormap cm,
                    const std::function<bool(double, double)>& inDomain) {
    std::vector<unsigned char> img =
        renderScalarPPMImage(fem, mesh, elem2dof, field, W, H, xmin, xmax,
                             ymin, ymax, vmin, vmax, cm, inDomain);
    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "writeScalarPPM: cannot open " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

void writeSchlierenPPM(const std::string& path, FEM& fem, const Mesh& mesh,
                       const MatrixXi& elem2dof, const VectorXd& rho, int W, int H,
                       double xmin, double xmax, double ymin, double ymax, double beta) {
    int locDof = fem.locDof, NT = static_cast<int>(mesh.elem.rows());
    MatrixXi midx = numSplit3(fem.ord);
    double dx = (xmax - xmin) / W, dy = (ymax - ymin) / H;
    auto inWin = [&](int t) {
        Vector2d c = (mesh.node.row(mesh.elem(t,0)) + mesh.node.row(mesh.elem(t,1)) + mesh.node.row(mesh.elem(t,2))).transpose() / 3.0;
        double m = 0.05 * std::max(xmax - xmin, ymax - ymin);
        return c.x() > xmin - m && c.x() < xmax + m && c.y() > ymin - m && c.y() < ymax + m;
    };
    VectorXd fe(locDof), mco(locDof);
    int ord = fem.ord;
    auto gradAt = [&](const MatrixXd& Dl, double l1, double l2, double l3, double& gx, double& gy) {
        double pw1[16], pw2[16], pw3[16];
        pw1[0] = pw2[0] = pw3[0] = 1.0;
        for (int e = 1; e <= ord; ++e) { pw1[e] = pw1[e-1]*l1; pw2[e] = pw2[e-1]*l2; pw3[e] = pw3[e-1]*l3; }
        double d1 = 0, d2 = 0, d3 = 0;
        for (int j = 0; j < locDof; ++j) {
            int a = midx(0, j), b = midx(1, j), c = midx(2, j);
            if (a > 0) d1 += mco(j) * a * pw1[a-1] * pw2[b]   * pw3[c];
            if (b > 0) d2 += mco(j) * b * pw1[a]   * pw2[b-1] * pw3[c];
            if (c > 0) d3 += mco(j) * c * pw1[a]   * pw2[b]   * pw3[c-1];
        }
        gx = Dl(0,0)*d1 + Dl(0,1)*d2 + Dl(0,2)*d3;
        gy = Dl(1,0)*d1 + Dl(1,1)*d2 + Dl(1,2)*d3;
    };
    double gmax = 1e-30;
    Vector3d lc(1.0/3, 1.0/3, 1.0/3);
    for (int t = 0; t < NT; ++t) {
        if (!inWin(t)) continue;
        for (int i = 0; i < locDof; ++i) fe(i) = rho(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        double gx, gy; gradAt(fem.Dlam[t], lc(0), lc(1), lc(2), gx, gy);
        gmax = std::max(gmax, std::hypot(gx, gy));
    }
    std::vector<unsigned char> img((size_t)W * H * 3, 255);
    for (int t = 0; t < NT; ++t) {
        if (!inWin(t)) continue;
        Vector2d P1 = mesh.node.row(mesh.elem(t,0)), P2 = mesh.node.row(mesh.elem(t,1)), P3 = mesh.node.row(mesh.elem(t,2));
        for (int i = 0; i < locDof; ++i) fe(i) = rho(elem2dof(t, i));
        mco.noalias() = fem.coef * fe;
        const MatrixXd& Dl = fem.Dlam[t];
        Vector2d v0 = P2 - P1, v1 = P3 - P1;
        double det = v0.x()*v1.y() - v1.x()*v0.y(); if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;
        double txmin=std::min({P1.x(),P2.x(),P3.x()}), txmax=std::max({P1.x(),P2.x(),P3.x()});
        double tymin=std::min({P1.y(),P2.y(),P3.y()}), tymax=std::max({P1.y(),P2.y(),P3.y()});
        int col0=std::max(0,(int)std::floor((txmin-xmin)/dx)-1), col1=std::min(W-1,(int)std::ceil((txmax-xmin)/dx)+1);
        int row0=std::max(0,(int)std::floor((ymax-tymax)/dy)-1), row1=std::min(H-1,(int)std::ceil((ymax-tymin)/dy)+1);
        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double x = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double y = ymax - (row + 0.5) * dy;
                double px = x - P1.x(), py = y - P1.y();
                double l2 = (px*v1.y() - v1.x()*py)*invd, l3 = (v0.x()*py - px*v0.y())*invd, l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                double gx, gy; gradAt(Dl, l1, l2, l3, gx, gy);
                double s = std::exp(-beta * std::hypot(gx, gy) / gmax);
                unsigned char v = (unsigned char)std::lround(255.0 * std::min(1.0, std::max(0.0, s)));
                size_t idx = ((size_t)row*W + col)*3;
                img[idx]=v; img[idx+1]=v; img[idx+2]=v;
            }
        }
    }
    std::ofstream out(path, std::ios::binary);
    if (!out) { std::cerr << "writeSchlierenPPM: cannot open " << path << "\n"; return; }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), (std::streamsize)img.size());
}

} // namespace euler
