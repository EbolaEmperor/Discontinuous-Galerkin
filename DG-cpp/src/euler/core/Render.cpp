#include "Render.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>

namespace euler {

using namespace Eigen;

bool readPPM(const std::string& path, int& W, int& H, std::vector<unsigned char>& img) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string magic;
    int maxv = 0;
    in >> magic >> W >> H >> maxv;
    in.get();
    if (magic != "P6" || W <= 0 || H <= 0) return false;
    img.resize(static_cast<size_t>(W) * H * 3);
    in.read(reinterpret_cast<char*>(img.data()), static_cast<std::streamsize>(img.size()));
    return true;
}

void writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img) {
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
}

std::vector<unsigned char> downsampleImage(const std::vector<unsigned char>& src,
                                           int srcW, int srcH, int factor,
                                           int& dstW, int& dstH) {
    if (factor <= 1) {
        dstW = srcW;
        dstH = srcH;
        return src;
    }
    if (srcW <= 0 || srcH <= 0 || srcW % factor != 0 || srcH % factor != 0) {
        dstW = 0;
        dstH = 0;
        return {};
    }

    dstW = srcW / factor;
    dstH = srcH / factor;
    std::vector<unsigned char> dst(static_cast<size_t>(dstW) * dstH * 3, 255);
    int samples = factor * factor;
    for (int y = 0; y < dstH; ++y) {
        for (int x = 0; x < dstW; ++x) {
            unsigned int acc[3] = {0, 0, 0};
            for (int sy = 0; sy < factor; ++sy) {
                int srcY = y * factor + sy;
                for (int sx = 0; sx < factor; ++sx) {
                    int srcX = x * factor + sx;
                    size_t srcIdx = (static_cast<size_t>(srcY) * srcW + srcX) * 3;
                    acc[0] += src[srcIdx];
                    acc[1] += src[srcIdx + 1];
                    acc[2] += src[srcIdx + 2];
                }
            }
            size_t dstIdx = (static_cast<size_t>(y) * dstW + x) * 3;
            dst[dstIdx] = static_cast<unsigned char>((acc[0] + samples / 2) / samples);
            dst[dstIdx + 1] = static_cast<unsigned char>((acc[1] + samples / 2) / samples);
            dst[dstIdx + 2] = static_cast<unsigned char>((acc[2] + samples / 2) / samples);
        }
    }
    return dst;
}

bool downsamplePPM(const std::string& srcPath, const std::string& dstPath, int factor) {
    int srcW = 0, srcH = 0;
    std::vector<unsigned char> src;
    if (!readPPM(srcPath, srcW, srcH, src)) return false;
    int dstW = 0, dstH = 0;
    std::vector<unsigned char> dst = downsampleImage(src, srcW, srcH, factor, dstW, dstH);
    if (dst.empty() || dstW <= 0 || dstH <= 0) return false;
    writePPM(dstPath, dstW, dstH, dst);
    return true;
}

void drawLine(std::vector<unsigned char>& img, int W, int H, int x0, int y0, int x1, int y1) {
    int dx = std::abs(x1 - x0), dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1, err = dx + dy;
    while (true) {
        if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) {
            size_t k = (static_cast<size_t>(y0) * W + x0) * 3;
            img[k] = static_cast<unsigned char>(0.55 * 30 + 0.45 * img[k]);
            img[k + 1] = static_cast<unsigned char>(0.55 * 220 + 0.45 * img[k + 1]);
            img[k + 2] = static_cast<unsigned char>(0.55 * 255 + 0.45 * img[k + 2]);
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

void overlayMesh(std::vector<unsigned char>& img, int W, int H,
                 const Mesh& mesh, double xa, double xb, double ya, double yb) {
    auto toPix = [&](const Vector2d& p, int& px, int& py) {
        px = static_cast<int>(std::lround((p.x() - xa) / (xb - xa) * W));
        py = static_cast<int>(std::lround((yb - p.y()) / (yb - ya) * H));
    };
    for (int t = 0; t < mesh.elem.rows(); ++t) {
        for (int k = 0; k < 3; ++k) {
            Vector2d A = mesh.node.row(mesh.elem(t, k));
            Vector2d B = mesh.node.row(mesh.elem(t, (k + 1) % 3));
            int ax, ay, bx, by;
            toPix(A, ax, ay);
            toPix(B, bx, by);
            drawLine(img, W, H, ax, ay, bx, by);
        }
    }
}

void overlayMesh(const std::string& path, const Mesh& mesh,
                 double xa, double xb, double ya, double yb) {
    int W = 0, H = 0;
    std::vector<unsigned char> img;
    if (!readPPM(path, W, H, img)) return;
    overlayMesh(img, W, H, mesh, xa, xb, ya, yb);
    writePPM(path, W, H, img);
}

} // namespace euler
