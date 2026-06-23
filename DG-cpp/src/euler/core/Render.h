#ifndef EULER_RENDER_H
#define EULER_RENDER_H

#include "Mesh.h"

#include <string>
#include <vector>

namespace euler {

bool readPPM(const std::string& path, int& W, int& H, std::vector<unsigned char>& img);
void writePPM(const std::string& path, int W, int H, const std::vector<unsigned char>& img);
std::vector<unsigned char> downsampleImage(const std::vector<unsigned char>& src,
                                           int srcW, int srcH, int factor,
                                           int& dstW, int& dstH);
bool downsamplePPM(const std::string& srcPath, const std::string& dstPath, int factor);
void drawLine(std::vector<unsigned char>& img, int W, int H,
              int x0, int y0, int x1, int y1);
void overlayMesh(std::vector<unsigned char>& img, int W, int H,
                 const Mesh& mesh, double xa, double xb, double ya, double yb);
void overlayMesh(const std::string& path, const Mesh& mesh,
                 double xa, double xb, double ya, double yb);

} // namespace euler

#endif
