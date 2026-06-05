#include "CahnHilliard.h"
#include "Quadrature.h"

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

// ---------------------------------------------------------------------------
// Quadrature helper: a rule that integrates degree `deg` exactly (capped at the
// 14th-order rule available in Quadrature).
// ---------------------------------------------------------------------------
namespace {
void quadDeg(int deg, MatrixXd& quadL, VectorXd& w) {
    int order = std::min(14, std::max(1, deg));
    Quadrature::quadpts2_my(order, quadL, w);
}
}

SparseMatrix<double> assembleMass_DG2D(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;

    MatrixXd quadL;
    VectorXd w;
    // phi_i * phi_j has degree 2*ord; the FEM default rule (2*(ord+1)) is enough.
    fem.quad2d(quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<MatrixXd> phi_q(nq);
    for (int q = 0; q < nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q)); // 1 x locDof
    }

    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(NT) * locDof * locDof);

    MatrixXd Me(locDof, locDof);
    for (int t = 0; t < NT; ++t) {
        Me.setZero();
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) {
            const MatrixXd& phi = phi_q[q]; // 1 x locDof
            Me.noalias() += (w(q) * area) * (phi.transpose() * phi);
        }
        for (int i = 0; i < locDof; ++i) {
            int gi = elem2dof(t, i);
            for (int j = 0; j < locDof; ++j) {
                trip.emplace_back(gi, elem2dof(t, j), Me(i, j));
            }
        }
    }

    SparseMatrix<double> M(nDof, nDof);
    M.setFromTriplets(trip.begin(), trip.end());
    return M;
}

VectorXd assembleNonlinearCH(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                             const VectorXd& c) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    VectorXd b = VectorXd::Zero(nDof);

    MatrixXd quadL;
    VectorXd w;
    // Integrand F'(c_h)*phi_i = (c_h^3 - c_h)*phi_i has degree 4*ord
    // (c_h^3 is 3*ord, times the degree-ord test function); integrate it accurately.
    quadDeg(4 * fem.ord + 2, quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<MatrixXd> phi_q(nq);
    for (int q = 0; q < nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q)); // 1 x locDof
    }

    VectorXd celem(locDof);
    VectorXd be(locDof);
    for (int t = 0; t < NT; ++t) {
        for (int i = 0; i < locDof; ++i) celem(i) = c(elem2dof(t, i));
        double area = fem.area(t);
        be.setZero();
        for (int q = 0; q < nq; ++q) {
            const MatrixXd& phi = phi_q[q]; // 1 x locDof
            double ch = phi.row(0).dot(celem);
            double f = ch * ch * ch - ch; // F'(c) = c^3 - c
            be.noalias() += (w(q) * area * f) * phi.transpose();
        }
        for (int i = 0; i < locDof; ++i) b(elem2dof(t, i)) += be(i);
    }
    return b;
}

double computeEnergyCH(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                       const VectorXd& c, const SparseMatrix<double>& A, double eps2) {
    int NT = mesh.elem.rows();
    int locDof = fem.locDof;

    // Gradient/interface energy: (eps^2/2) * c^T A c.
    double gradEnergy = 0.5 * eps2 * (c.dot(A * c));

    // Bulk (double-well) energy: sum_K \int_K F(c_h),  F(c) = (1/4)(c^2-1)^2 (degree 4*ord).
    MatrixXd quadL;
    VectorXd w;
    quadDeg(4 * fem.ord + 2, quadL, w);
    int nq = static_cast<int>(w.size());

    std::vector<MatrixXd> phi_q(nq);
    for (int q = 0; q < nq; ++q) {
        phi_q[q] = fem.computeBasisValue_all(quadL.row(q));
    }

    double bulk = 0.0;
    VectorXd celem(locDof);
    for (int t = 0; t < NT; ++t) {
        for (int i = 0; i < locDof; ++i) celem(i) = c(elem2dof(t, i));
        double area = fem.area(t);
        for (int q = 0; q < nq; ++q) {
            double ch = phi_q[q].row(0).dot(celem);
            double Fc = 0.25 * (ch * ch - 1.0) * (ch * ch - 1.0);
            bulk += w(q) * area * Fc;
        }
    }
    return gradEnergy + bulk;
}

double computeMassCH(const SparseMatrix<double>& Mm, const VectorXd& c) {
    // 1^T Mm c = sum over all entries of (Mm c).
    return (Mm * c).sum();
}

VectorXd initSpinodal(int nDof, double cbar, double amp, unsigned seed) {
    VectorXd c(nDof);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < nDof; ++i) {
        c(i) = cbar + amp * dist(gen);
    }
    return c;
}

// ---------------------------------------------------------------------------
// Visualisation
// ---------------------------------------------------------------------------
namespace {
// Diverging "coolwarm" colormap: cmin -> blue, midpoint -> light grey, cmax -> red.
void coolwarm(double v, double cmin, double cmax,
              unsigned char& r, unsigned char& g, unsigned char& b) {
    double t = (v - cmin) / (cmax - cmin);
    t = std::min(1.0, std::max(0.0, t));
    // Control colours.
    const double lo[3] = { 59.0, 76.0, 192.0};
    const double mid[3] = {242.0, 242.0, 242.0};
    const double hi[3] = {180.0, 4.0, 38.0};
    double col[3];
    if (t < 0.5) {
        double s = t / 0.5;
        for (int k = 0; k < 3; ++k) col[k] = lo[k] + s * (mid[k] - lo[k]);
    } else {
        double s = (t - 0.5) / 0.5;
        for (int k = 0; k < 3; ++k) col[k] = mid[k] + s * (hi[k] - mid[k]);
    }
    r = static_cast<unsigned char>(std::lround(col[0]));
    g = static_cast<unsigned char>(std::lround(col[1]));
    b = static_cast<unsigned char>(std::lround(col[2]));
}
}

void writeFramePPM(const std::string& path, const Mesh& mesh,
                   const MatrixXi& elem2dof, const VectorXd& c,
                   int Npix, double cmin, double cmax) {
    // visualisation uses linear (vertex) interpolation per element
    const int W = Npix, H = Npix;

    // Use a SQUARE sampling window centred on the mesh bounding box, so a non-square
    // domain (e.g. a hexagon) is rendered with correct aspect (uncovered corners stay
    // background). For a unit square / inscribed disk this is the bounding box itself.
    double bxmin = mesh.node.col(0).minCoeff();
    double bxmax = mesh.node.col(0).maxCoeff();
    double bymin = mesh.node.col(1).minCoeff();
    double bymax = mesh.node.col(1).maxCoeff();
    double cx = 0.5 * (bxmin + bxmax), cy = 0.5 * (bymin + bymax);
    double half = 0.5 * std::max(bxmax - bxmin, bymax - bymin);
    double xmin = cx - half, xmax = cx + half;
    double ymin = cy - half, ymax = cy + half;
    double dx = (xmax - xmin) / W;
    double dy = (ymax - ymin) / H;

    // Background = midpoint colour (in case of any uncovered pixel).
    std::vector<unsigned char> img(static_cast<size_t>(W) * H * 3, 242);

    int NT = mesh.elem.rows();
    for (int t = 0; t < NT; ++t) {
        Vector2d P1 = mesh.node.row(mesh.elem(t, 0));
        Vector2d P2 = mesh.node.row(mesh.elem(t, 1));
        Vector2d P3 = mesh.node.row(mesh.elem(t, 2));
        // This element's nodal values at the three vertices (first 3 local dofs).
        double c1 = c(elem2dof(t, 0));
        double c2 = c(elem2dof(t, 1));
        double c3 = c(elem2dof(t, 2));

        Vector2d v0 = P2 - P1;
        Vector2d v1 = P3 - P1;
        double det = v0.x() * v1.y() - v1.x() * v0.y();
        if (std::abs(det) < 1e-300) continue;
        double invd = 1.0 / det;

        double txmin = std::min({P1.x(), P2.x(), P3.x()});
        double txmax = std::max({P1.x(), P2.x(), P3.x()});
        double tymin = std::min({P1.y(), P2.y(), P3.y()});
        double tymax = std::max({P1.y(), P2.y(), P3.y()});

        int col0 = std::max(0, static_cast<int>(std::floor((txmin - xmin) / dx)) - 1);
        int col1 = std::min(W - 1, static_cast<int>(std::ceil((txmax - xmin) / dx)) + 1);
        // y runs downwards in image space: row 0 is the top (y = ymax).
        int row0 = std::max(0, static_cast<int>(std::floor((ymax - tymax) / dy)) - 1);
        int row1 = std::min(H - 1, static_cast<int>(std::ceil((ymax - tymin) / dy)) + 1);

        const double tol = -1e-9;
        for (int col = col0; col <= col1; ++col) {
            double x = xmin + (col + 0.5) * dx;
            for (int row = row0; row <= row1; ++row) {
                double y = ymax - (row + 0.5) * dy;
                double px = x - P1.x();
                double py = y - P1.y();
                double l2 = (px * v1.y() - v1.x() * py) * invd;
                double l3 = (v0.x() * py - px * v0.y()) * invd;
                double l1 = 1.0 - l2 - l3;
                if (l1 < tol || l2 < tol || l3 < tol) continue;
                double val = l1 * c1 + l2 * c2 + l3 * c3;
                unsigned char r, g, b;
                coolwarm(val, cmin, cmax, r, g, b);
                size_t idx = (static_cast<size_t>(row) * W + col) * 3;
                img[idx] = r; img[idx + 1] = g; img[idx + 2] = b;
            }
        }
    }

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        std::cerr << "writeFramePPM: cannot open " << path << "\n";
        return;
    }
    out << "P6\n" << W << " " << H << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data()), static_cast<std::streamsize>(img.size()));
    if (!out) std::cerr << "writeFramePPM: write failed for " << path << "\n";
}

// ---------------------------------------------------------------------------
// Generic scalar load  L(i) = sum_K \int_K f(x,y) phi_i
// ---------------------------------------------------------------------------
VectorXd assembleScalarLoad(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                            const std::function<double(double, double)>& f) {
    int NT = mesh.elem.rows();
    int nDof = elem2dof.maxCoeff() + 1;
    int locDof = fem.locDof;
    VectorXd L = VectorXd::Zero(nDof);

    MatrixXd quadL;
    VectorXd w;
    quadDeg(2 * fem.ord + 8, quadL, w); // f is smooth but non-polynomial: integrate accurately
    int nq = static_cast<int>(w.size());
    std::vector<MatrixXd> phi_q(nq);
    for (int q = 0; q < nq; ++q) phi_q[q] = fem.computeBasisValue_all(quadL.row(q));

    VectorXd Le(locDof);
    for (int t = 0; t < NT; ++t) {
        Vector2d p1 = mesh.node.row(mesh.elem(t, 0));
        Vector2d p2 = mesh.node.row(mesh.elem(t, 1));
        Vector2d p3 = mesh.node.row(mesh.elem(t, 2));
        double area = fem.area(t);
        Le.setZero();
        for (int q = 0; q < nq; ++q) {
            Vector3d lam = quadL.row(q).transpose();
            Vector2d pt = lam(0) * p1 + lam(1) * p2 + lam(2) * p3;
            double fv = f(pt(0), pt(1));
            Le.noalias() += (w(q) * area * fv) * phi_q[q].transpose();
        }
        for (int i = 0; i < locDof; ++i) L(elem2dof(t, i)) += Le(i);
    }
    return L;
}

// ---------------------------------------------------------------------------
// CHIntegrator
// ---------------------------------------------------------------------------
namespace {
// Insert scale * M into a triplet list at block offset (rowOff, colOff).
// Duplicate (row,col) entries are summed by setFromTriplets (used for -(eps^2 A + S Mm)).
void addBlockTo(std::vector<Triplet<double>>& trip, const SparseMatrix<double>& M,
                double scale, int rowOff, int colOff) {
    for (int k = 0; k < M.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
            trip.emplace_back(rowOff + static_cast<int>(it.row()),
                              colOff + static_cast<int>(it.col()),
                              scale * it.value());
}
}

CHIntegrator::CHIntegrator(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                           const SparseMatrix<double>& Mm, const SparseMatrix<double>& A,
                           double tau, double mob, double eps2, double S, int timeOrder)
    : fem_(fem), mesh_(mesh), elem2dof_(elem2dof), Mm_(Mm), A_(A),
      tau_(tau), mob_(mob), eps2_(eps2), S_(S), order_(timeOrder),
      nDof_(static_cast<int>(Mm.rows())), n_(0) {
    // J1: backward-Euler block (used by order 1 and as the SBDF2 bootstrap), (1,1)=(1/tau)Mm.
    buildAndFactor(1.0, J1_, lu1_);
    // J2: SBDF2 block, (1,1)=(3/2/tau)Mm.
    if (order_ == 2) buildAndFactor(1.5, J2_, lu2_);
}

void CHIntegrator::buildAndFactor(double a0, SparseMatrix<double>& J,
                                  SparseLU<SparseMatrix<double>>& lu) {
    int n = nDof_;
    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(3) * Mm_.nonZeros() + static_cast<size_t>(2) * A_.nonZeros());
    addBlockTo(trip, Mm_, a0 / tau_, 0, 0); // (1,1): (a0/tau) Mm
    addBlockTo(trip, A_,  mob_,      0, n); // (1,2): M * A
    addBlockTo(trip, A_,  -eps2_,    n, 0); // (2,1): -(eps^2 A + S Mm)
    addBlockTo(trip, Mm_, -S_,       n, 0);
    addBlockTo(trip, Mm_, 1.0,       n, n); // (2,2): Mm
    J.resize(2 * n, 2 * n);
    J.setFromTriplets(trip.begin(), trip.end());
    J.makeCompressed();
    lu.analyzePattern(J);
    lu.factorize(J);
    if (lu.info() != Success)
        std::cerr << "CHIntegrator: block-matrix factorisation failed\n";
}

void CHIntegrator::setInitial(const VectorXd& c0) {
    c_ = c0;
    cprev_ = c0;
    n_ = 0;
}

bool CHIntegrator::step(const VectorXd& sourceLoad) {
    int n = nDof_;
    VectorXd rhs(2 * n);
    bool useBDF2 = (order_ == 2 && n_ >= 1); // first step bootstraps with backward Euler

    if (useBDF2) {
        VectorXd cstar = 2.0 * c_ - cprev_;                       // 2nd-order extrapolation
        VectorXd b = assembleNonlinearCH(fem_, mesh_, elem2dof_, cstar);
        rhs.head(n) = (1.0 / (2.0 * tau_)) * (Mm_ * (4.0 * c_ - cprev_));
        rhs.tail(n) = b - S_ * (Mm_ * cstar);
    } else {
        VectorXd b = assembleNonlinearCH(fem_, mesh_, elem2dof_, c_);
        VectorXd Mmc = Mm_ * c_;
        rhs.head(n) = (1.0 / tau_) * Mmc;
        rhs.tail(n) = b - S_ * Mmc;
    }
    if (sourceLoad.size() == n) rhs.head(n) += sourceLoad; // forcing enters the c-equation

    VectorXd x;
    if (useBDF2) {
        x = lu2_.solve(rhs);
        if (lu2_.info() != Success) return false;
    } else {
        x = lu1_.solve(rhs);
        if (lu1_.info() != Success) return false;
    }
    cprev_ = c_;
    c_ = x.head(n);
    ++n_;
    return true;
}

// ---------------------------------------------------------------------------
// CHAdaptiveStepper
// ---------------------------------------------------------------------------
CHAdaptiveStepper::CHAdaptiveStepper(FEM& fem, const Mesh& mesh, const MatrixXi& elem2dof,
                                     const SparseMatrix<double>& Mm, const SparseMatrix<double>& A,
                                     double mob, double eps2, double S, double tauMin, double tauMax)
    : fem_(fem), mesh_(mesh), elem2dof_(elem2dof), Mm_(Mm), A_(A),
      mob_(mob), eps2_(eps2), S_(S), tauMin_(tauMin), tauMax_(tauMax),
      nDof_(static_cast<int>(Mm.rows())), lastTau_(tauMin), lastRel_(0.0) {
    // Largest level whose tau = tauMin*2^level does not exceed tauMax (floor, so tau stays in range).
    maxLevel_ = std::max(0, static_cast<int>(std::floor(std::log2(tauMax_ / tauMin_) + 1e-9)));
}

void CHAdaptiveStepper::setInitial(const VectorXd& c0) { c_ = c0; }

int CHAdaptiveStepper::levelFor(double tau) const {
    if (tau <= tauMin_) return 0;
    int lvl = static_cast<int>(std::lround(std::log2(tau / tauMin_)));
    if (lvl < 0) lvl = 0;
    if (lvl > maxLevel_) lvl = maxLevel_;
    return lvl;
}

SparseLU<SparseMatrix<double>>* CHAdaptiveStepper::solverForLevel(int lvl) {
    auto it = cache_.find(lvl);
    if (it != cache_.end()) return &it->second;

    double tau = tauMin_ * std::pow(2.0, lvl);  // backward-Euler block at this tau
    int n = nDof_;
    std::vector<Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(3) * Mm_.nonZeros() + static_cast<size_t>(2) * A_.nonZeros());
    addBlockTo(trip, Mm_, 1.0 / tau, 0, 0);
    addBlockTo(trip, A_,  mob_,      0, n);
    addBlockTo(trip, A_,  -eps2_,    n, 0);
    addBlockTo(trip, Mm_, -S_,       n, 0);
    addBlockTo(trip, Mm_, 1.0,       n, n);
    SparseMatrix<double> J(2 * n, 2 * n);
    J.setFromTriplets(trip.begin(), trip.end());
    J.makeCompressed();

    SparseLU<SparseMatrix<double>>& lu = cache_[lvl];  // default-construct in place (map node is stable)
    lu.analyzePattern(J);
    lu.factorize(J);
    if (lu.info() != Success)
        std::cerr << "CHAdaptiveStepper: factorisation failed at tau-level " << lvl << " (tau=" << tau << ")\n";
    return &lu;
}

bool CHAdaptiveStepper::step(double tauDesired, const VectorXd& src) {
    int lvl = levelFor(tauDesired);
    double tau = tauMin_ * std::pow(2.0, lvl);
    SparseLU<SparseMatrix<double>>* lu = solverForLevel(lvl);

    int n = nDof_;
    VectorXd b = assembleNonlinearCH(fem_, mesh_, elem2dof_, c_); // explicit F'(c^n)
    VectorXd Mmc = Mm_ * c_;
    VectorXd rhs(2 * n);
    rhs.head(n) = (1.0 / tau) * Mmc;
    if (src.size() == n) rhs.head(n) += src;
    rhs.tail(n) = b - S_ * Mmc;

    VectorXd x = lu->solve(rhs);
    if (lu->info() != Success) return false;

    VectorXd cnew = x.head(n);
    VectorXd d = cnew - c_;
    double num = std::sqrt(std::max(0.0, d.dot(Mm_ * d)));
    double den = std::sqrt(std::max(1e-300, cnew.dot(Mm_ * cnew)));
    lastRel_ = num / den;
    lastTau_ = tau;
    c_ = cnew;
    return true;
}
