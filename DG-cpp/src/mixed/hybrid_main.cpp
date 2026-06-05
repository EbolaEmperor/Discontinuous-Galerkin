#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>

#include "Mesh.h"
#include "FEM.h"
#include "VecFEM.h"
#include "AssemblyHDG.h"   // l2ErrorScalar / l2ErrorVector
#include "PostProcess.h"   // solveLocalPoisson (superconvergent u*)
#include "HybridHDG.h"
#include "ExactSolution.h"

using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::VectorXd;

int main(int argc, char** argv) {
    // ----------------------- configuration -----------------------
    // Optional CLI overrides (all have sensible defaults):
    //   poisson_hdg_hybrid [k] [tau] [Nref] [solver]
    int    ord    = 2;     // polynomial degree k (u, sigma both in P_k)
    double tau    = 1.0;   // HDG stabilization (O(1) -> sigma,u order k+1)
    double h0     = 0.5;   // initial mesh size
    int    Nref   = 5;     // refinement levels
    int    solver = 3;     // 0=LDLT 1=LLT 2=CG 3=CHOLMOD
    ExactSolution sol(0.3);
    if (argc > 1) ord    = std::atoi(argv[1]);
    if (argc > 2) tau    = std::atof(argv[2]);
    if (argc > 3) Nref   = std::atoi(argv[3]);
    if (argc > 4) solver = std::atoi(argv[4]);
    // --------------------------------------------------------------

    std::cout << "Hybridized HDG Poisson  (sigma=grad u, -div sigma = f)\n";
    std::cout << "ord(k)=" << ord << "  tau=" << tau << "\n";
    std::cout << "Global system: SPD Schur complement in the trace lambda only.\n\n";

    std::vector<double> errSigma(Nref), errU(Nref), errUStar(Nref);
    std::vector<int>    ndof(Nref);

    double h = h0;
    for (int lv = 0; lv < Nref; ++lv) {
        Mesh mesh;
        mesh.getMesh(h);
        MatrixXi edge, edge2side;
        mesh.getEdge2Side(edge, edge2side);

        VecFEM femQ(ord, mesh);
        FEM    femU(ord, mesh);
        MatrixXi elem2dofSigma, elem2dofU;
        int nSig, nU;
        femQ.getDOF(mesh, elem2dofSigma, nSig);
        femU.getDOF(mesh, elem2dofU, nU);

        HybridHDGResult r = solveHybridHDG(femQ, femU, mesh, elem2dofSigma,
                                           elem2dofU, edge, edge2side, tau, sol, solver);

        errSigma[lv] = l2ErrorVector(femQ, mesh, elem2dofSigma, r.sigmah, sol);
        errU[lv]     = l2ErrorScalar(femU, mesh, elem2dofU, r.uh, sol);

        // ---- local postprocessing -> superconvergent u* in P_{k+1} ----
        FEM femStar(ord + 1, mesh);
        MatrixXi elem2dofStar;
        int nStar;
        femStar.getDOF(mesh, elem2dofStar, nStar);
        auto f = [&sol](const Vector2d& p) {
            MatrixXd pt(1, 2); pt << p(0), p(1);
            return -sol.laplace_u_exact(pt)(0);
        };
        auto tp0 = std::chrono::high_resolution_clock::now();
        VectorXd ustar = solveLocalPoisson(femStar, femQ, femU, mesh, elem2dofStar,
                                           elem2dofSigma, elem2dofU, r.sigmah, r.uh, f);
        auto tp1 = std::chrono::high_resolution_clock::now();
        double tPost = std::chrono::duration<double>(tp1 - tp0).count();
        errUStar[lv] = l2ErrorScalar(femStar, mesh, elem2dofStar, ustar, sol);

        // assemble (condensation) + recover are both part of the HDG machinery
        double tAsm   = r.tAssemble + r.tRecover;
        double tTotal = tAsm + r.tSolve + tPost;

        ndof[lv] = r.nLambdaFree;
        if (lv == 0) std::cout << "Solver: " << r.solverName << "\n\n";
        std::cout << "Lv " << lv << " (h=" << std::fixed << std::setprecision(4) << h << ")\n";
        std::cout << "  trace dofs: " << r.nLambdaFree << " free / " << r.nLambda
                  << " total   (nSig=" << nSig << ", nU=" << nU << ")\n";
        std::cout << std::setprecision(4)
                  << "  time[s]: assemble=" << tAsm << " (cond=" << r.tAssemble
                  << "+recover=" << r.tRecover << ")  solve=" << r.tSolve
                  << "  post=" << tPost << "  TOTAL=" << tTotal << "\n";
        std::cout << std::scientific << std::setprecision(3)
                  << "  |sigma-grad u|_L2 = " << errSigma[lv] << "\n"
                  << "  |u - u_h|_L2      = " << errU[lv] << "\n"
                  << "  |u - u*|_L2       = " << errUStar[lv] << "\n"
                  << std::defaultfloat;
        h /= 2.0;
    }

    std::cout << "\nConvergence rates (k=" << ord << "):\n";
    std::cout << "  expect: sigma ~ k+1=" << ord + 1
              << ", u ~ k+1=" << ord + 1
              << ", u* ~ k+2=" << ord + 2 << "\n";
    for (int i = 0; i < Nref - 1; ++i) {
        double rS  = std::log2(errSigma[i]  / errSigma[i + 1]);
        double rU  = std::log2(errU[i]      / errU[i + 1]);
        double rUs = std::log2(errUStar[i]  / errUStar[i + 1]);
        std::cout << "  Lv " << i << "->" << i + 1 << std::fixed << std::setprecision(2)
                  << "  sigma=" << rS << "  u=" << rU << "  u*=" << rUs
                  << std::defaultfloat << "\n";
    }
    return 0;
}
