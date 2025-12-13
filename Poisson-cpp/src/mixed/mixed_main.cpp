#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#ifdef EIGEN_USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

#include "Mesh.h"
#include "FEM.h"
#include "VecFEM.h"
#include "AssemblyHDG.h"
#include "PostProcess.h"
#include "DGAssembly.h"
#include "ExactSolution.h"

using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::Vector2d;
using Eigen::VectorXd;

int main() {
    int ord = 4;
    double alpha = 1.0;
    double h0 = 0.5;
    int Nref = 5;
    int solver_type = 3; // 0=LDLT,1=LLT,2=CG,3=Cholmod(if available)

    ExactSolution sol(0.3);

    std::vector<double> hlist(Nref);
    std::vector<double> errSigma(Nref), errU(Nref), errUStar(Nref);

    std::cout << "Mixed HDG Poisson (C++)\n";
    std::cout << "ord=" << ord << " alpha=" << alpha << "\n";

    for (int lv = 0; lv < Nref; ++lv) {
        double h = h0;
        hlist[lv] = h;
        std::cout << "Lv " << lv << " (h=" << h << ")\n";

        Mesh mesh;
        mesh.getMesh(h);
        MatrixXi edge, edge2side;
        mesh.getEdge2Side(edge, edge2side);

        VecFEM femSigma(ord, mesh);
        FEM femU(ord, mesh);
        MatrixXi elem2dofSigma, elem2dofU;
        int nSig, nU;
        femSigma.getDOF(mesh, elem2dofSigma, nSig);
        femU.getDOF(mesh, elem2dofU, nU);
        int nTot = nSig + nU;

        SparseMatrix<double> M = assembleMass(femSigma, mesh, elem2dofSigma);
        SparseMatrix<double> B1 = assembleDivMass(femSigma, femU, mesh, elem2dofSigma, elem2dofU);
        SparseMatrix<double> B2 = assembleGradMass(femU, femSigma, mesh, elem2dofU, elem2dofSigma);

        std::vector<Triplet<double>> trip;
        trip.reserve(M.nonZeros() + B1.nonZeros() + B2.nonZeros());
        auto addBlock = [&](const SparseMatrix<double>& mat, int roff, int coff) {
            for (int k = 0; k < mat.outerSize(); ++k) {
                for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                    trip.emplace_back(it.row() + roff, it.col() + coff, it.value());
                }
            }
        };
        addBlock(M, 0, 0);
        addBlock(B1, 0, nSig);
        addBlock(B2, nSig, 0);

        SparseMatrix<double> K(nTot, nTot);
        K.setFromTriplets(trip.begin(), trip.end());

        SparseMatrix<double> P = assembleIP_HDG(femSigma, femU, mesh, elem2dofSigma, elem2dofU, edge, edge2side, alpha);
        SparseMatrix<double> A = K + P;

        VectorXd F = VectorXd::Zero(nTot);
        F.tail(nU) = assembleLoadVector(femU, mesh, elem2dofU, sol);
        F += assembleWeakBDC(femSigma, femU, mesh, elem2dofSigma, elem2dofU, edge, edge2side, alpha, sol);

        VectorXd solVec;
        if (solver_type == 0) {
            Eigen::SimplicialLDLT<SparseMatrix<double>> solver;
            solver.compute(A);
            solVec = solver.solve(F);
        } else if (solver_type == 1) {
            Eigen::SimplicialLLT<SparseMatrix<double>> solver;
            solver.compute(A);
            solVec = solver.solve(F);
        } else if (solver_type == 2) {
            Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
            solver.compute(A);
            solVec = solver.solve(F);
        } else {
#ifdef EIGEN_USE_CHOLMOD
            Eigen::CholmodSupernodalLLT<SparseMatrix<double>> solver;
            solver.compute(A);
            solVec = solver.solve(F);
#else
            Eigen::SimplicialLLT<SparseMatrix<double>> solver;
            solver.compute(A);
            solVec = solver.solve(F);
#endif
        }

        VectorXd sigmah = solVec.head(nSig);
        VectorXd uh = solVec.tail(nU);

        errSigma[lv] = l2ErrorVector(femSigma, mesh, elem2dofSigma, sigmah, sol);
        errU[lv] = l2ErrorScalar(femU, mesh, elem2dofU, uh, sol);

        FEM femStar(ord + 1, mesh);
        MatrixXi elem2dofStar;
        int nStar;
        femStar.getDOF(mesh, elem2dofStar, nStar);
        auto f = [&sol](const Vector2d& p) {
            MatrixXd pt(1, 2);
            pt << p(0), p(1);
            return -sol.laplace_u_exact(pt)(0);
        };
        VectorXd ustar = solveLocalPoisson(femStar, femSigma, femU, mesh, elem2dofStar, elem2dofSigma, elem2dofU, sigmah, uh, f);
        errUStar[lv] = l2ErrorScalar(femStar, mesh, elem2dofStar, ustar, sol);

        std::cout << "  nSig=" << nSig << " nU=" << nU;
        std::cout << " |sigma|_L2=" << std::scientific << errSigma[lv]
                  << " |u|_L2=" << errU[lv]
                  << " |u*|_L2=" << errUStar[lv] << std::defaultfloat << "\n";

        h0 /= 2.0;
    }

    std::cout << "\nRates:\n";
    for (int i = 0; i < Nref - 1; ++i) {
        double rSigma = std::log2(errSigma[i] / errSigma[i + 1]);
        double rU = std::log2(errU[i] / errU[i + 1]);
        double rUStar = std::log2(errUStar[i] / errUStar[i + 1]);
        std::cout << "Lv " << i << "->" << i + 1
                  << " sigma " << std::fixed << std::setprecision(2) << rSigma
                  << " u " << rU << " u* " << rUStar << std::defaultfloat << "\n";
    }

    return 0;
}
