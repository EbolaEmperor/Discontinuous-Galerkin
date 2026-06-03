#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <string>

#include "Mesh.h"
#include "Argyris.h"
#include "ExactSolution.h"

#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#ifdef EIGEN_USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

using namespace Eigen;
using namespace std;
using namespace std::chrono;

int main(int argc, char** argv) {
    // Argyris (C^1, P5) for  Delta^2 u = f,  clamped BC from the exact solution.
    // Optional CLI:  biharmonic_argyris [Nref] [solver]
    int Nref   = 5;
    double h0  = 0.5;
    int solver = 3;          // 0=LDLT 1=LLT 2=CG 3=CHOLMOD
    if (argc > 1) Nref   = atoi(argv[1]);
    if (argc > 2) solver = atoi(argv[2]);

    ExactSolution sol(0.3);

    cout << "Biharmonic, Argyris C1 element (P5, 21 dof/triangle)\n";
    cout << "Conforming H2 discretization: a(u,v)=∫ D^2u:D^2v, no penalty.\n";
    cout << "Domain: unit square, clamped BC (u, grad u, Hessian, d_n u pinned to exact).\n\n";

    vector<double> hlist(Nref), errL2(Nref), errH1(Nref), errH2(Nref);
    string solverName;

    double h = h0;
    for (int lv = 0; lv < Nref; ++lv) {
        hlist[lv] = h;
        Mesh mesh;
        mesh.getMesh(h);

        auto ta0 = high_resolution_clock::now();
        ArgyrisFEM fem(mesh);
        SparseMatrix<double> A = fem.assembleStiffness();
        VectorXd F = fem.assembleLoad(sol);

        VectorXd c;
        vector<int> freeDof;
        fem.strongBC(sol, c, freeDof);
        F = F - A * c;

        int nFree = static_cast<int>(freeDof.size());
        VectorXi g2f = VectorXi::Constant(fem.nDof, -1);
        for (int i = 0; i < nFree; ++i) g2f(freeDof[i]) = i;

        vector<Triplet<double>> trip;
        trip.reserve(A.nonZeros());
        for (int k = 0; k < A.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                int fr = g2f(it.row()), fc = g2f(it.col());
                if (fr != -1 && fc != -1) trip.emplace_back(fr, fc, it.value());
            }
        SparseMatrix<double> Af(nFree, nFree);
        Af.setFromTriplets(trip.begin(), trip.end());
        VectorXd Ff(nFree);
        for (int i = 0; i < nFree; ++i) Ff(i) = F(freeDof[i]);
        auto ta1 = high_resolution_clock::now();

        auto ts0 = high_resolution_clock::now();
        VectorXd uf;
        if (solver == 0) {
            SimplicialLDLT<SparseMatrix<double>> s; s.compute(Af); uf = s.solve(Ff);
            solverName = "SimplicialLDLT";
        } else if (solver == 1) {
            SimplicialLLT<SparseMatrix<double>> s; s.compute(Af); uf = s.solve(Ff);
            solverName = "SimplicialLLT";
        } else if (solver == 2) {
            ConjugateGradient<SparseMatrix<double>, Lower | Upper> s; s.compute(Af); uf = s.solve(Ff);
            solverName = "ConjugateGradient";
        } else {
#ifdef EIGEN_USE_CHOLMOD
            CholmodSupernodalLLT<SparseMatrix<double>> s; s.compute(Af);
            if (s.info() != Success) { cout << "CHOLMOD factorization failed\n"; return 1; }
            uf = s.solve(Ff);
            solverName = "CholmodSupernodalLLT";
#else
            SimplicialLLT<SparseMatrix<double>> s; s.compute(Af); uf = s.solve(Ff);
            solverName = "SimplicialLLT (CHOLMOD unavailable)";
#endif
        }
        for (int i = 0; i < nFree; ++i) c(freeDof[i]) = uf(i);
        auto ts1 = high_resolution_clock::now();

        auto te0 = high_resolution_clock::now();
        fem.errors(c, sol, errL2[lv], errH1[lv], errH2[lv]);
        auto te1 = high_resolution_clock::now();

        if (lv == 0) cout << "Solver: " << solverName << "\n\n";
        cout << "Lv " << lv << " (h=" << fixed << setprecision(4) << h << ")\n";
        cout << "  nDof=" << fem.nDof << " (free=" << nFree << ")\n";
        cout << setprecision(4)
             << "  time[s]: assemble=" << duration<double>(ta1 - ta0).count()
             << "  solve=" << duration<double>(ts1 - ts0).count()
             << "  error=" << duration<double>(te1 - te0).count() << "\n";
        cout << scientific << setprecision(3)
             << "  L2=" << errL2[lv] << "  H1=" << errH1[lv]
             << "  H2(semi)=" << errH2[lv] << defaultfloat << "\n";
        h /= 2.0;
    }

    cout << "\nConvergence rates (Argyris P5; expect L2~6, H1~5, H2~4):\n";
    for (int i = 0; i < Nref - 1; ++i) {
        double rL2 = log2(errL2[i] / errL2[i + 1]);
        double rH1 = log2(errH1[i] / errH1[i + 1]);
        double rH2 = log2(errH2[i] / errH2[i + 1]);
        cout << "  Lv " << i << "->" << i + 1 << fixed << setprecision(2)
             << "  L2=" << rL2 << "  H1=" << rH1 << "  H2=" << rH2
             << defaultfloat << "\n";
    }
    return 0;
}
