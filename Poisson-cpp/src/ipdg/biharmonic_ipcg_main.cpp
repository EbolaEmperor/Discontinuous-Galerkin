#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <string>

#include "Mesh.h"
#include "FEM.h"
#include "DGAssembly.h"
#include "ExactSolution.h"

#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#ifdef EIGEN_USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

using namespace Eigen;
using namespace std;
using namespace std::chrono;

static MatrixXd regularPolygon(int n) {
    MatrixXd v(n, 2);
    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * i / n;
        v(i, 0) = cos(theta);
        v(i, 1) = sin(theta);
    }
    return v;
}

int main() {
    int ord = 4;
    if (ord < 2) {
        cout << "Biharmonic C0-IPCG requires ord >= 2\n";
        return 1;
    }
    int Nref = 5;
    double h0 = 0.5;
    double sigma = 1.5 * ord * (ord + 1);
    string ip_type = "SIPG"; // options: SIPG, NIPG, IIPG
    double beta = (ip_type == "NIPG") ? -1.0 : (ip_type == "IIPG") ? 0.0 : 1.0;

    int solver_type = 3; // 0=LDLT, 1=LLT, 2=CG, 3=Cholmod (if available)

    ExactSolution sol(0.3);

    cout << "Biharmonic C0-IPCG (C++)\n";
    cout << "ord = " << ord << " (sigma=" << sigma << ", beta=" << beta << ")\n";
    cout << "IP type = " << ip_type << endl;
    if (solver_type == 0)
        cout << "Solver: SimplicialLDLT\n";
    else if (solver_type == 1)
        cout << "Solver: SimplicialLLT\n";
    else if (solver_type == 2)
        cout << "Solver: ConjugateGradient\n";
    else if (solver_type == 3) {
#ifdef EIGEN_USE_CHOLMOD
        cout << "Solver: CholmodSupernodalLLT\n";
#else
        cout << "Solver: Cholmod requested but not found. Fallback to SimplicialLLT.\n";
        solver_type = 1;
#endif
    }
    cout << fixed << setprecision(4);

    struct CaseCfg {
        string name;
        bool isPolygon;
        MatrixXd vertices;
        double h0;
    };

    vector<CaseCfg> cases;
    cases.push_back({"Unit square", false, MatrixXd(), h0});
    // cases.push_back({"Polygon (regular hexagon)", true, regularPolygon(6), h0});

    for (const auto& cfg : cases) {
        cout << "\nDomain: " << cfg.name << endl;
        vector<double> hlist(Nref);
        vector<double> errL2(Nref);
        vector<double> errH1(Nref);

        double h = cfg.h0;
        for (int lv = 0; lv < Nref; ++lv) {
            hlist[lv] = h;
            cout << "Lv " << lv << " (h=" << h << "):" << endl;

            auto t_init_start = high_resolution_clock::now();
            Mesh mesh;
            if (cfg.isPolygon) {
                mesh.getPolygonMesh(cfg.vertices, h);
            } else {
                mesh.getMesh(h);
            }

            FEM fem(ord, mesh);
            MatrixXi elem2dof;
            MatrixXd dofCoords;
            int nDof = 0;
            fem.getConformingDOF(mesh, elem2dof, nDof, dofCoords);

            MatrixXi edge, edge2side;
            mesh.getEdge2Side(edge, edge2side);
            auto t_init_end = high_resolution_clock::now();

            auto t_assemble_start = high_resolution_clock::now();
            SparseMatrix<double> K = assembleK_Bihar2D(fem, mesh, elem2dof);
            SparseMatrix<double> P = assembleIP_Bihar2D(fem, mesh, elem2dof, edge, edge2side, sigma, beta);
            SparseMatrix<double> A = K + P;

            VectorXd F = assembleLoadVectorBihar(fem, mesh, elem2dof, sol);
            F += modifyLoadVectorNitsche_Bihar2D(fem, mesh, elem2dof, edge, edge2side, sigma, beta, sol);

            VectorXd c;
            vector<int> freeDof;
            const MatrixXd* polyPtr = cfg.isPolygon ? &cfg.vertices : nullptr;
            interpStrongBDC(fem, mesh, elem2dof, sol, c, freeDof, polyPtr);
            F = F - A * c;
            auto t_assemble_end = high_resolution_clock::now();

            auto t_solve_start = high_resolution_clock::now();
            int nFree = static_cast<int>(freeDof.size());
            SparseMatrix<double> A_free(nFree, nFree);
            VectorXd F_free(nFree);

            VectorXi g2f = VectorXi::Constant(nDof, -1);
            for (int i = 0; i < nFree; ++i) g2f(freeDof[i]) = i;

            vector<Triplet<double>> trip;
            trip.reserve(A.nonZeros());
            for (int k = 0; k < A.outerSize(); ++k) {
                for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                    int fr = g2f(it.row());
                    int fc = g2f(it.col());
                    if (fr != -1 && fc != -1) {
                        trip.emplace_back(fr, fc, it.value());
                    }
                }
            }
            A_free.setFromTriplets(trip.begin(), trip.end());
            for (int i = 0; i < nFree; ++i) F_free(i) = F(freeDof[i]);

            VectorXd u_free;
            if (solver_type == 0) {
                SimplicialLDLT<SparseMatrix<double>> solver;
                solver.compute(A_free);
                if (solver.info() != Success) { cout << "Decomposition failed\n"; return 1; }
                u_free = solver.solve(F_free);
                if (solver.info() != Success) { cout << "Solve failed\n"; return 1; }
            } else if (solver_type == 1) {
                SimplicialLLT<SparseMatrix<double>> solver;
                solver.compute(A_free);
                if (solver.info() != Success) { cout << "Decomposition failed\n"; return 1; }
                u_free = solver.solve(F_free);
                if (solver.info() != Success) { cout << "Solve failed\n"; return 1; }
            } else if (solver_type == 2) {
                ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
                solver.compute(A_free);
                if (solver.info() != Success) { cout << "Decomposition failed\n"; return 1; }
                u_free = solver.solve(F_free);
                if (solver.info() != Success) { cout << "Solve failed\n"; return 1; }
            } else if (solver_type == 3) {
#ifdef EIGEN_USE_CHOLMOD
                CholmodSupernodalLLT<SparseMatrix<double>> solver;
                solver.compute(A_free);
                if (solver.info() != Success) { cout << "Decomposition failed\n"; return 1; }
                u_free = solver.solve(F_free);
                if (solver.info() != Success) { cout << "Solve failed\n"; return 1; }
#endif
            }

            for (int i = 0; i < nFree; ++i) c(freeDof[i]) = u_free(i);
            auto t_solve_end = high_resolution_clock::now();

            auto t_error_start = high_resolution_clock::now();
            getH1Err(fem, mesh, elem2dof, c, sol, errH1[lv], errL2[lv]);
            auto t_error_end = high_resolution_clock::now();

            duration<double> d_init = t_init_end - t_init_start;
            duration<double> d_assemble = t_assemble_end - t_assemble_start;
            duration<double> d_solve = t_solve_end - t_solve_start;
            duration<double> d_error = t_error_end - t_error_start;

            cout << "  Init: " << d_init.count() << "s\n";
            cout << "  Assemble: " << d_assemble.count() << "s\n";
            cout << "  Solve: " << d_solve.count() << "s\n";
            cout << "  Error: " << d_error.count() << "s\n";
            cout << "  nDof=" << nDof << " L2=" << scientific << errL2[lv]
                 << " H1=" << errH1[lv] << defaultfloat << endl;

            h /= 2.0;
        }

        cout << "\nRates (" << cfg.name << "):" << endl;
        for (int i = 0; i < Nref - 1; ++i) {
            double rL2 = log2(errL2[i] / errL2[i + 1]);
            double rH1 = log2(errH1[i] / errH1[i + 1]);
            cout << "Lv " << i << "->" << i + 1 << ": L2 rate=" << fixed << setprecision(2)
                 << rL2 << " H1 rate=" << rH1 << defaultfloat << endl;
        }
    }

    return 0;
}
