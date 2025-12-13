#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "Mesh.h"
#include "FEM.h"
#include "DGAssembly.h"
#include "ExactSolution.h"

// Includes for Solvers
#include <Eigen/SparseCholesky> // For SimplicialLLT
#include <Eigen/IterativeLinearSolvers> // For ConjugateGradient
#ifdef EIGEN_USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

using namespace Eigen;
using namespace std;
using namespace std::chrono;

int main() {
    int ord = 4;
    double h0 = 0.5;
    int Nref = 6;
    double sigma = 3 * ord * (ord + 1);
    double beta = 1; // SIPG
    
    // Choose solver: 0 = SimplicialLDLT, 1 = SimplicialLLT, 2 = ConjugateGradient, 3 = Cholmod (if available)
    int solver_type = 3; 

    ExactSolution sol(0.3);
    
    vector<double> hlist(Nref);
    vector<double> errL2(Nref);
    vector<double> errH1(Nref);
    
    cout << "Poisson DG (C++)" << endl;
    cout << "ord = " << ord << endl;
    cout << "sigma = " << sigma << ", beta = " << beta << endl;
    if(solver_type == 0) cout << "Solver: SimplicialLDLT" << endl;
    else if(solver_type == 1) cout << "Solver: SimplicialLLT (Cholesky)" << endl;
    else if(solver_type == 2) cout << "Solver: ConjugateGradient" << endl;
    else if(solver_type == 3) {
        #ifdef EIGEN_USE_CHOLMOD
        cout << "Solver: CholmodSupernodalLLT" << endl;
        #else
        cout << "Solver: Cholmod requested but not found. Fallback to SimplicialLLT." << endl;
        solver_type = 1;
        #endif
    }
    cout << fixed << setprecision(4);
    
    for (int lv = 0; lv < Nref; ++lv) {
        hlist[lv] = h0;
        
        cout << "Lv " << lv << " (h=" << h0 << "):" << endl;

        // 1. Init
        auto t_init_start = high_resolution_clock::now();
        Mesh mesh;
        mesh.getMesh(h0);
        
        FEM fem(ord, mesh);
        
        MatrixXi elem2dof;
        int nDof;
        fem.getDOF(mesh, elem2dof, nDof);
        
        MatrixXi edge, edge2side;
        mesh.getEdge2Side(edge, edge2side);
        auto t_init_end = high_resolution_clock::now();
        
        // 2. Assemble
        auto t_assemble_start = high_resolution_clock::now();
        SparseMatrix<double> K = assembleK_Poi2D(fem, mesh, elem2dof);
        SparseMatrix<double> P = assembleIP_Poi2D(fem, mesh, elem2dof, edge, edge2side, sigma, beta);
        
        SparseMatrix<double> A = K + P;
        
        VectorXd F = assembleLoadVector(fem, mesh, elem2dof, sol);
        
        VectorXd c;
        vector<int> freeDof;
        interpStrongBDC(fem, mesh, elem2dof, sol, c, freeDof);
        
        F = F - A * c;
        auto t_assemble_end = high_resolution_clock::now();
        
        // 3. Solve
        auto t_solve_start = high_resolution_clock::now();
        
        // Extract free system
        int nFree = freeDof.size();
        SparseMatrix<double> A_free(nFree, nFree);
        VectorXd F_free(nFree);
        
        // Map global to free
        VectorXi g2f = VectorXi::Constant(nDof, -1);
        for(int i=0; i<nFree; ++i) g2f(freeDof[i]) = i;
        
        vector<Triplet<double>> trip;
        // Iterate A efficiently
        for (int k=0; k<A.outerSize(); ++k) {
            for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it) {
                int r = it.row();
                int c_idx = it.col();
                int fr = g2f(r);
                int fc = g2f(c_idx);
                if (fr != -1 && fc != -1) {
                    trip.push_back(Triplet<double>(fr, fc, it.value()));
                }
            }
        }
        A_free.setFromTriplets(trip.begin(), trip.end());
        
        for(int i=0; i<nFree; ++i) F_free(i) = F(freeDof[i]);
        
        // Solve
        VectorXd u_free;
        if (solver_type == 0) {
            SimplicialLDLT<SparseMatrix<double>> solver;
            solver.compute(A_free);
            if(solver.info() != Success) { cout << "Decomposition failed" << endl; return 1; }
            u_free = solver.solve(F_free);
            if(solver.info() != Success) { cout << "Solve failed" << endl; return 1; }
        } else if (solver_type == 1) {
            SimplicialLLT<SparseMatrix<double>> solver; // Cholesky
            solver.compute(A_free);
            if(solver.info() != Success) { cout << "Decomposition failed" << endl; return 1; }
            u_free = solver.solve(F_free);
            if(solver.info() != Success) { cout << "Solve failed" << endl; return 1; }
        } else if (solver_type == 2) {
            ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
            solver.compute(A_free);
            if(solver.info() != Success) { cout << "Decomposition failed" << endl; return 1; }
            u_free = solver.solve(F_free);
            if(solver.info() != Success) { cout << "Solve failed" << endl; return 1; }
        } else if (solver_type == 3) {
            #ifdef EIGEN_USE_CHOLMOD
            CholmodSupernodalLLT<SparseMatrix<double>> solver;
            solver.compute(A_free);
            if(solver.info() != Success) { cout << "Decomposition failed" << endl; return 1; }
            u_free = solver.solve(F_free);
            if(solver.info() != Success) { cout << "Solve failed" << endl; return 1; }
            #endif
        }

        for(int i=0; i<nFree; ++i) c(freeDof[i]) = u_free(i);
        auto t_solve_end = high_resolution_clock::now();
        
        // 4. Error
        auto t_error_start = high_resolution_clock::now();
        getH1Err(fem, mesh, elem2dof, c, sol, errH1[lv], errL2[lv]);
        auto t_error_end = high_resolution_clock::now();
        
        // Print times
        duration<double> d_init = t_init_end - t_init_start;
        duration<double> d_assemble = t_assemble_end - t_assemble_start;
        duration<double> d_solve = t_solve_end - t_solve_start;
        duration<double> d_error = t_error_end - t_error_start;
        
        cout << "  Init: " << d_init.count() << "s" << endl;
        cout << "  Assemble: " << d_assemble.count() << "s" << endl;
        cout << "  Solve: " << d_solve.count() << "s" << endl;
        cout << "  Error: " << d_error.count() << "s" << endl;
        
        cout << "  nDof=" << nDof << " L2=" << scientific << errL2[lv] << " H1=" << errH1[lv] << defaultfloat << endl;
             
        h0 /= 2.0;
    }
    
    // Compute rates
    cout << "\nRates:" << endl;
    for(int i=0; i<Nref-1; ++i) {
        double rL2 = log2(errL2[i] / errL2[i+1]);
        double rH1 = log2(errH1[i] / errH1[i+1]);
        cout << "Lv " << i << "-" << i+1 << ": L2 rate=" << fixed << setprecision(2) << rL2 
             << " H1 rate=" << rH1 << endl;
    }

    return 0;
}
