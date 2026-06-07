// TEMPORARY diagnostic: temporal self-convergence (Richardson) of the DG-NS solver
// under the ACTUAL cylinder-flow boundary conditions (inflow Dirichlet vel + KIO
// Neumann p; no-slip cylinder + KIO Neumann p; slip walls + dp/dn=0; do-nothing
// outflow + Dirichlet p=0).  No exact solution exists, so we measure the temporal
// error as ||u_h(dt) - u_h(dt_ref)||_L2 on a FIXED mesh (spatial error cancels).
// A coarse mesh + short time suffice: temporal order is a property of the scheme +
// BC layout, not of the resolution.
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "FEM.h"
#include "Mesh.h"
#include "MeshGen.h"
#include "NavierStokes.h"

using namespace Eigen;
using namespace std;
using namespace ns;

struct Ctx {
    FEM* fem; Mesh* mesh; const MatrixXi* e2d; const MatrixXi* edge; const MatrixXi* e2s;
    VectorXi tag, isCyl; double nu, Uinf, sigma, gradDiv, ppeDivDamping, rampTime; int pressureMode; VectorXd u0, v0;
    SparseMatrix<double> M, Gx, Gy; SimplicialLDLT<SparseMatrix<double>> luM;
    double setupSec = 0.0, stepSec = 0.0;
    long long steps = 0, cases = 0;
};

static BCData makeBC(const Ctx& c) {
    int NE = c.edge->rows();
    BCData bc;
    bc.bcU = VectorXi::Zero(NE); bc.bcV = VectorXi::Zero(NE); bc.bcP = VectorXi::Zero(NE);
    for (int e = 0; e < NE; ++e) {
        switch (c.tag(e)) {
            case BD_INFLOW:  bc.bcU(e)=BCN_DIRICHLET; bc.bcV(e)=BCN_DIRICHLET; bc.bcP(e)=BCN_NEUMANN_HO;   break;
            case BD_CYL:     bc.bcU(e)=BCN_DIRICHLET; bc.bcV(e)=BCN_DIRICHLET; bc.bcP(e)=BCN_NEUMANN_HO;    break;
            case BD_WALL:    bc.bcU(e)=BCN_NEUMANN;   bc.bcV(e)=BCN_DIRICHLET; bc.bcP(e)=BCN_NEUMANN_ZERO;  break;
            case BD_OUTFLOW: bc.bcU(e)=BCN_NEUMANN;   bc.bcV(e)=BCN_NEUMANN;   bc.bcP(e)=BCN_NEUMANN_HO;    break;
            default: break;
        }
    }
    auto ramp = [](double t, double tr) {
        if (tr <= 0.0 || t >= tr) return 1.0;
        double s = std::max(0.0, t / tr);
        return s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
    };
    auto rampDer = [](double t, double tr) {
        if (tr <= 0.0 || t <= 0.0 || t >= tr) return 0.0;
        double s = t / tr;
        return (30.0 * s * s - 60.0 * s * s * s + 30.0 * s * s * s * s) / tr;
    };
    double U = c.Uinf, tr = c.rampTime;
    bc.velDir = [U, tr, ramp](double x, double y, double t) -> Vector2d {
        if (std::hypot(x, y) < 0.5 + 0.3 * 0.3) return Vector2d(0.0, 0.0);  // no-slip cylinder (r=0.5)
        return Vector2d(U * ramp(t, tr), 0.0);
    };
    bc.velAcc  = [U, tr, rampDer](double, double, double t) {
        return Vector2d(U * rampDer(t, tr), 0.0);
    };
    bc.presDir = [](double, double, double) { return 0.0; };
    return bc;
}

static double velL2Diff(const Ctx& c, const VectorXd& ua, const VectorXd& va,
                        const VectorXd& ub, const VectorXd& vb) {
    int NT = c.mesh->elem.rows(), locDof = c.fem->locDof;
    MatrixXd quadL; VectorXd w; c.fem->quad2d(quadL, w);
    int nq = (int)w.size();
    std::vector<RowVectorXd> phi(nq);
    for (int q = 0; q < nq; ++q) phi[q] = c.fem->computeBasisValue_all(quadL.row(q)).row(0);
    double err = 0.0; VectorXd du(locDof), dv(locDof);
    for (int e = 0; e < NT; ++e) {
        for (int i = 0; i < locDof; ++i) { du(i) = ua((*c.e2d)(e,i)) - ub((*c.e2d)(e,i));
                                           dv(i) = va((*c.e2d)(e,i)) - vb((*c.e2d)(e,i)); }
        double area = c.fem->area(e);
        for (int q = 0; q < nq; ++q) {
            double a = phi[q].dot(du), b = phi[q].dot(dv);
            err += w(q) * area * (a*a + b*b);
        }
    }
    return std::sqrt(err);
}

static double divL2(const Ctx& c, const VectorXd& u, const VectorXd& v) {
    VectorXd d = c.luM.solve(c.Gx * u + c.Gy * v);
    return std::sqrt(d.dot(c.M * d));
}

static void runCase(Ctx& c, double dt, double T, VectorXd& uOut, VectorXd& vOut) {
    BCData bc = makeBC(c);
    auto tSetup0 = std::chrono::steady_clock::now();
    NSIntegrator integ(*c.fem, *c.mesh, *c.e2d, *c.edge, *c.e2s, bc, c.nu, dt, c.sigma,
                       1.0, c.gradDiv, c.pressureMode, c.ppeDivDamping);
    auto tSetup1 = std::chrono::steady_clock::now();
    integ.setInitial(c.u0, c.v0);
    int nsteps = (int)std::lround(T / dt);
    auto tStep0 = std::chrono::steady_clock::now();
    for (int s = 1; s <= nsteps; ++s) if (!integ.step(s * dt)) { std::cerr << "  solve fail @"<<s<<"\n"; uOut=vOut=VectorXd(); return; }
    auto tStep1 = std::chrono::steady_clock::now();
    c.setupSec += std::chrono::duration<double>(tSetup1 - tSetup0).count();
    c.stepSec += std::chrono::duration<double>(tStep1 - tStep0).count();
    c.steps += nsteps;
    ++c.cases;
    uOut = integ.u(); vOut = integ.v();
}

int main(int argc, char** argv) {
    cout << "DG-NS TEMPORAL order under the CYLINDER-FLOW boundary conditions\n";
    // Temporal order is Re-independent (set by the BC/operator structure), so we use a
    // mild Re where the coarse body-fitted mesh is well-resolved and stable.
    int ord = 2; double Re = 20.0, Uinf = 1.0, D = 1.0, nu = Uinf * D / Re;   // nu=0.05

    // moderate/coarse body-fitted mesh (order is resolution-independent)
    CylinderGeom geom{-3.0, 8.0, -3.0, 3.0, 0.0, 0.0, 0.5};
    Mesh mesh; generateCylinderMesh(mesh, geom, 0.16, 4.0, 0.30, 800, false);
    FEM fem(ord, mesh);
    MatrixXi e2d; int nDof; fem.getDOF(mesh, e2d, nDof);
    MatrixXi edge, e2s; mesh.getEdge2Side(edge, e2s);

    Ctx c; c.fem=&fem; c.mesh=&mesh; c.e2d=&e2d; c.edge=&edge; c.e2s=&e2s;
    c.nu=nu; c.Uinf=Uinf; c.sigma = 8.0*(ord+1)*(ord+1);
    c.gradDiv = (argc > 1) ? std::stod(argv[1]) : 0.0;
    c.ppeDivDamping = (argc > 3) ? std::stod(argv[3]) : 10.0;
    c.rampTime = 0.5;
    c.pressureMode = NSPRESSURE_DIRECT_PPE;
    c.tag = classifyEdges(mesh, edge, e2s, geom);
    c.M = assembleScalarMassDG(fem, mesh, e2d);
    assembleWeakGrad(fem, mesh, e2d, edge, e2s, c.Gx, c.Gy);
    c.luM.compute(c.M);
    cout << "  nElem=" << mesh.elem.rows() << "  nDof=" << nDof << "  Re=" << Re << "  nu=" << nu << "  dP" << ord
         << "  pressure=DirectPPE  gradDiv=" << c.gradDiv << "  ppeDivDamp=" << c.ppeDivDamping
         << "  rampTime=" << c.rampTime << "\n";

    // initial condition: uniform free stream + asymmetric wake kick (as in ns_main)
    auto project = [&](const std::function<double(double,double)>& f) -> VectorXd {
        int NT = mesh.elem.rows(), locDof = fem.locDof;
        MatrixXd quadL; VectorXd w; fem.quad2d(quadL, w); int nq=(int)w.size();
        std::vector<RowVectorXd> phi(nq);
        for (int q=0;q<nq;++q) phi[q]=fem.computeBasisValue_all(quadL.row(q)).row(0);
        VectorXd b = VectorXd::Zero(nDof);
        for (int e=0;e<NT;++e){ Vector2d p1=mesh.node.row(mesh.elem(e,0)),p2=mesh.node.row(mesh.elem(e,1)),p3=mesh.node.row(mesh.elem(e,2));
            double area=fem.area(e); VectorXd le=VectorXd::Zero(locDof);
            for(int q=0;q<nq;++q){Vector3d lam=quadL.row(q).transpose();Vector2d P=lam(0)*p1+lam(1)*p2+lam(2)*p3;
                le.noalias()+=(w(q)*area*f(P.x(),P.y()))*phi[q].transpose();}
            for(int i=0;i<locDof;++i) b(e2d(e,i))+=le(i);} return c.luM.solve(b); };
    c.u0 = project([&](double,double){ return 0.0; });
    c.v0 = project([&](double,double){ return 0.0; });

    double T = 1.5;
    int Nref = (argc > 2) ? std::stoi(argv[2]) : 19200;
    double dtRef = T / Nref;
    VectorXd uRef, vRef; runCase(c, dtRef, T, uRef, vRef);
    cout << "  T=" << T << "  dt_ref=" << scientific << setprecision(2) << dtRef
         << "  ||div u_ref||=" << divL2(c, uRef, vRef) << "\n\n";
    cout << "   dt          ||u(dt)-u_ref||    rate     ||div u(dt)||\n";

    std::vector<double> dts;
    for (int k = 0; k < 6; ++k) dts.push_back(T / (100.0 * std::pow(2.0, k)));  // T/100 .. T/3200
    double prev=-1, prevDt=0;
    for (double dt : dts) {
        VectorXd uu, vv; runCase(c, dt, T, uu, vv);
        if (uu.size()==0) { cout << "  " << scientific<<setprecision(2)<<dt << "   (unstable/failed)\n"; continue; }
        double d = velL2Diff(c, uu, vv, uRef, vRef);
        double dv = divL2(c, uu, vv);
        cout << "  " << scientific << setprecision(2) << dt << "    " << setprecision(3) << d;
        if (prev>0) cout << "        " << fixed << setprecision(2) << std::log(prev/d)/std::log(prevDt/dt);
        else        cout << "            ";
        cout << "      " << scientific << setprecision(3) << dv << "\n";
        prev=d; prevDt=dt;
    }
    cout << "\n  timing: setup=" << fixed << setprecision(2) << c.setupSec << "s over " << c.cases
         << " cases, step-loop=" << c.stepSec << "s over " << c.steps
         << " steps, avg step=" << scientific << setprecision(3)
         << (c.steps ? c.stepSec / (double)c.steps : 0.0) << "s\n";
    return 0;
}
