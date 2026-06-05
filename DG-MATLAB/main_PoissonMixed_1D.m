clc
clear
close all

ord = 5;
fem1 = DPk1D(ord);  % FEM for sigma
fem2 = DPk1D(ord);  % FEM for u
fem_star = DPk1D(ord+1);  % FEM for u^* (super-convergence of uh)
alpha = 2 * ord * (ord+1);  % Penalty coefficient

Nref = 6;
h0 = 0.5;
x0 = 0;
x1 = 1;
bdc.left = "N";
bdc.right = "D";
% "D" for Dirichlet BDC, "N" for Neumann BDC.

expscale = 0.5;
u_exact = @(x) sin(pi * x) - exp(x) * expscale;
du_exact = @(x) pi * cos(pi * x) - exp(x) * expscale;
f = @(x) (pi^2) * sin(pi * x) + exp(x) * expscale;

hlist = zeros(1, Nref);
err_u = zeros(1, Nref);
err_ustar = zeros(1, Nref);
err_sigma = zeros(1, Nref);

bdval = [u_exact(0), u_exact(1)];
if bdc.left == "N", bdval(1) = -du_exact(x0); end
if bdc.right == "N", bdval(2) = du_exact(x1); end

for cycle = 1 : Nref
    hlist(cycle) = h0;
    grid = x0 : h0 : x1;
    assert(length(grid) >= 3);

    NT = length(grid) - 1;
    nDof1 = fem1.locDof * NT;
    nDof2 = fem2.locDof * NT;

    %% Solve mixed-form Poisson equation using HDG method.
    M = assembleMass_DPk1D(fem1, grid);
    B1 = assembleMixed_HDG1D(fem1, fem2, grid);
    B2 = assembleMixed_HDG1D(fem2, fem1, grid);
    K = [M, B1; B2, sparse(nDof2,nDof2)];

    P = assembleIP_HDG1D(fem1, fem2, grid, alpha, bdc);
    A = K + P;
    
    F = [zeros(nDof1, 1); 
         assembleLoad_DPk1D(fem2, grid, f)];
    F = F + assembleWeakBDC_HDG1D(fem1, fem2, grid, alpha, bdc, bdval);

    sol = A \ F;
    sigmah = sol(1:nDof1);
    uh = sol(nDof1+1:end);

    err_sigma(cycle) = getL2Err_DPk1D(fem1, grid, sigmah, du_exact);
    err_u(cycle) = getL2Err_DPk1D(fem2, grid, uh, u_exact);

    %% See super-convergence of uh by solving piece-wise Poisson equations.
    locDofstar = fem_star.locDof;
    uh_star = zeros(locDofstar * NT, 1);
    for i = 1 : NT
        loc_sigmah = sigmah((i-1) * fem1.locDof+ (1:fem1.locDof));
        loc_uh = uh((i-1) * fem2.locDof + (1:fem2.locDof));
        loc_uhStar = solveLocalPoisson_HDG1D(fem_star, fem1, fem2, loc_sigmah, loc_uh, grid(i), grid(i+1), f);
        uh_star((i-1) * locDofstar + (1:locDofstar)) = loc_uhStar;
    end

    err_ustar(cycle) = getL2Err_DPk1D(fem_star, grid, uh_star, u_exact);

    h0 = h0 / 2;
end

h = figure;
set(h, "Position", [100, 100, 900, 1000]);

subplot(2, 2, 1);
opt.nPoint = 20;
opt.u_exact = u_exact;
opt.varname = "u^*";
plotSol_DPk1D(fem_star, grid, uh_star, opt);

subplot(2, 2, 2);
opt.nPoint = 20;
opt.u_exact = du_exact;
opt.varname = "\sigma";
plotSol_DPk1D(fem1, grid, sigmah, opt);

subplot(2, 2, 3);
showrateh_mdf(hlist, err_ustar, Nref-1, '-o', "$||u^*-u_h||_{L^2}$");
title("$||u-u^*_h||_{L^2}$", "Interpreter", "latex");

subplot(2, 2, 4);
showrateh_mdf(hlist, err_sigma, Nref-1, '-o', "$||u'-\sigma_h||_{L^2}$");
title("$||u'-\sigma_h||_{L^2}$", "Interpreter", "latex");
