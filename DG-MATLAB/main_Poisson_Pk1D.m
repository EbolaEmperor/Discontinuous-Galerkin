%% Lagrange Pk element for Poisson equation

clc
clear
close all

ord = 4;
Nref = 7;
h0 = 0.5;
expscale = 1/2;
u_exact = @(x) sin(pi * x) - exp(x) * expscale;
du_exact = @(x) pi * cos(pi * x) - exp(x) * expscale;
f = @(x) (pi^2) * sin(pi * x) + exp(x) * expscale;
x0 = 0;
x1 = 1;

fem = Pk1D(ord);
hlist = zeros(1, Nref);
errL2 = zeros(1, Nref);
errH1 = zeros(1, Nref);

for cycle = 1 : Nref
    hlist(cycle) = h0;
    grid = x0 : h0 : x1;
    assert(length(grid) >= 3);

    A = assembleK_Poi1D(fem, grid);
    F = assembleLoad_1D(fem, grid, f);
    
    % strong Dirichlet BDC
    N = size(A, 1);
    sol = zeros(N, 1);
    sol(1) = u_exact(x0);
    sol(end) = u_exact(x1);
    F = F - A * sol;
    sol(2:end-1) = A(2:end-1,2:end-1) \ F(2:end-1);

    [errH1(cycle), errL2(cycle), ~] = ...
        getH1Err_1D(fem, grid, sol, u_exact, du_exact);

    h0 = h0 / 2;
end

h = figure;
set(h, "Position", [100, 300, 1500, 400]);

subplot(1, 3, 1);
opt.nPoint = 20;
opt.u_exact = u_exact;
plotSol_1D(fem, grid, sol, opt);

subplot(1, 3, 2);
showrateh_mdf(hlist, errL2, Nref-1, '-o', "$||u-u_h||_{L^2}$");
title("$||u-u_h||_{L^2}$", "Interpreter", "latex");

subplot(1, 3, 3);
showrateh_mdf(hlist, errH1, Nref-1, '-o', "$||u-u_h||_{H^1}$");
title("$||u-u_h||_{H^1}$", "Interpreter", "latex");
