%% Bratu equation solver
%   u'' + lambda * exp(u) = 0
%   u(0) = u(1) = 0

clear; clc; close all;

% lambda = 3.5138;  % 存在两个解的临界 lambda 值
lambda = 0.1;
n = 50;
h = 1 / (n + 1);

%% Finite-Element Solver (two branches)
ord = 4;
fem = Pk1D(ord);
grid = 0 : h : 1;
NT = numel(grid) - 1;
nDof = fem.nDof(NT);
A = assembleK_Poi1D(fem, grid);

betas = solveBratuBeta(lambda);
betas = sort(betas(:));
if numel(betas) < 2
    error("Expected two Bratu branches for this lambda, but found %d beta root(s).", numel(betas));
end
[u_exact_1, du_exact_1] = bratuExact1D(betas(1));
[u_exact_2, du_exact_2] = bratuExact1D(betas(end));

optN.tol = 1e-12;
optN.maxIt = 50;
optN.verbose = true;

u0 = zeros(nDof, 1);
[u1, ~] = solveBratu_FEM_Newton(fem, grid, A, lambda, u0, optN);

optD = optN;
optD.deflate_p = 2;
optD.deflate_alpha = 1.0;
optD.maxIt = 80;
[u2, ~] = solveBratuSecondSolution_Deflated(fem, grid, A, lambda, u1, optD);
if isempty(u2)
    error("Deflated solve failed to produce a second solution.");
end

% Match numerical solutions to the two exact branches
[~, errL2_11, ~] = getH1Err_1D(fem, grid, u1, u_exact_1, du_exact_1);
[~, errL2_12, ~] = getH1Err_1D(fem, grid, u1, u_exact_2, du_exact_2);
[~, errL2_21, ~] = getH1Err_1D(fem, grid, u2, u_exact_1, du_exact_1);
[~, errL2_22, ~] = getH1Err_1D(fem, grid, u2, u_exact_2, du_exact_2);

if errL2_11 + errL2_22 <= errL2_12 + errL2_21
    u_small = u1;  u_large = u2;
    u_exact_small = u_exact_1; du_exact_small = du_exact_1;
    u_exact_large = u_exact_2; du_exact_large = du_exact_2;
else
    u_small = u2;  u_large = u1;
    u_exact_small = u_exact_1; du_exact_small = du_exact_1;
    u_exact_large = u_exact_2; du_exact_large = du_exact_2;
end

fprintf("Branch solve summary:\n");
[errH1s, errL2s, ~] = getH1Err_1D(fem, grid, u_small, u_exact_small, du_exact_small);
[errH1l, errL2l, ~] = getH1Err_1D(fem, grid, u_large, u_exact_large, du_exact_large);
fprintf("  small branch: errL2=%.3e  errH1=%.3e\n", errL2s, errH1s);
fprintf("  large branch: errL2=%.3e  errH1=%.3e\n", errL2l, errH1l);

set(gcf, "Position", [100, 100, 1200, 420]);
subplot(1, 2, 1);
opt.nPoint = 20;
opt.u_exact = u_exact_small;
plotSol_1D(fem, grid, u_small, opt);
title(sprintf("Branch 1 (small): errH1=%.3e", errH1s));

subplot(1, 2, 2);
opt.u_exact = u_exact_large;
plotSol_1D(fem, grid, u_large, opt);
title(sprintf("Branch 2 (large): errH1=%.3e", errH1l));
