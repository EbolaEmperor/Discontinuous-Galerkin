% clc
% clear
% close all

ord = 2;
Nref = 6;
h0 = 0.25;
sigma = 5 * ord * (ord + 1);
u_exact = @(x) sin(pi * x) + exp(x);
du_exact = @(x) pi * cos(pi * x) + exp(x);
f = @(x) (pi^2) * sin(pi * x) - exp(x);
x0 = 0;
x1 = 1;

fem = DPk1D(ord);
hlist = zeros(1, Nref);
errL2 = zeros(1, Nref);
errH1 = zeros(1, Nref);

for cycle = 1 : Nref
    hlist(cycle) = h0;
    grid = x0 : h0 : x1;
    assert(length(grid) >= 3);

    K = assembleStiffness1D(fem, grid);
    P = assembleInnerPenalty1D(fem, grid, sigma);
    % Nitsche weak BDC
    % B = assembleBdryNitsche1D(fem, grid, sigma);
    A = K + P;

    N = size(A, 1);
    % strong BDC
    A(1,1) = 1;
    A(1,2:end) = 0;
    A(end,end) = 1;
    A(end, 1:end-1) = 0;

    F = assembleLoadVector1D(fem, grid, f);
    F(1) = u_exact(x0);
    F(end) = u_exact(x1);
    sol = A \ F;

    [errH1(cycle), errL2(cycle), ~] = ...
        getH1Err1D(fem, grid, sol, u_exact, du_exact);

    h0 = h0 / 2;
end

h = figure;
set(h, "Position", [100, 300, 1500, 400]);

subplot(1, 3, 1);
opt.nPoint = 20;
opt.u_exact = u_exact;
plotSol1D(fem, grid, sol, opt);

subplot(1, 3, 2);
showrateh_mdf(hlist, errL2, Nref-1, '-o', "$||u-u_h||_{L^2}$");
title("$||u-u_h||_{L^2}$", "Interpreter", "latex");

subplot(1, 3, 3);
showrateh_mdf(hlist, errH1, Nref-1, '-o', "$||u-u_h||_{H^1}$");
title("$||u-u_h||_{H^1}$", "Interpreter", "latex");
