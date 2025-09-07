clc
clear
close all

ord = 2;
Nref = 5;
h0 = 0.25;
sigma = 5 * ord * (ord + 1);
u_exact = @(x) sin(pi * x);
f = @(x) (pi^2) * sin(pi * x);
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
    B = assembleBdry1D(fem, grid, sigma);
    A = K + P + B;

    F = assembleLoadVector1D(fem, grid, f);
    sol = A \ F;

    h0 = h0 / 2;
end

opt.nPoint = 20;
opt.u_exact = u_exact;
plotSol1D(fem, grid, sol, opt);