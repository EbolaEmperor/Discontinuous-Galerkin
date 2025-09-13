clc
clear
close all

delpath
addpath("Poisson-IPDG-1D");
addpath("tools");

ord = 2;
Nref = 7;
h0 = 0.5;
sigma = 5 * ord * (ord + 1);
u_exact = @(x) sin(pi * x) - exp(x)/2;
du_exact = @(x) pi * cos(pi * x) - exp(x)/2;
f = @(x) (pi^2) * sin(pi * x) + exp(x)/2;
x0 = 0;
x1 = 1;
IP_type = "SIPG";
% You can also try NIPG or IIPG

fem = DPk(ord);
hlist = zeros(1, Nref);
errL2 = zeros(1, Nref);
errH1 = zeros(1, Nref);

beta = 1;
if IP_type == "NIPG", beta = -1; end
if IP_type == "IIPG", beta = 0; end

for cycle = 1 : Nref
    hlist(cycle) = h0;
    grid = x0 : h0 : x1;
    assert(length(grid) >= 3);

    K = assembleStiffness(fem, grid);
    P = assembleInnerPenalty(fem, grid, sigma, beta);
    A = K + P;

    F = assembleLoadVector(fem, grid, f);

    % Nitsche weak Dirichlet BDC
    [Ahat, Fhat] = assembleBdryNitsche(fem, grid, sigma, beta, u_exact(x0), u_exact(x1));
    A = A + Ahat;
    F = F + Fhat;
    sol = A \ F;
    
    % strong Dirichlet BDC
    % N = size(A, 1);
    % sol = zeros(N, 1);
    % sol(1) = u_exact(x0);
    % sol(end) = u_exact(x1);
    % F = F - A * sol;
    % sol(2:end-1) = A(2:end-1,2:end-1) \ F(2:end-1);

    [errH1(cycle), errL2(cycle), ~] = ...
        getH1Err(fem, grid, sol, u_exact, du_exact);

    h0 = h0 / 2;
end

h = figure;
set(h, "Position", [100, 300, 1500, 400]);

subplot(1, 3, 1);
opt.nPoint = 20;
opt.u_exact = u_exact;
plotSol(fem, grid, sol, opt);

subplot(1, 3, 2);
showrateh_mdf(hlist, errL2, Nref-1, '-o', "$||u-u_h||_{L^2}$");
title("$||u-u_h||_{L^2}$", "Interpreter", "latex");

subplot(1, 3, 3);
showrateh_mdf(hlist, errH1, Nref-1, '-o', "$||u-u_h||_{H^1}$");
title("$||u-u_h||_{H^1}$", "Interpreter", "latex");
