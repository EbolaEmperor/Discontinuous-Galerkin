clc
clear
close all

delpath
addpath("Biharmonic-IPDG-1D");
addpath("tools");

ord = 3;
Nref = 5;
h0 = 0.25;
sigma = 3 * ord * (ord+1);
expscale = 1/3;
u_exact = @(x) sin(pi * x).^2 - exp(x) * expscale;
du_exact = @(x) pi * sin(2 * pi * x) - exp(x) * expscale;
d2u_exact = @(x) 2*pi^2 * cos(2 * pi * x) - exp(x) * expscale;
d3u_exact = @(x) -4*pi^3 * sin(2 * pi * x) - exp(x) * expscale;
d4u_exact = @(x) -8*pi^4 * cos(2 * pi * x) - exp(x) * expscale;
x0 = 0;
x1 = 1;
IP_type = "SIPG";
% You can also try NIPG or IIPG

assert(ord >= 2);
fem = Pk(ord);
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

    F = assembleLoadVector(fem, grid, d4u_exact);

    % Nitsche weak BDC for u'(x0) and u'(x1)
    [Ahat, Fhat] = assembleBdryNitsche(fem, grid, sigma, beta, du_exact(x0), du_exact(x1));
    A = A + Ahat;
    F = F + Fhat;
    
    % strong BDC for u(x0) and u(x1)
    N = size(A, 1);
    sol = zeros(N, 1);
    sol(1) = u_exact(x0);
    sol(end) = u_exact(x1);
    F = F - A * sol;

    sol(2:end-1) = A(2:end-1,2:end-1) \ F(2:end-1);

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
