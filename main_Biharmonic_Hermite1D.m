clc
clear
close all

ord = 2;
Nref = 4;
h0 = 0.25;
expscale = 1 / 3;
u_exact = @(x) sin(pi * x).^2 - exp(x) * expscale;
du_exact = @(x) pi * sin(2 * pi * x) - exp(x) * expscale;
d2u_exact = @(x) 2*pi^2 * cos(2 * pi * x) - exp(x) * expscale;
d3u_exact = @(x) -4*pi^3 * sin(2 * pi * x) - exp(x) * expscale;
d4u_exact = @(x) -8*pi^4 * cos(2 * pi * x) - exp(x) * expscale;
x0 = 0;
x1 = 1;

assert(ord >= 2);
hlist = zeros(1, Nref);
errL2 = zeros(1, Nref);
errH1 = zeros(1, Nref);


for cycle = 1 : Nref
    hlist(cycle) = h0;
    grid = x0 : h0 : x1;
    assert(length(grid) >= 3);

    fem = Hermite1D(ord, h0);
    A = assembleK_Bihar1D(fem, grid);
    F = assembleLoad_1D(fem, grid, d4u_exact);

    % strong BDC for u(x0), u(x1) and u'(x0), u'(x1)
    N = size(A, 1);
    sol = zeros(N, 1);
    sol(1) = u_exact(x0);
    sol(2) = du_exact(x0);
    sol(N-1) = u_exact(x1);
    sol(N) = du_exact(x1);
    F = F - A * sol;

    frD = 3:N-2;
    sol(frD) = A(frD,frD) \ F(frD);

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
