%% Lagrange Pk element for Poisson equation

clc
clear
close all
tic
ord = 4;
h0 = 0.5;
domain = Polygon(6);
Nref = 6;
fun = sinsin(0.3);

u_exact = fun.u_exact;
grad_u_exact = fun.grad_u_exact;
lap_u_exact = fun.laplace_u_exact;

hlist = zeros(Nref, 1);
errL2 = zeros(Nref, 1);
errH1 = zeros(Nref, 1);

flag = 0;

for lv = 1 : Nref
    hlist(lv) = h0;
    [node, elem] = domain.getMesh(h0);
    
    fem = Pk(ord, node, elem);
    [elem2dof, nDof] = fem.getDOF(elem);
    
    A = assembleK_Poi2D(fem, node, elem, elem2dof);
    F = -assembleLoadVector(fem, node, elem, elem2dof, lap_u_exact);

    [c, freeDof] = interpStrongBDC(fem, node, elem, elem2dof, domain, u_exact);
    F = F - A * c;
    c(freeDof) = A(freeDof, freeDof) \ F(freeDof);

    [errH1(lv), errL2(lv), ~] = ...
        getH1Err(fem, node, elem, elem2dof, c, u_exact, grad_u_exact);

    if (h0 <= 1/16 || lv == Nref) && flag == 0
        flag = 1;
        h = figure;
        set(h, "Position", [100, 300, 1500, 400]);

        subplot(1, 3, 1);
        plotSol(fem, node, elem, c, elem2dof);
        ss = sprintf("$u_h\\;\\left(h=\\frac{1}{%d}\\right)$", 1/h0);
        title(ss, "Interpreter", "latex");
    end

    h0 = h0 / 2;
end

subplot(1, 3, 2);
showrateh_mdf(hlist, errL2, Nref-1, '-o', "$||u-u_h||_{L^2}$");
title("$||u-u_h||_{L^2}$", "Interpreter", "latex");

subplot(1, 3, 3);
showrateh_mdf(hlist, errH1, Nref-1, '-o', "$||u-u_h||_{H^1}$");
title("$||u-u_h||_{H^1}$", "Interpreter", "latex");
toc