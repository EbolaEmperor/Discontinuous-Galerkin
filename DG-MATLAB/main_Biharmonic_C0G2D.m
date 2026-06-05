clc
clear
close all
tic

ord = 4;
h0 = 0.5;
domain = square();
Nref = 5;
sigma = 1.5 * ord * (ord+1);
fun = sinsin(0.3);
IP_type = "SIPG";
% You can also try NIPG or IIPG

assert(ord >= 2);
u_exact = fun.u_exact;
grad_u_exact = fun.grad_u_exact;
lap_u_exact = fun.laplace_u_exact;
lap_lap_u_exact = fun.lap_lap_u_exact;

beta = 1;
if IP_type == "NIPG", beta = -1; end
if IP_type == "IIPG", beta = 0; end

hlist = zeros(Nref, 1);
errL2 = zeros(Nref, 1);
errH1 = zeros(Nref, 1);
flag = 0;

for lv = 1 : Nref
    hlist(lv) = h0;
    [node, elem] = domain.getMesh(h0);
    
    opt.useR = true;
    fem = Pk(ord, node, elem, opt);
    [elem2dof, nDof] = fem.getDOF(elem);
    [edge, edge2side] = getEdge2Side(node, elem);
    
    K = assembleK_Bihar2D(fem, node, elem, elem2dof);
    P = assembleIP_Bihar2D(fem, node, elem, elem2dof, edge, edge2side, sigma, beta);
    A = K + P;

    F = assembleLoadVector(fem, node, elem, elem2dof, lap_lap_u_exact);
    F_upd = modifyLoadVectorNitsche_Bihar2D(fem, node, elem, elem2dof, edge, edge2side, sigma, beta, grad_u_exact);
    F = F + F_upd;

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