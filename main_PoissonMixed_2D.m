clc
clear
close all

ord = 4;
alpha = 1;

Nref = 5;
h0 = 0.5;
domain = square();
fun = sinsin(0.3);

u_exact = fun.u_exact;
grad_u_exact = fun.grad_u_exact;
lap_u_exact = fun.laplace_u_exact;
f = @(x) -lap_u_exact(x);

hlist = zeros(1, Nref);
err_u = zeros(1, Nref);
err_ustar = zeros(1, Nref);
err_sigma = zeros(1, Nref);

for lv = 1 : Nref
    hlist(lv) = h0;
    [node, elem] = domain.getMesh(h0);
    [edge, edge2side] = getEdge2Side(node, elem);

    %% Solve mixed-form Poisson equation using HDG method.
    fem1 = vecDPk(ord, node, elem);
    [elem2dof1, nDof1] = fem1.getDOF(elem);

    fem2 = DPk(ord, node, elem);
    [elem2dof2, nDof2] = fem2.getDOF(elem);
    
    M  = assembleMass(fem1, node, elem, elem2dof1);
    B1 = assembleDivMass(fem1, fem2, node, elem, elem2dof1, elem2dof2);
    B2 = assembleGradMass(fem2, fem1, node, elem, elem2dof2, elem2dof1);
    K  = [M, B1; B2, sparse(nDof2,nDof2)];

    P = assembleIP_HDG2D(fem1, fem2, node, elem, elem2dof1, elem2dof2, edge, edge2side, alpha);
    A = K + P;

    F = [zeros(nDof1,1);
         assembleLoadVector(fem2, node, elem, elem2dof2, f)];
    % Apply non-homo Dirichlet BDC
    F = F + assembleWeakBDC_HDG2D(fem1, fem2, node, elem, elem2dof1, elem2dof2, edge, edge2side, alpha, u_exact);

    sol = A \ F;
    sigmah = sol(1:nDof1);
    uh     = sol(nDof1+1:end);

    err_sigma(lv) = getL2Err(fem1, node, elem, elem2dof1, sigmah, grad_u_exact);
    err_u(lv)     = getL2Err(fem2, node, elem, elem2dof2, uh, u_exact);

    %% Post-process: Solve local Poisson equations to obtain the super-convergence.
    fem_star = DPk(ord+1, node, elem);
    [elem2dof_star, ~] = fem_star.getDOF(elem);

    uh_star = solveLocalPoisson2D_HDG2D( ...
        fem_star, fem1, fem2, node, elem, ...
        elem2dof_star, elem2dof1, elem2dof2, ...
        sigmah, uh, f);

    err_ustar(lv) = getL2Err(fem_star, node, elem, elem2dof_star, uh_star, u_exact);

    h0 = h0 / 2;
end

h = figure;
set(h, "Position", [100, 100, 1100, 900]);

subplot(2, 2, 1);
plotSol(fem_star, node, elem, uh_star, elem2dof_star);
title("$u_h$", "Interpreter", "latex");

subplot(2, 2, 2);
plotSolVec(fem1, node, elem, sigmah, elem2dof1);
title("$\sigma_h$", "Interpreter", "latex");

subplot(2, 2, 3);
showrateh_mdf(hlist, err_ustar, Nref-1, '-o', "$\|u-u_h^\star\|_{L^2}$");
title("$\|u-u_h^\star\|_{L^2}$", "Interpreter", "latex");

subplot(2, 2, 4);
showrateh_mdf(hlist, err_sigma, Nref-1, '-o', "$\|\nabla u-\sigma_h\|_{L^2}$");
title("$\|\nabla u-\sigma_h\|_{L^2}$", "Interpreter", "latex");
