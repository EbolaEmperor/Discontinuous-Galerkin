clc
clear
close all

delpath
addpath("PoissonMixed-HDG-2D");
addpath("common-2D");
addpath("tools");

ord = 1;
alpha = ord * (ord+1);

Nref = 6;
h0 = 0.5;
domain = square();
fun = sinsin(0);

u_exact = fun.u_exact;
grad_u_exact = fun.grad_u_exact;
lap_u_exact = fun.laplace_u_exact;
f = @(x) -lap_u_exact(x);

hlist = zeros(1, Nref);
err_u = zeros(1, Nref);
err_ustar = zeros(1, Nref);
err_sigma = zeros(1, Nref);

for cycle = 1 : Nref
    hlist(lv) = h0;
    [node, elem] = domain.getMesh(h0);
    [edge, edge2side] = getEdge2Side(node, elem);

    fem1 = vecDPk(ord, node, elem);
    [elem2dof1, nDof1] = fem1.getDOF(elem);

    fem2 = DPk(ord, node, elem);
    [elem2dof2, nDof2] = fem2.getDOF(elem);
    
    M = assembleMass(fem1, node, elem, elem2dof1);
    B1 = assembleDivMass(fem1, fem2, node, elem, elem2dof1, elem2dof2);
    B2 = assembleGradMass(fem2, fem1, node, elem, elem2dof2, elem2dof1);
    K = [M, B1; B2, sparse(nDof2,nDof2)];

    h0 = h0 / 2;
end

% h = figure;
% set(h, "Position", [100, 100, 900, 1000]);
