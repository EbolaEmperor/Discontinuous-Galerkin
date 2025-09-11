clc
clear
close all

% delpath
% addpath("Poisson-IPDG-2D");
% addpath("tools");

ord = 1;
h0 = 1;
domain = square();

[node, elem] = domain.getMesh(h0);

fem = DPk(ord, node, elem);
[elem2dof, nDof] = fem.getDOF(elem);

K = assembleStiffness(fem, node, elem, elem2dof);
disp(full(K))
