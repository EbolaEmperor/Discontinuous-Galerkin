%% Bratu equation: FD iterative solver (standalone)
%   u'' + lambda * exp(u) = 0
%   u(0) = u(1) = 0
clear; clc;
lambda = 1;
n = 100;

set(gcf, "Position", [100, 100, 1200, 420]);
opt.axSol = subplot(1, 2, 1);
opt.axConv = subplot(1, 2, 2);
opt.nIter = 7;
opt.method = "direct";
opt.alpha = 10;
opt.recJacobian = true;

solveBratu_FDIter(lambda, n, opt);
