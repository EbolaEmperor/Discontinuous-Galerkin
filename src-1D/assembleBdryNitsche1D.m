function [B, F] = assembleBdryNitsche1D(fem, grid, sigma, g0, g1)

NT = length(grid) - 1;
nDof = fem.locDof;
B = sparse(NT*nDof, NT*nDof);
F = sparse(NT*nDof, 1);

%% 左端点的边界条件
h0 = grid(2) - grid(1);
phi = fem.getSpan([1,0]);
dphi = -fem.diffSpan([1,0], h0);
idx = 1:nDof;
B(idx,idx) = - dphi' * phi - phi' * dphi + (sigma/h0) * (phi' * phi);
F(idx) = - dphi * g0 + (sigma/h0) * phi * g0;

%% 右端点的边界条件
h1 = grid(end) - grid(end-1);
phi = fem.getSpan([0,1]);
dphi = fem.diffSpan([0,1], h1);
idx = (NT-1)*nDof + (1:nDof);
B(idx,idx) = - dphi' * phi - phi' * dphi + (sigma/h1) * (phi' * phi);
F(idx) = - dphi * g1 + (sigma/h1) * phi * g1;

end