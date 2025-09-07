function B = assembleBdry1D(fem, grid, sigma)

NT = length(grid) - 1;
nDof = fem.locDof;
B = sparse(NT*nDof, NT*nDof);

%% 左端点的边界条件
h0 = grid(2) - grid(1);
phi = fem.getSpan([1,0]);
dphi = -fem.diffSpan([1,0], h0);
idx = 1:nDof;
B(idx,idx) = - dphi' * phi - phi' * dphi + (sigma/h0^2) * (phi' * phi);

%% 右端点的边界条件
h1 = grid(end) - grid(end-1);
phi = fem.getSpan([0,1]);
dphi = fem.diffSpan([0,1], h1);
idx = (NT-1)*nDof + (1:nDof);
B(idx,idx) = - dphi' * phi - phi' * dphi + (sigma/h1^2) * (phi' * phi);

end