function [B, F] = assembleBdryNitsche_Bihar1D(fem, grid, sigma, beta, g0, g1)

NT = length(grid) - 1;
k = fem.ord;
B = sparse(NT*k+1, NT*k+1);
F = sparse(NT*k+1, 1);

%% 左端点的边界条件
h0 = grid(2) - grid(1);
dphi = fem.diffSpan([1,0], h0);
d2phi = -fem.diff2Span([1,0], h0);
idx = 1:k+1;
B(idx,idx) = - beta * d2phi' * dphi - dphi' * d2phi + (sigma/h0) * (dphi' * dphi);
F(idx) = - beta * d2phi * g0 + (sigma/h0) * dphi * g0;

%% 右端点的边界条件
h1 = grid(end) - grid(end-1);
dphi = fem.diffSpan([0,1], h1);
d2phi = fem.diff2Span([0,1], h1);
idx = (NT-1)*k + (1:k+1);
B(idx,idx) = - beta * d2phi' * dphi - dphi' * d2phi + (sigma/h1) * (dphi' * dphi);
F(idx) = - beta * d2phi * g1 + (sigma/h1) * dphi * g1;

end