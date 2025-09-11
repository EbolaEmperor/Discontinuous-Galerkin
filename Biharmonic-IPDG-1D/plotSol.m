function plotSol(fem, grid, sol, opt)

if nargin < 4, opt.nPoint = 10; end

nPoint = opt.nPoint;
k = fem.ord;
nDof = fem.locDof;
locGrid = 0 : 1/(nPoint-1) : 1;
phi = zeros(nDof, nPoint);
for i = 1 : nPoint
    lam = [1-locGrid(i), locGrid(i)];
    phi(:,i) = fem.getSpan(lam)';
end

NT = length(grid) - 1;
x = zeros(1, NT * nPoint);
y = zeros(1, NT * nPoint);

for i = 1 : NT
    dofIdx = (i-1)*k + (1:k+1);
    idx = (i-1)*nPoint + (1:nPoint);
    x(idx) = grid(i) + locGrid * (grid(i+1) - grid(i));
    y(idx) = sum(phi .* sol(dofIdx));
end

plot(x, y, "-", "LineWidth", 1.5);

if isfield(opt, "u_exact")
    u = opt.u_exact(x);
    hold on;
    plot(x, u, "--", "LineWidth", 1.5);
    legend("IPCG-C0Pk", "exact", "Location", "best");
else
    legend("IPCG-C0Pk", "Location", "best");
end

end