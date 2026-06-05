function plotSol_1D(fem, grid, sol, opt)

if nargin < 4, opt.nPoint = 10; end

nPoint = opt.nPoint;
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
    dofIdx = fem.dofMap(i);
    idx = (i-1)*nPoint + (1:nPoint);
    h = grid(i+1) - grid(i);
    x(idx) = grid(i) + locGrid * h;
    y(idx) = sum(phi .* sol(dofIdx));
end

plot(x, y, "-", "LineWidth", 1.5);
str = sprintf("IPDG-C0P%d", fem.ord);

if isfield(opt, "u_exact")
    u = opt.u_exact(x);
    hold on;
    plot(x, u, "--", "LineWidth", 1.5);
    legend(str, "exact", "Location", "best");
else
    legend(str, "Location", "best");
end

end