function plotSol_DPk1D(fem, grid, sol, opt)

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
    dofIdx = (i-1)*nDof + (1:nDof);
    idx = (i-1)*nPoint + (1:nPoint);
    x(idx) = grid(i) + locGrid * (grid(i+1) - grid(i));
    y(idx) = sum(phi .* sol(dofIdx));
end

plot(x, y, "-", "LineWidth", 1.5);

if isfield(opt, "varname")
    varname = sprintf("$%s_h$, $DP(%d)$", opt.varname, fem.ord);
else
    varname = "IPDG";
end

if isfield(opt, "u_exact")
    u = opt.u_exact(x);
    hold on;
    plot(x, u, "--", "LineWidth", 1.5);
    legend(varname, "exact", "Location", "best", "Interpreter", "latex");
else
    legend(varname, "Location", "best", "Interpreter", "latex");
end

end