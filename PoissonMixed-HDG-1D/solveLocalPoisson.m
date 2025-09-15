function uh_star = solveLocalPoisson(fem_star, fem1, fem2, sigmah, uh, x0, x1, f)

[quadL, w] = quadpts1(fem_star.ord * 2);
nq = numel(w);

nDofstar = fem_star.locDof;

h = x1 - x0;
K = zeros(nDofstar, nDofstar);
F = zeros(nDofstar, 1);
int_mean_star = zeros(1, nDofstar);
int_mean = 0;

for i = 1 : nq
    lam = quadL(i,:);
    pts = dot(lam, [x0, x1]);

    phi2 = fem2.getSpan(lam);
    phi = fem_star.getSpan(lam);
    dphi = fem_star.diffSpan(lam, h);
    f_q = f(pts);

    K = K + dphi(:) * dphi(:)' * w(i) * h;
    F = F + f_q * phi(:) * w(i) * h;

    int_mean_star = int_mean_star + phi * w(i);
    int_mean = int_mean + phi2 * uh * w(i);
end

sigma0 = - fem1.getSpan([1,0]) * sigmah;
v0 = fem_star.getSpan([1,0]);
F = F + sigma0 * v0(:);

sigma1 = fem1.getSpan([0,1]) * sigmah;
v1 = fem_star.getSpan([0,1]);
F = F + sigma1 * v1(:);

uh_star = [K; int_mean_star] \ [F; int_mean];

end