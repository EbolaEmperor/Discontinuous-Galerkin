function betas = solveBratuBeta(lambda)
%SOLVEBRATUBETA Find beta values for the 1D Bratu exact solution branches.
%   For the classical Bratu problem:
%       u'' + lambda * exp(u) = 0,  u(0)=u(1)=0,
%   the exact solution is parameterized by beta>0 with
%       lambda = 2*beta^2 / cosh(beta/2)^2.
%   For 0 < lambda < lambda*, there are two betas (two branches).

    assert(lambda > 0);
    g = @(b) 2 * b.^2 ./ (cosh(b / 2).^2) - lambda;

    bMax = 50;
    nScan = 4000;
    bs = linspace(1e-10, bMax, nScan);
    vals = g(bs);
    s = sign(vals);
    s(s == 0) = 1;
    idx = find(s(1:end-1) .* s(2:end) < 0);

    betas = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        a = bs(idx(k));
        b = bs(idx(k) + 1);
        betas(k) = fzero(g, [a, b]);
    end
end
