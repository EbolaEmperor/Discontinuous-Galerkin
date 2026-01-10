function alpha = solveBratuAlpha(lambda)
%SOLVEBRATUALPHA Solve the semi-analytic parameter for 1D Bratu (compat).
%   Bratu: u'' + lambda*exp(u) = 0, u(0)=u(1)=0.
%   This returns the small-branch parameter (beta) satisfying:
%       lambda = 2*beta^2 / cosh(beta/2)^2, beta > 0.
%   Historical note: older code used a trig-based "alpha" for a different sign convention.

    assert(lambda > 0);
    betas = solveBratuBeta(lambda);
    betas = sort(betas(:));
    if isempty(betas)
        error("No real beta found for this lambda.");
    end
    alpha = betas(1);
end
