function [u_exact, du_exact] = bratuExact1D(beta)
%BRATUEXACT1D Exact solution (and derivative) for 1D Bratu (classical).
%   Bratu: u'' + lambda*exp(u) = 0, u(0)=u(1)=0.
%   u(x)  = 2*log( cosh(beta/2) / cosh(beta*(x-1/2)) )
%   u'(x) = -2*beta*tanh(beta*(x-1/2))

    u_exact = @(x) 2 * log(cosh(beta / 2) ./ cosh(beta * (x - 0.5)));
    du_exact = @(x) -2 * beta * tanh(beta * (x - 0.5));
end
