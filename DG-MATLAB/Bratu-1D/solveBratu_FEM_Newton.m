function [u, info] = solveBratu_FEM_Newton(fem, grid, A, lambda, u0, opt)
%SOLVEBRATU_FEM_NEWTON Newton solver for the 1D Bratu FEM nonlinear system.
%   Bratu: u'' + lambda * exp(u) = 0, u(0)=u(1)=0.
%   Residual (weak form): A*u - lambda*F(u) = 0, where F_i = \int exp(u_h)*phi_i dx.

    nDof = size(A, 1);
    if nargin < 6 || isempty(opt)
        opt = struct();
    end
    opt = fillDefaults(opt, struct( ...
        "I", 2:(nDof - 1), ...
        "tol", 1e-12, ...
        "maxIt", 50, ...
        "verbose", true, ...
        "lsMax", 12, ...
        "lsC1", 1e-4 ...
    ));

    u = u0;
    u(1) = 0;
    u(end) = 0;

    info.resHist = zeros(opt.maxIt + 1, 1);
    info.stepHist = zeros(opt.maxIt, 1);
    info.converged = false;

    for it = 0 : opt.maxIt
        [Fnl, Jnl] = assemble_Bratu_terms(fem, grid, u);
        R = A * u - lambda * Fnl;         % u'' + lambda*exp(u) = 0  (i.e. -u'' = lambda*exp(u))
        res = norm(R(opt.I));
        info.resHist(it + 1) = res;
        if opt.verbose
            fprintf("it=%2d  ||R||_2=%.3e\n", it, res);
        end
        if res < opt.tol
            info.converged = true;
            break;
        end
        if it == opt.maxIt
            break;
        end

        J = A - lambda * Jnl;
        du = zeros(nDof, 1);
        du(opt.I) = -(J(opt.I, opt.I) \ R(opt.I));

        step = 1.0;
        accepted = false;
        for ls = 1 : opt.lsMax
            u_try = u + step * du;
            [F_try, ~] = assemble_Bratu_terms(fem, grid, u_try);
            R_try = A * u_try - lambda * F_try;
            if norm(R_try(opt.I)) <= (1 - opt.lsC1 * step) * res
                u = u_try;
                u(1) = 0;
                u(end) = 0;
                accepted = true;
                break;
            end
            step = step / 2;
        end

        info.stepHist(it + 1) = step;
        if ~accepted
            warning("Newton line-search failed to decrease residual; stopping.");
            break;
        end
    end

    info.it = it;
    info.res = info.resHist(it + 1);
    info.I = opt.I;
end
