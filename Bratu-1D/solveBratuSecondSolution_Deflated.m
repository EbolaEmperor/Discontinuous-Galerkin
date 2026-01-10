function [u2, info2] = solveBratuSecondSolution_Deflated(fem, grid, A, lambda, u1, opt)
%SOLVEBRATUSECONDSOLUTION_DEFLATED Find a second Bratu solution via deflation.
%   Uses deflated Newton with a small set of deterministic initial guesses.

    nDof = size(A, 1);
    if nargin < 6 || isempty(opt)
        opt = struct();
    end
    opt = fillDefaults(opt, struct( ...
        "I", 2:(nDof - 1), ...
        "minDist", 1e-10, ...
        "initScales", [0, 2, 5, 10, 20, 50, 100] / lambda ...
    ));

    info2 = struct();
    u2 = [];

    for s = opt.initScales
        u0 = s * u1;
        if s == 0
            u0 = zeros(nDof, 1);
        end

        [u_try, info_try] = solveBratu_FEM_DeflatedNewton(fem, grid, A, lambda, u0, u1, opt);
        if info_try.converged && norm(u_try(opt.I) - u1(opt.I)) > opt.minDist
            u2 = u_try;
            info2 = info_try;
            return;
        end
    end

    warning("Failed to find a distinct second solution with deflated Newton.");
end
