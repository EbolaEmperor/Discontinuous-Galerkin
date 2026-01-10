function [u, info] = solveBratu_FEM_DeflatedNewton(fem, grid, A, lambda, u0, uKnown, opt)
%SOLVEBRATU_FEM_DEFLATEDNEWTON Deflated Newton to avoid known Bratu solutions.
%   Solve the classical Bratu nonlinear system using a deflated residual:
%       R(u)   = A*u - lambda*F(u)   (Bratu: u'' + lambda*exp(u) = 0)
%       R~(u)  = M(u) * R(u),
%   where M(u) penalizes convergence to solutions in uKnown.

    nDof = size(A, 1);
    if nargin < 7 || isempty(opt)
        opt = struct();
    end
    opt = fillDefaults(opt, struct( ...
        "I", 2:(nDof - 1), ...
        "tol", 1e-12, ...
        "maxIt", 60, ...
        "verbose", true, ...
        "lsMax", 12, ...
        "lsC1", 1e-4, ...
        "deflate_p", 2, ...
        "deflate_alpha", 1.0, ...
        "deflate_eps", 1e-14 ...
    ));

    u = u0;

    if isempty(uKnown)
        error("uKnown must contain at least one known solution for deflation.");
    end
    if isvector(uKnown)
        uKnown = uKnown(:);
    end
    if size(uKnown, 1) ~= nDof
        error("uKnown must have size nDof-by-k.");
    end

    info.resHist = zeros(opt.maxIt + 1, 1);
    info.converged = false;

    for it = 0 : opt.maxIt
        [Fnl, Jnl] = assemble_Bratu_terms(fem, grid, u);
        R = A * u - lambda * Fnl;
        J = A - lambda * Jnl;

        [M, gI] = deflationFactor(u, uKnown, opt.I, opt.deflate_p, opt.deflate_alpha, opt.deflate_eps);
        Rdef = M * R(opt.I);
        res = norm(Rdef);
        info.resHist(it + 1) = res;
        if opt.verbose
            fprintf("it=%2d  ||R_def||_2=%.3e\n", it, res);
        end
        if res < opt.tol
            info.converged = true;
            break;
        end

        JII = J(opt.I, opt.I);
        Jdef = M * full(JII) + (R(opt.I) * gI.');
        duI = -(Jdef \ Rdef);

        du = zeros(nDof, 1);
        du(opt.I) = duI;

        step = 1.0;
        accepted = false;
        for ls = 1 : opt.lsMax
            u_try = u + step * du;
            [F_try, ~] = assemble_Bratu_terms(fem, grid, u_try);
            R_try = A * u_try - lambda * F_try;
            [M_try, ~] = deflationFactor(u_try, uKnown, opt.I, opt.deflate_p, opt.deflate_alpha, opt.deflate_eps);
            Rdef_try = M_try * R_try(opt.I);
            if norm(Rdef_try) <= (1 - opt.lsC1 * step) * res
                u = u_try;
                u(1) = 0;
                u(end) = 0;
                accepted = true;
                break;
            end
            step = step / 2;
        end

        if ~accepted
            warning("Deflated Newton line-search failed to decrease deflated residual; stopping.");
            break;
        end
    end

    info.it = it;
    info.res = info.resHist(it + 1);
    info.I = opt.I;
    info.deflate_p = opt.deflate_p;
    info.deflate_alpha = opt.deflate_alpha;
end

function [M, gI] = deflationFactor(u, uKnown, I, p, alpha, epsDist)
    % M(u) = prod_i ( 1/||u-u_i||^p + alpha )
    % gI = dM/du restricted to interior dofs I (column vector)

    k = size(uKnown, 2);
    mi = zeros(k, 1);
    dists = zeros(k, 1);
    for j = 1:k
        d = u(I) - uKnown(I, j);
        dist = norm(d);
        dist = max(dist, epsDist);
        dists(j) = dist;
        mi(j) = dist.^(-p) + alpha;
    end

    M = prod(mi);

    gI = zeros(numel(I), 1);
    for j = 1:k
        d = u(I) - uKnown(I, j);
        dist = dists(j);
        dmi = -p * d * dist.^(-p-2);
        gI = gI + (M / mi(j)) * dmi;
    end
end
