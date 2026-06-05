function out = solveBratu_FDIter(lambda, n, opt)
%SOLVEBRATU_FDITER Fixed-point FD iteration for 1D Bratu equation.
%   Solve: u'' + lambda*exp(u) = 0 on (0,1), with u(0)=u(1)=0.

    opt = fillDefaults(opt, struct( ...
        "nIter", 10, ...
        "axSol", [], ...
        "axConv", [], ...
        "method", "direct", ...
        "alpha", 1, ...
        "branch", 1, ...
        "recJacobian", false ...
    ));

    h = 1 / (n + 1);

    % Semi-analytic exact solution (for validation):
    %   u(x)  = 2*log( cosh(beta/2) / cosh(beta*(x-1/2)) )
    % where beta>0 solves lambda = 2*beta^2 / cosh(beta/2)^2.
    betas = solveBratuBeta(lambda);
    betas = sort(betas(:));
    if isempty(betas)
        error("No real beta found for this lambda; no Bratu exact solution branch available.");
    end
    beta_exact = betas(opt.branch);
    [u_exact, du_exact] = bratuExact1D(beta_exact);

    x_in = h * (1:n).';
    x_full = h * (0:n+1);

    u = opt.alpha * x_in .* (1-x_in);
    % u = u_exact(x_in);
    errIter = zeros(opt.nIter + 1, 1);
    errIter(1) = fdNodeL2Error(u, x_in, h, u_exact);
    A = (diag(2*ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1)) / (h^2);
    if opt.recJacobian
        jacobian = zeros(n, n, opt.nIter);
    end
    for k = 1 : opt.nIter
        if opt.method == "direct"
            if opt.recJacobian
                jacobian(:,:,k) = lambda * inv(A) .* exp(u)';
                lam = svd(jacobian(:,:,k));
                fprintf("Jacobian: lambda in [%f, %f]\n", lam(end), lam(1));
            end
            u = A \ (lambda * exp(u));
        end
        if opt.method == "jacobi"
            r1 = lambda * exp(u);
            r2 = [0; u(1:end-1)] / (h^2);
            r3 = [u(2:end); 0] / (h^2);
            u = h^2/2 * (r1 + r2 + r3);
        end
        errIter(k+1) = fdNodeL2Error(u, x_in, h, u_exact);
    end

    [pIter, qIter] = fitConvergenceOrder(errIter);
    fprintf("FD iteration fit: q≈%.6g, p≈%.6g (based on node L2-error)\n", qIter, pIter);

    if ~isempty(opt.axSol) && ishandle(opt.axSol)
        plot(opt.axSol, x_full, [0; u; 0], "LineWidth", 1.5);
        hold(opt.axSol, "on");
        plot(opt.axSol, x_full, u_exact(x_full), "--", "LineWidth", 1.5);
        hold(opt.axSol, "off");
        legend(opt.axSol, "FD-iter", "exact", "Location", "best");
        title(opt.axSol, "Iterative Solver");
    end

    if ~isempty(opt.axConv) && ishandle(opt.axConv)
        it = 0:opt.nIter;
        semilogy(opt.axConv, it, errIter, "--", "LineWidth", 1.5);
        xlabel(opt.axConv, "iteration");
        legend(opt.axConv, "FD error to exact (node L2)", "Location", "best");
        title(opt.axConv, sprintf("Convergence: q≈%.3g, p≈%.3g", qIter, pIter));
    end

    out.lambda = lambda;
    out.n = n;
    out.nIter = opt.nIter;
    out.h = h;
    out.x_in = x_in;
    out.x_full = x_full;
    out.u = u;
    out.errIter = errIter;
    out.pIter = pIter;
    out.qIter = qIter;
    out.beta_exact = beta_exact;
    out.u_exact = u_exact;
    out.du_exact = du_exact;
    if opt.recJacobian
        out.jacobian = jacobian;
    end
end
