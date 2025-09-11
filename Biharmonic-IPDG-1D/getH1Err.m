%% 以下为无向量化代码
function [errH1, errL2, errSemiH1] = getH1Err(fem, grid, uh, u, du)
    NT = length(grid) - 1;
    k = fem.ord;
    nDof = fem.locDof;
    errL2 = 0;
    errSemiH1 = 0;

    [quadL, w] = quadpts1(2 * k + 2);
    w = w(:)';
    nq = numel(w);
    phi = zeros(nDof, nq);
    dphi = zeros(nDof, nq);
    for i = 1 : nq
        lam = quadL(i, :);
        phi(:,i) = fem.getSpan(lam)';
        dphi(:,i) = fem.diffSpan(lam, 1)';
    end

    for t = 1 : NT
        vtx = grid(t:t+1);
        h = vtx(2) - vtx(1);
        uh_t = uh((t-1)*k + (1:k+1))';
        u_ex = zeros(1, nq);
        du_ex = zeros(1, nq);
        for i = 1 : nq
            lam = quadL(i, :);
            pts = dot(lam, vtx);
            u_ex(i) = u(pts);
            du_ex(i) = du(pts);
        end
        errL2 = errL2 + sum( (uh_t * phi - u_ex).^2 .* w ) * h;
        errSemiH1 = errSemiH1 + sum( (uh_t * dphi/h - du_ex).^2 .* w ) * h;
    end

    errH1 = sqrt(errL2 + errSemiH1);
    errL2 = sqrt(errL2);
    errSemiH1 = sqrt(errSemiH1);
end