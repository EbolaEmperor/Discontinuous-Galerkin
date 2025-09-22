function F = assembleLoad_Bihar1D(fem, grid, f)
    NT = length(grid) - 1;
    nDof = fem.locDof;
    k = fem.ord;
    F = zeros(NT*k+1, 1);

    [quadL, w] = quadpts1(2 * k + 2);
    w = w(:)';
    nq = numel(w);
    phi = zeros(nDof, nq);
    for i = 1 : nq
        lam = quadL(i, :);
        phi(:,i) = fem.getSpan(lam)';
    end

    for t = 1 : NT
        vtx = grid(t:t+1);
        fq = zeros(1, nq);
        for i = 1 : nq
            lam = quadL(i, :);
            pts = dot(lam, vtx);
            fq(i) = f(pts);
        end
        idx = (t-1) * k + (1:k+1);
        F(idx) = F(idx) + sum(phi .* fq .* w, 2) * (vtx(2) - vtx(1));
    end
end