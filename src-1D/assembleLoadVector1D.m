function F = assembleLoadVector1D(fem, grid, f)
    NT = length(grid) - 1;
    nDof = fem.locDof;
    F = zeros(nDof, NT);

    [quadL, w] = quadpts1(2 * (nDof - 1));
    w = w(:)';
    nq = numel(w);
    phi = zeros(nDof, nq);
    for i = 1 : nq
        lam = quadL(i, :);
        phi(:,i) = fem.getSpan(lam)';
    end

    parfor t = 1 : NT
        vtx = grid(t:t+1);
        fq = zeros(1, nq);
        for i = 1 : nq
            lam = quadL(i, :);
            pts = dot(lam, vtx);
            fq(i) = f(pts);
        end
        F(:,t) = sum(phi .* fq .* w, 2) * (vtx(2) - vtx(1));
    end

    F = F(:);
end