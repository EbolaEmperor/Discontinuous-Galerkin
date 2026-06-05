function F = assembleLoad_1D(fem, grid, f)
    NT = length(grid) - 1;
    nDof = fem.locDof;
    F = zeros(fem.nDof(NT), 1);

    [quadL, w] = quadpts1_my(2 * fem.locDof);
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
        idx = fem.dofMap(t);
        F(idx) = F(idx) + sum(phi .* fq .* w, 2) * (vtx(2) - vtx(1));
    end
end