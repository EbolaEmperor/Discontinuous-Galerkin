function K = assembleMixed_HDG1D(fem1, fem2, grid)

    [quadL, w] = quadpts1(fem1.ord + fem2.ord);
    nq = numel(w);
    Kloc = zeros(fem1.locDof, fem2.locDof);
    for i = 1 : nq
        lam = quadL(i,:);
        dphi = fem1.diffSpan(lam,1);
        psi = fem2.getSpan(lam);
        Kloc = Kloc + dphi(:) * psi(:)' * w(i);
    end

    NT = length(grid) - 1;
    K = kron(spdiags(ones(NT,1), 0, NT, NT), Kloc);

end