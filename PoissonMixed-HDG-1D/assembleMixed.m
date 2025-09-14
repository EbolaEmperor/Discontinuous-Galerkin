function K = assembleMixed(fem1, fem2, grid)

    [quadL, w] = quadpts1(fem1.ord + fem2.ord);
    nq = numel(w);
    Kloc = zeros(fem1.locDof, fem2.locDof);
    for i = 1 : nq
        lam = quadL(i,:);
        dphi = fem1.diffSpan(lam, 1);
        phi = fem2.getSpan(lam, 1);
        Kloc = Kloc + dphi(:) * phi(:)' * w(i);
    end

    NT = length(grid) - 1;
    hinv = 1 ./ (grid(2:end) - grid(1:end-1));
    K = kron(spdiags(hinv(:), 0, NT, NT), Kloc);

end