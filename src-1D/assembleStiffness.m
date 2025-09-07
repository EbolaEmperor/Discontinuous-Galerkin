function K = assembleStiffness(fem, grid)

    [quadL, w] = quadpts1(2 * fem.ord);
    nq = numel(w);
    Kloc = zeros(fem.locDof);
    for i = 1 : nq
        lam = quadL(i,:);
        dphi = fem.diffSpan(lam, 1);
        Kloc = Kloc + dphi(:) * dphi(:)';
    end

    NT = length(grid) - 1;
    hinv = 1 ./ (grid(2:end) - grid(1:end-1));
    K = kron(spdiags(hinv(:), 0, NT, NT), Kloc);

end