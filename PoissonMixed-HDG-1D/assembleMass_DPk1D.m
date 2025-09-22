function K = assembleMass_DPk1D(fem, grid)

    [quadL, w] = quadpts1(2 * fem.ord);
    nq = numel(w);
    Kloc = zeros(fem.locDof);
    for i = 1 : nq
        lam = quadL(i,:);
        phi = fem.getSpan(lam);
        Kloc = Kloc + phi(:) * phi(:)' * w(i);
    end

    NT = length(grid) - 1;
    h = grid(2:end) - grid(1:end-1);
    K = kron(spdiags(h(:), 0, NT, NT), Kloc);

end