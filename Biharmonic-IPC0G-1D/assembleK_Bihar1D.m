function K = assembleK_Bihar1D(fem, grid)
    [quadL, w] = quadpts1(2 * (fem.locDof - 3));
    NT   = numel(grid) - 1;
    ndof = fem.nDof(NT);
    h = grid(2:end) - grid(1:end-1);

    D2ref = fem.diff2Span(quadL, 1);
    baseK = (D2ref.' * (D2ref .* w));

    K = sparse(ndof,ndof);
    for t = 1:NT
        idx = fem.dofMap(t);
        K(idx,idx) = K(idx,idx) + baseK / (h(t)^3);
    end
end