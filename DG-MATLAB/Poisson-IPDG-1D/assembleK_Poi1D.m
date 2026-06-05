function K = assembleK_Poi1D(fem, grid)
    [quadL, w] = quadpts1_my(2 * (fem.locDof - 2));
    NT   = numel(grid) - 1;
    ndof = fem.nDof(NT);
    h = grid(2:end) - grid(1:end-1);

    Dref = fem.diffSpan(quadL, 1);
    baseK = (Dref.' * (Dref .* w));

    K = sparse(ndof,ndof);
    for t = 1:NT
        idx = fem.dofMap(t);
        K(idx,idx) = K(idx,idx) + baseK / h(t);
    end
end