function K = assembleStiffness(fem, grid)
    [quadL, w] = quadpts1(2 * (fem.ord-2));
    NT   = numel(grid) - 1;
    k    = fem.ord;
    ndof = NT * k + 1;
    h = grid(2:end) - grid(1:end-1);

    D2ref = fem.diff2Span(quadL, 1);
    baseK = (D2ref.' * (D2ref .* w));

    K = sparse(ndof,ndof);
    for t = 1:NT
        idx = (t-1)*k + (1:k+1);
        K(idx,idx) = K(idx,idx) + baseK / (h(t)^3);
    end
end