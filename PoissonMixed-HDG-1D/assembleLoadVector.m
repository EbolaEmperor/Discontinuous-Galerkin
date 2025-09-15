function F = assembleLoadVector(fem, grid, f)
    NT   = length(grid) - 1;
    nDof = fem.locDof;

    [quadL, w] = quadpts1(2 * (nDof - 1));
    w  = w(:);
    nq = numel(w);

    phi = zeros(nDof, nq);
    for i = 1:nq
        phi(:, i) = fem.getSpan(quadL(i, :)).';
    end

    h = diff(grid);
    v1v2 = [grid(1:NT) ; grid(2:NT+1)];
    pts  = quadL * v1v2;

    try
        fq = f(pts);
    catch
        fq = arrayfun(f, pts);
    end

    weighted = bsxfun(@times, fq, w);
    Fmat = phi * weighted;
    Fmat = bsxfun(@times, Fmat, h);
    F = Fmat(:);
end
