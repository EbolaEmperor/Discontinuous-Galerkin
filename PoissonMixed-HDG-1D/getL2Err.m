function errL2 = getL2Err(fem, grid, uh, u)
    NT   = length(grid) - 1;
    nDof = fem.locDof;

    [quadL, w] = quadpts1(2 * nDof);
    w  = w(:);
    nq = numel(w);

    phi  = zeros(nDof, nq);
    for i = 1:nq
        phi(:, i)  = fem.getSpan(quadL(i, :)).';
    end

    v1v2 = [grid(1:NT); grid(2:NT+1)];
    h    = diff(grid).';
    pts = quadL * v1v2;

    try
        Uex = u(pts).';
    catch
        Uex = arrayfun(u, pts).';
    end
    
    uh_mat = reshape(uh, [nDof, NT]);
    Uh_val = (uh_mat.' * phi);
    
    diff0   = Uh_val - Uex;
    errL2   = sqrt( sum( (diff0.^2) * w .* h ) );
end
