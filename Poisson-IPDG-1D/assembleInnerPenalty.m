function P = assembleInnerPenalty(fem, grid, sigma, beta)

    dphi0 = fem.diffSpan([1,0], 1);
    dphi1 = fem.diffSpan([0,1], 1);
    phi0 = fem.getSpan([1,0]);
    phi1 = fem.getSpan([0,1]);

    nDof = fem.locDof;
    NT = length(grid) - 1;
    h = grid(2:end) - grid(1:end-1);
    P = sparse(NT*nDof, NT*nDof);
    blank = zeros(1, nDof);

    for eid = 1 : NT-1
        hmean = (h(eid) + h(eid+1)) / 2;

        phi_left = [phi1, blank];
        phi_right = [blank, phi0];
        dphi_left = [dphi1, blank] / h(eid);
        dphi_right = [blank, dphi0] / h(eid+1);

        dphi_jump = phi_right - phi_left;
        d2phi_mean = 0.5 * (dphi_left + dphi_right);

        idx = (eid-1)*nDof + (1 : 2*nDof);
        P(idx, idx) = P(idx, idx) + (dphi_jump' * d2phi_mean) ...
                                  + beta * (d2phi_mean' * dphi_jump) ...
                  + (sigma/hmean) * (dphi_jump' * dphi_jump);
    end

end