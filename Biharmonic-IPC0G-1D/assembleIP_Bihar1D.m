function P = assembleIP_Bihar1D(fem, grid, sigma, beta)
    d2phi0 = fem.diff2Span([1,0], 1);
    d2phi1 = fem.diff2Span([0,1], 1);
    dphi0 = fem.diffSpan([1,0], 1);
    dphi1 = fem.diffSpan([0,1], 1);

    k = fem.ord;
    NT = length(grid) - 1;
    h = grid(2:end) - grid(1:end-1);
    P = sparse(NT*k+1, NT*k+1);
    blank = zeros(1, k);

    for eid = 1 : NT-1
        hmean = (h(eid) + h(eid+1)) / 2;

        dphi_left = [dphi1, blank] / h(eid);
        dphi_right = [blank, dphi0] / h(eid+1);
        d2phi_left = [d2phi1, blank] / (h(eid)^2);
        d2phi_right = [blank, d2phi0] / (h(eid+1)^2);

        dphi_jump = dphi_right - dphi_left;
        d2phi_mean = 0.5 * (d2phi_left + d2phi_right);

        idx = (eid-1)*k + (1 : 2*k+1);
        P(idx, idx) = P(idx, idx) + (dphi_jump' * d2phi_mean) ...
                                  + beta * (d2phi_mean' * dphi_jump) ...
                  + (sigma/hmean) * (dphi_jump' * dphi_jump);
    end
end