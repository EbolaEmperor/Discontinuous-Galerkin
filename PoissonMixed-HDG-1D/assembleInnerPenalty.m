function P = assembleInnerPenalty(fem1, fem2, grid, sigma)

    u0 = fem1.getSpan([1,0]);
    u1 = fem1.getSpan([0,1]);
    sigma0 = fem.getSpan([1,0]);
    sigma1 = fem.getSpan([0,1]);

    nDof1 = fem1.locDof;
    nDof2 = fem2.locDof;
    NT = length(grid) - 1;
    P11 = sparse(NT*nDof1, NT*nDof1);
    P12 = sparse(NT*nDof1, NT*nDof2);
    P21 = sparse(NT*nDof2, NT*nDof1);
    P22 = sparse(NT*nDof2, NT*nDof2);
    blank1 = zeros(1, nDof1);
    blank2 = zeros(1, nDof2);

    h = grid(2:end) - grid(1:end-1);

    for eid = 1 : NT-1
        hmean = (h(eid) + h(eid+1)) / 2;

        % phi_left = [phi1, blank];
        % phi_right = [blank, phi0];
        % dphi_left = [dphi1, blank] / h(eid);
        % dphi_right = [blank, dphi0] / h(eid+1);
        % 
        % dphi_jump = phi_right - phi_left;
        % d2phi_mean = 0.5 * (dphi_left + dphi_right);
        % 
        % idx = (eid-1)*nDof + (1 : 2*nDof);
        % P(idx, idx) = P(idx, idx) + (dphi_jump' * d2phi_mean) ...
        %                           + (d2phi_mean' * dphi_jump) ...
        %           + (sigma/hmean) * (dphi_jump' * dphi_jump);
    end

end