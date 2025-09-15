function P = assembleMixedPenalty(fem1, fem2, grid, alpha, bdc)

    if nargin < 5
        bdc.left = "D";
        bdc.right = "D";
    end

    sigma0 = fem1.getSpan([1,0]);
    sigma1 = fem1.getSpan([0,1]);
    u0 = fem2.getSpan([1,0]);
    u1 = fem2.getSpan([0,1]);

    nDof1 = fem1.locDof;
    nDof2 = fem2.locDof;
    NT = length(grid) - 1;
    P11 = sparse(NT*nDof1, NT*nDof1);
    P12 = sparse(NT*nDof1, NT*nDof2);
    P21 = sparse(NT*nDof2, NT*nDof1);
    P22 = sparse(NT*nDof2, NT*nDof2);
    blank1 = zeros(1, nDof1);
    blank2 = zeros(1, nDof2);

    for eid = 1 : NT-1
        sigma_left = [sigma1, blank1];
        sigma_right = [blank1, sigma0];
        u_left = [u1, blank2];
        u_right = [blank2, u0];
        
        u_jump = u_left - u_right;
        u_mean = 0.5 * (u_left + u_right);
        sigma_jump = sigma_left - sigma_right;
        sigma_mean = 0.5 * (sigma_left + sigma_right);

        idx1 = (eid-1) * nDof1 + (1 : 2*nDof1);
        idx2 = (eid-1) * nDof2 + (1 : 2*nDof2);

        P11(idx1, idx1) = P11(idx1, idx1) + ( 1/(2*alpha) * sigma_jump' ) * sigma_jump;
        P12(idx1, idx2) = P12(idx1, idx2) - sigma_jump' * u_mean;
        P21(idx2, idx1) = P21(idx2, idx1) - u_jump' * sigma_mean;
        P22(idx2, idx2) = P22(idx2, idx2) + ( (alpha/2) * u_jump' ) * u_jump;
    end

    % left boundary
    idx1 = 1 : nDof1;
    idx2 = 1 : nDof2;
    if bdc.left == "D"
        P21(idx2, idx1) = P21(idx2, idx1) + u0' * sigma0;
        P22(idx2, idx2) = P22(idx2, idx2) + ( alpha * u0' ) * u0;
    else
        assert(bdc.left == "N");
        P11(idx1, idx1) = P11(idx1, idx1) + ( 1/alpha * sigma0' ) * sigma0;
        P12(idx1, idx2) = P12(idx1, idx2) + sigma0' * u0;
    end

    % right boundary
    idx1 = (NT-1) * nDof1 + (1 : nDof1);
    idx2 = (NT-1) * nDof2 + (1 : nDof2);
    if bdc.right == "D"
        P21(idx2, idx1) = P21(idx2, idx1) - u1' * sigma1;
        P22(idx2, idx2) = P22(idx2, idx2) + ( alpha * u1' ) * u1;
    else
        assert(bdc.right == "N");
        P11(idx1, idx1) = P11(idx1, idx1) + ( 1/alpha * sigma1' ) * sigma1;
        P12(idx1, idx2) = P12(idx1, idx2) - sigma1' * u1;
    end

    P = [P11, P12; P21, P22];

end