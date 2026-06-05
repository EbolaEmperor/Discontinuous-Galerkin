function F = assembleWeakBDC_HDG1D(fem1, fem2, grid, alpha, bdc, bdval)

    sigma0 = fem1.getSpan([1,0]);
    sigma1 = fem1.getSpan([0,1]);
    u0 = fem2.getSpan([1,0]);
    u1 = fem2.getSpan([0,1]);

    NT = length(grid) - 1;
    nDof1 = fem1.locDof;
    nDof2 = fem2.locDof;
    F1 = zeros(nDof1 * NT, 1);
    F2 = zeros(nDof2 * NT, 1);

    idx1 = 1 : nDof1;
    idx2 = 1 : nDof2;

    if bdc.left == "D"
        F1(idx1) = - sigma0' * bdval(1);
        F2(idx2) = ( alpha * u0' ) * bdval(1);
    else
        F1(idx1) = - ( 1/alpha * sigma0' ) * bdval(1);
        F2(idx2) = u0' * bdval(1);
    end

    idx1 = (NT-1) * nDof1 + (1 : nDof1);
    idx2 = (NT-1) * nDof2 + (1 : nDof2);
    
    if bdc.right == "D"
        F1(idx1) = sigma1' * bdval(2);
        F2(idx2) = ( alpha * u1' ) * bdval(2);
    else
        F1(idx1) = ( 1/alpha * sigma1' ) * bdval(2);
        F2(idx2) = u1' * bdval(2);
    end

    F = [F1; F2];

end