function B = assembleDivMass(fem1, fem2, node, elem, elem2dof1, elem2dof2)
    NT    = size(elem,1);
    nDof1  = max(elem2dof1(:));
    nDof2  = max(elem2dof2(:));
    iiC  = cell(NT,1);
    jjC  = cell(NT,1);
    ssBC = cell(NT,1);

    for t = 1:NT
        idx1 = elem2dof1(t,:);
        idx2 = elem2dof2(t,:);
        p    = node(elem(t,:),:);
        Bt   = assembleDivMassLocal(fem1, fem2, p, t);
        [I,J] = ndgrid(idx1, idx2);
        iiC{t}  = I(:);
        jjC{t}  = J(:);
        ssBC{t} = Bt(:);
    end

    ii  = vertcat(iiC{:});
    jj  = vertcat(jjC{:});
    ssB = vertcat(ssBC{:});
    B = sparse(ii, jj, ssB, nDof1, nDof2);
end


function B_t = assembleDivMassLocal(fem1, fem2, p, tid)
    area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
    [quadL, w] = quadpts(fem1.ord + fem2.ord - 1);
    nq = numel(w);
    B_t = zeros(fem1.locDof, fem2.locDof);
    for q = 1:nq
        lambda = quadL(q,:);
        psi_div = fem1.computeBasisDiv_all(tid, lambda);
        phi_val = fem2.computeBasisValue_all(tid, lambda);
        B_t = B_t + w(q) * area * (psi_div' * phi_val);
    end
end