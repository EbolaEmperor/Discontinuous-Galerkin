function A = assembleStiffness(fem, node, elem, elem2dof)
    NT    = size(elem,1);
    nDof  = max(elem2dof(:));

    iiC  = cell(NT,1);
    jjC  = cell(NT,1);
    ssAC = cell(NT,1);

    for t = 1 : NT
        idx = elem2dof(t,:);
        p = node(elem(t,:),:);
        At = assembleLocalStiffness(fem, p, t);

        [I,J] = ndgrid(idx,idx);
        iiC{t}  = I(:);
        jjC{t}  = J(:);
        ssAC{t} = At(:);
    end

    ii  = vertcat(iiC{:});
    jj  = vertcat(jjC{:});
    ssA = vertcat(ssAC{:});
    A = sparse(ii,jj,ssA,nDof,nDof);
end

function At = assembleLocalStiffness(fem, p, tid)
    area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
    nDof = fem.locDof;
    At = zeros(nDof, nDof);
    [quadL, w] = fem.quad2d();
    nq = numel(w);

    for i = 1 : nq
        lam = quadL(i,:);
        grad_phi = fem.computeBasisGrad_all(tid, lam);
        At = At + grad_phi' * (grad_phi * (w(i) * area));
    end
end