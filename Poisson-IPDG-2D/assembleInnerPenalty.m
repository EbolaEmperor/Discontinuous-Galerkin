function P = assembleInnerPenalty(fem, node, elem, elem2dof, edge, edge2side, sigma)
    NE    = size(edge,1);
    nDof  = max(elem2dof(:));

    iiC  = cell(NE,1);
    jjC  = cell(NE,1);
    ssPC = cell(NE,1);

    for e = 1 : NE
        tid1 = edge2side(e, 1);
        tid2 = edge2side(e, 2);
        if tid1 == 0 || tid2 == 0, continue; end
        idx = [elem2dof(tid1,:), elem2dof(tid2,:)];

        Pt = assembleLocalInnerPenalty(fem, edge(e,:), node(edge(e,:),:),...
                                       tid1, elem(tid1,:),...
                                       tid2, elem(tid2,:),...
                                       sigma);

        [I,J] = ndgrid(idx,idx);
        iiC{e}  = I(:);
        jjC{e}  = J(:);
        ssPC{e} = Pt(:);
    end

    ii  = vertcat(iiC{:});
    jj  = vertcat(jjC{:});
    ssP = vertcat(ssPC{:});
    P = sparse(ii, jj, ssP, nDof, nDof);
end

function Pt = assembleLocalInnerPenalty(fem, eid, ep, tid1, elem1, tid2, elem2, sigma)
    nDof = fem.locDof;
    Pt = zeros(nDof*2, nDof*2);
    blank = zeros(1, nDof);

    u = ep(1,:);
    v = ep(2,:);
    e = v - u;
    he = norm(e);
    n = [e(2), -e(1)] / he;

    uid = eid(1);
    vid = eid(2);
    K1id = [find(elem1 == uid), find(elem1 == vid)];
    K2id = [find(elem2 == uid), find(elem2 == vid)];

    [quadL, w] = fem.quad1d();
    nq = numel(w);
    for i = 1 : nq
        elam = quadL(i,:);
        lam1 = zeros(1, 3);
        lam1(K1id) = elam;
        lam2 = zeros(1, 3);
        lam2(K2id) = elam;

        phi_left = [fem.computeBasisValue_all(tid1, lam1), blank];
        phi_right = [blank, fem.computeBasisValue_all(tid2, lam2)];

        ngrad_phi_left = [n * fem.computeBasisGrad_all(tid1, lam1), blank];
        ngrad_phi_right = [blank, n * fem.computeBasisGrad_all(tid2, lam2)];

        phi_jump = phi_left - phi_right;
        ngrad_phi_mean = 0.5 * (ngrad_phi_left + ngrad_phi_right);

        Pt = Pt + w(i) * he * ( ...
                - ngrad_phi_mean' * phi_jump ...
                - phi_jump' * ngrad_phi_mean ...
                + ((sigma/he) * phi_jump') * phi_jump );
    end
end