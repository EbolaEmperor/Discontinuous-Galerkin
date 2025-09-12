function P = assembleInnerPenalty(fem, node, elem, elem2dof, edge, edge2side, sigma)
    NE    = size(edge,1);
    nDof  = max(elem2dof(:));

    iiC  = cell(NE,1);
    jjC  = cell(NE,1);
    ssPC = cell(NE,1);

    for e = 1 : NE
        tid1 = edge2side(e, 1);
        tid2 = edge2side(e, 2);
        if tid1 == 0 || tid2 == 0
            idloc = find([tid1, tid2] ~= 0);
            idgl = edge2side(e, idloc);
            idx = elem2dof(idgl,:);
            idx = unique(sort(idx));
            Pt = assembleLocalBdryPenalty(fem, edge(e,:), node(edge(e,:),:), idx,...
                                          idgl, elem(idgl,:), elem2dof(idgl,:), node(elem(idgl,:),:),...
                                          idloc, sigma);
        else
            idx = [elem2dof(tid1,:), elem2dof(tid2,:)];
            idx = unique(sort(idx));
            Pt = assembleLocalInnerPenalty(fem, edge(e,:), node(edge(e,:),:), idx,...
                                           tid1, elem(tid1,:), elem2dof(tid1,:),...
                                           tid2, elem(tid2,:), elem2dof(tid2,:),...
                                           sigma);
        end

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

function Pt = assembleLocalInnerPenalty(fem, eid, ep, dof, tid1, elem1, dof1, tid2, elem2, dof2, sigma)
    nDof = numel(dof);
    Pt = zeros(nDof, nDof);
    [~,locdof1] = ismember(dof1, dof);
    [~,locdof2] = ismember(dof2, dof);

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

        dn_phi_l = zeros(1, nDof);
        dn_phi_r = zeros(1, nDof);
        dn_phi_l(locdof1) = n * fem.computeBasisGrad_all(tid1, lam1);
        dn_phi_r(locdof2) = n * fem.computeBasisGrad_all(tid2, lam2);

        d2n_phi_l = zeros(1, nDof);
        d2n_phi_r = zeros(1, nDof);
        d2n_phi_l(locdof1) = fem.computeBasisDirectedDiff2_all(tid1, lam1, n);
        d2n_phi_r(locdof2) = fem.computeBasisDirectedDiff2_all(tid2, lam2, n);

        dn_phi_jump = dn_phi_l - dn_phi_r;
        d2n_phi_mean = 0.5 * (d2n_phi_l + d2n_phi_r);

        Pt = Pt + w(i) * he * ( ...
                - d2n_phi_mean' * dn_phi_jump ...
                - dn_phi_jump' * d2n_phi_mean ...
                + ((sigma/he) * dn_phi_jump') * dn_phi_jump );
    end
end

function Pt = assembleLocalBdryPenalty(fem, eid, ep, dof, tid1, elem1, dof1, vtx, id, sigma)
    nDof = numel(dof);
    Pt = zeros(nDof, nDof);
    [~,locdof1] = ismember(dof1, dof);

    u = ep(1,:);
    v = ep(2,:);
    e = v - u;
    he = norm(e);
    n = [e(2), -e(1)] / he;
    
    uid = eid(1);
    vid = eid(2);
    K1id = [find(elem1 == uid), find(elem1 == vid)];

    wid = 6 - sum(K1id);
    w = vtx(wid, :);
    if cross2d(e, w - u) < 0, n = -n; end

    [quadL, w] = fem.quad1d();
    nq = numel(w);
    for i = 1 : nq
        elam = quadL(i,:);
        lam1 = zeros(1, 3);
        lam1(K1id) = elam;

        dn_phi_jump = zeros(1, nDof);
        dn_phi_jump(locdof1) = n * fem.computeBasisGrad_all(tid1, lam1);

        d2n_phi_mean = zeros(1, nDof);
        d2n_phi_mean(locdof1) = fem.computeBasisDirectedDiff2_all(tid1, lam1, n);

        Pt = Pt + w(i) * he * ( ...
                - d2n_phi_mean' * dn_phi_jump ...
                - dn_phi_jump' * d2n_phi_mean ...
                + ((sigma/he) * dn_phi_jump') * dn_phi_jump );
    end
end