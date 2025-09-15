function P = assembleInnerPenalty(fem1, fem2, node, elem, elem2dof1, elem2dof2, edge, edge2side, alpha)

    NE    = size(edge, 1);
    nDof1 = max(elem2dof1(:));   % for sigma
    nDof2 = max(elem2dof2(:));   % for u
    Ntot  = nDof1 + nDof2;

    iiC = cell(NE,1);
    jjC = cell(NE,1);
    ssC = cell(NE,1);

    for e = 1:NE
        tid1 = edge2side(e, 1);
        tid2 = edge2side(e, 2);

        if tid1 == 0 && tid2 == 0
            continue;
        elseif tid1 == 0 || tid2 == 0
            tid  = tid1 + tid2;                 % 非零的那个
            idx1 = elem2dof1(tid,:);            % 该侧 sigma 的 DOF
            idx2 = elem2dof2(tid,:);            % 该侧 u 的 DOF

            [Pt21, Pt22] = assembleLocalBdryPenalty_D( ...
                fem1, fem2, edge(e,:), node(edge(e,:),:), ...
                tid, elem(tid,:), node(elem(tid,:),:), alpha);

            [i21, j21, s21] = pack_block_(idx2 + nDof1, idx1,         Pt21);
            [i22, j22, s22] = pack_block_(idx2 + nDof1, idx2 + nDof1, Pt22);

            iiC{e} = [i21; i22];
            jjC{e} = [j21; j22];
            ssC{e} = [s21; s22];

        else
            idx1 = [elem2dof1(tid1,:), elem2dof1(tid2,:)];
            idx2 = [elem2dof2(tid1,:), elem2dof2(tid2,:)];

            [Pt11, Pt12, Pt21, Pt22] = ...
                assembleLocalInnerPenalty( ...
                    fem1, fem2, edge(e,:), node(edge(e,:),:), ...
                    tid1, elem(tid1,:), ...
                    tid2, elem(tid2,:), ...
                    alpha);

            [i11, j11, s11] = pack_block_(idx1,          idx1,          Pt11);
            [i12, j12, s12] = pack_block_(idx1,          idx2 + nDof1,  Pt12);
            [i21, j21, s21] = pack_block_(idx2 + nDof1,  idx1,          Pt21);
            [i22, j22, s22] = pack_block_(idx2 + nDof1,  idx2 + nDof1,  Pt22);

            iiC{e} = [i11; i12; i21; i22];
            jjC{e} = [j11; j12; j21; j22];
            ssC{e} = [s11; s12; s21; s22];
        end
    end

    ii = vertcat(iiC{:});
    jj = vertcat(jjC{:});
    ss = vertcat(ssC{:});
    P  = sparse(ii, jj, ss, Ntot, Ntot);
end

function [ii, jj, ss] = pack_block_(rows, cols, M)
    [I, J] = ndgrid(rows, cols);
    ii = I(:);
    jj = J(:);
    ss = M(:);
end


function [Pt11, Pt12, Pt21, Pt22] = ...
    assembleLocalInnerPenalty(fem1, fem2, eid, ep, tid1, elem1, tid2, elem2, alpha)
    
    nDof1 = fem1.locDof;
    nDof2 = fem2.locDof;
    
    Pt11 = zeros(nDof1*2, nDof1*2);
    Pt12 = zeros(nDof1*2, nDof2*2);
    Pt21 = zeros(nDof2*2, nDof1*2);
    Pt22 = zeros(nDof2*2, nDof2*2);

    blank_sigma = zeros(2, nDof1);
    blank_u = zeros(1, nDof2);

    u = ep(1,:);
    v = ep(2,:);
    e = v - u;
    he = norm(e);
    n = [e(2), -e(1)] / he;

    uid = eid(1);
    vid = eid(2);
    K1id = [find(elem1 == uid), find(elem1 == vid)];
    K2id = [find(elem2 == uid), find(elem2 == vid)];

    [quadL, w] = quadpts1(fem1.ord + fem2.ord);
    nq = numel(w);
    for i = 1 : nq
        elam = quadL(i,:);
        lam1 = zeros(1, 3);
        lam1(K1id) = elam;
        lam2 = zeros(1, 3);
        lam2(K2id) = elam;

        sigma_left = [fem1.computeBasisValue_all(tid1, lam1), blank_sigma];
        sigma_right = [blank_sigma, fem1.computeBasisValue_all(tid2, lam2)];
        sigma_jump = n * (sigma_left - sigma_right);
        sigma_mean = 0.5 * (sigma_left + sigma_right);

        u_left = [fem2.computeBasisValue_all(tid1, lam1), blank_u];
        u_right = [blank_u, fem2.computeBasisValue_all(tid2, lam2)];
        u_jump = n' * (u_left - u_right);
        u_mean = 0.5 * (u_left + u_right);

        Pt11 = Pt11 + ( 1/(2*alpha) * w(i) * he * sigma_jump' ) * sigma_jump;
        Pt12 = Pt12 - w(i) * he * sigma_jump' * u_mean;
        Pt21 = Pt21 - w(i) * he * u_jump' * sigma_mean;
        Pt22 = Pt22 + ( (alpha/2) * w(i) * he * u_jump' ) * u_jump;
    end
end



function [Pt21, Pt22] = ...
    assembleLocalBdryPenalty_D(fem1, fem2, eid, ep, tid, elem_t, vtx_t, alpha)

    nDof1 = fem1.locDof;
    nDof2 = fem2.locDof;

    Pt21 = zeros(nDof2, nDof1);
    Pt22 = zeros(nDof2, nDof2);

    u = ep(1,:); v = ep(2,:);
    e  = v - u;
    he = norm(e);
    n  = [e(2), -e(1)] / he;

    % 取第三顶点并判向，确保 n 为外法向（指向单元外）
    uid = eid(1); vid = eid(2);
    Kid = [find(elem_t == uid), find(elem_t == vid)]; % 边在三角的两个顶点编号（1..3）
    wid = 6 - sum(Kid);                               % 1+2+3=6
    wpt = vtx_t(wid, :);
    if cross2d(e, wpt - u) < 0, n = -n; end

    [quadL, w] = quadpts1(fem1.ord + fem2.ord);
    nq = numel(w);
    for i = 1:nq
        elam = quadL(i,:);
        lam  = zeros(1,3); lam(Kid) = elam;

        % sigma_in: 2×nDof1；n*sigma_in -> 1×nDof1
        sigma_in  = fem1.computeBasisValue_all(tid, lam);
        nsigma    = n * sigma_in;              % 1×nDof1

        % u_in: 1×nDof2
        u_in = fem2.computeBasisValue_all(tid, lam);

        Pt21 = Pt21 - (w(i) * he) * (u_in') * nsigma;      % (nDof2×1)*(1×nDof1)
        Pt22 = Pt22 + (w(i) * he) * alpha * (u_in') * u_in;
    end

end