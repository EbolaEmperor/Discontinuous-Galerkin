%% 向量化装配内边的 Inner Penalty，再用原来的代码逐单元装配边界的 Penalty
function P = assembleIP_Bihar2D(fem, node, elem, elem2dof, edge, edge2side, sigma, beta)
    nloc = fem.locDof;
    nDof = max(elem2dof(:));

    if ~isprop(fem,'R') || isempty(fem.R)
        Rall = buildR_from_Dlam(fem.Dlam);     % 3 x 6 x NT
    else
        Rall = fem.R;
    end

    tid1 = edge2side(:,1); tid2 = edge2side(:,2);
    maskI = (tid1>0) & (tid2>0);
    P_inner = sparse(nDof, nDof);
    if any(maskI)
        tid1I = tid1(maskI); 
        tid2I = tid2(maskI);
        eIdxI = find(maskI);
        NEi   = numel(tid1I);

        uI    = node(edge(eIdxI,1),:);
        vI    = node(edge(eIdxI,2),:);
        evecI = vI - uI;
        heI   = vecnorm(evecI,2,2);
        nvecI = [evecI(:,2), -evecI(:,1)] ./ heI;

        % 边在两侧三角形里的局部顶点编号（用于 bary 置入）
        K1 = zeros(NEi,2); K2 = zeros(NEi,2);
        for k = 1:NEi
            tri1 = elem(tid1I(k),:);
            tri2 = elem(tid2I(k),:);
            uid  = edge(eIdxI(k),1);
            vid  = edge(eIdxI(k),2);
            K1(k,:) = [find(tri1==uid), find(tri1==vid)];
            K2(k,:) = [find(tri2==uid), find(tri2==vid)];
        end
        mapCode = @(ab) (ab(:,1)==1 & ab(:,2)==2)*1 + ...
                        (ab(:,1)==2 & ab(:,2)==1)*2 + ...
                        (ab(:,1)==2 & ab(:,2)==3)*3 + ...
                        (ab(:,1)==3 & ab(:,2)==2)*4 + ...
                        (ab(:,1)==3 & ab(:,2)==1)*5 + ...
                        (ab(:,1)==1 & ab(:,2)==3)*6;
        code1 = mapCode(K1);
        code2 = mapCode(K2);

        D1 = fem.Dlam(:,:,tid1I);           % 2 x 3 x NEi
        D2 = fem.Dlam(:,:,tid2I);           % 2 x 3 x NEi
        R1 = Rall(:,:,tid1I);               % 3 x 6 x NEi
        R2 = Rall(:,:,tid2I);               % 3 x 6 x NEi

        % 全局索引（2*nloc block）
        idx12 = [elem2dof(tid1I,:), elem2dof(tid2I,:)];   % NEi x (2*nloc)
        [Iloc,Jloc] = ndgrid(1:2*nloc, 1:2*nloc);
        Iloc = Iloc(:); Jloc = Jloc(:);
        iiI = idx12(:, Iloc).';  iiI = iiI(:);
        jjI = idx12(:, Jloc).';  jjI = jjI(:);

        PtI_all = zeros((2*nloc)^2, NEi);

        % 1D 求积
        [elam, w1d] = fem.quad1d();
        nq = numel(w1d);

        % 方向二阶导权重（Voigt 第三分量含 sqrt(2)）
        nxI = nvecI(:,1); nyI = nvecI(:,2);
        srowI = permute([nxI.^2, nyI.^2, sqrt(2)*nxI.*nyI], [3 2 1]); % 1x3xNEi

        % 法向一阶导权重
        nrowI = permute(nvecI, [3 2 1]);   % 1x2xNEi

        for iq = 1:nq
            t = elam(iq,1); s = elam(iq,2);

            % 六种边置入
            lam12 = mkLam([t s],[1 2]);
            lam21 = mkLam([t s],[2 1]);
            lam23 = mkLam([t s],[2 3]);
            lam32 = mkLam([t s],[3 2]);
            lam31 = mkLam([t s],[3 1]);
            lam13 = mkLam([t s],[1 3]);

            % 一阶（对 λ）与二阶（bary Hessian 的 6 分量）
            dl12 = fem.computeBasisDlam_all(lam12);
            dl21 = fem.computeBasisDlam_all(lam21);
            dl23 = fem.computeBasisDlam_all(lam23);
            dl32 = fem.computeBasisDlam_all(lam32);
            dl31 = fem.computeBasisDlam_all(lam31);
            dl13 = fem.computeBasisDlam_all(lam13);

            H12  = fem.computeBasisHlam_all(lam12);
            H21  = fem.computeBasisHlam_all(lam21);
            H23  = fem.computeBasisHlam_all(lam23);
            H32  = fem.computeBasisHlam_all(lam32);
            H31  = fem.computeBasisHlam_all(lam31);
            H13  = fem.computeBasisHlam_all(lam13);

            dplam = cat(3, dl12,dl21,dl23,dl32,dl31,dl13);   % 3 x nloc x 6
            Hlam  = cat(3, H12, H21, H23, H32, H31, H13);    % 6 x nloc x 6

            % 选形
            dplamL = dplam(:,:,code1);      % 3 x nloc x NEi
            dplamR = dplam(:,:,code2);      % 3 x nloc x NEi
            HlamL  = Hlam(:,:,code1);       % 6 x nloc x NEi
            HlamR  = Hlam(:,:,code2);       % 6 x nloc x NEi

            % n·grad
            gL = pagemtimes(D1, dplamL);           % 2 x nloc x NEi
            gR = pagemtimes(D2, dplamR);           % 2 x nloc x NEi
            nphiL = pagemtimes(nrowI, gL);         % 1 x nloc x NEi
            nphiR = pagemtimes(nrowI, gR);         % 1 x nloc x NEi

            % n^T H n
            HvoigtL = pagemtimes(R1, HlamL);       % 3 x nloc x NEi
            HvoigtR = pagemtimes(R2, HlamR);       % 3 x nloc x NEi
            d2nL    = pagemtimes(srowI, HvoigtL);  % 1 x nloc x NEi
            d2nR    = pagemtimes(srowI, HvoigtR);  % 1 x nloc x NEi

            % [∂_n] 与 {∂_n^2}
            dn_jump  = cat(2, nphiL, -nphiR);      % 1 x (2nloc) x NEi
            d2n_mean = 0.5 * cat(2, d2nL, d2nR);   % 1 x (2nloc) x NEi

            wscale = heI .* w1d(iq);
            wrow   = permute(wscale, [3 2 1]);     % 1 x 1 x NEi

            jt = permute(dn_jump,  [2 1 3]);       % (2nloc) x 1 x NEi
            mt = permute(d2n_mean, [2 1 3]);       % (2nloc) x 1 x NEi

            term1 = -beta  * pagemtimes(mt, dn_jump);                          % - beta * {∂_n^2 v} [∂_n u]
            term2 = -       pagemtimes(jt, d2n_mean);                          % - [∂_n v] {∂_n^2 u}
            term3 = (sigma) * pagemtimes(jt, dn_jump) ./ permute(heI,[2 3 1]); % (σ/h) [∂_n v][∂_n u]

            Pt_i  = (term1 + term2 + term3) .* wrow;                           % (2nloc)x(2nloc)xNEi
            PtI_all = PtI_all + reshape(Pt_i, (2*nloc)^2, NEi);
        end

        P_inner = sparse(iiI, jjI, PtI_all(:), nDof, nDof);
    end

    % 合并
    P = P_inner + assembleIP_Bihar2D_bdry(fem, node, elem, elem2dof, edge, edge2side, sigma, beta);
end


function P = assembleIP_Bihar2D_bdry(fem, node, elem, elem2dof, edge, edge2side, sigma, beta)
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
                                          idloc, sigma, beta);
        else
            continue
            % 废弃原本的逐单元装配内边 Inner Penalty 的步骤
            % idx = [elem2dof(tid1,:), elem2dof(tid2,:)];
            % idx = unique(sort(idx));
            % Pt = assembleLocalInnerPenalty(fem, edge(e,:), node(edge(e,:),:), idx,...
            %                                tid1, elem(tid1,:), elem2dof(tid1,:),...
            %                                tid2, elem(tid2,:), elem2dof(tid2,:),...
            %                                sigma, beta);
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

function Pt = assembleLocalInnerPenalty(fem, eid, ep, dof, tid1, elem1, dof1, tid2, elem2, dof2, sigma, beta)
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
                - beta * d2n_phi_mean' * dn_phi_jump ...
                - dn_phi_jump' * d2n_phi_mean ...
                + ((sigma/he) * dn_phi_jump') * dn_phi_jump );
    end
end

function Pt = assembleLocalBdryPenalty(fem, eid, ep, dof, tid1, elem1, dof1, vtx, id, sigma, beta)
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
                - beta * d2n_phi_mean' * dn_phi_jump ...
                - dn_phi_jump' * d2n_phi_mean ...
                + ((sigma/he) * dn_phi_jump') * dn_phi_jump );
    end
end
