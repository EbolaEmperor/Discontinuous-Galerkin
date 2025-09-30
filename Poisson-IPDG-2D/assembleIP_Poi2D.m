%% 向量化代码，一次性装配所有内边
function P = assembleIP_Poi2D(fem, node, elem, elem2dof, edge, edge2side, sigma, beta)
    nloc = fem.locDof;
    nDof = max(elem2dof(:));

    tid1 = edge2side(:,1); tid2 = edge2side(:,2);
    mask = (tid1>0) & (tid2>0);
    if ~any(mask), P = sparse(nDof,nDof); return; end
    tid1 = tid1(mask); tid2 = tid2(mask);
    eIdx = find(mask);

    u = node(edge(eIdx,1),:);
    v = node(edge(eIdx,2),:);
    evec = v - u;
    he   = vecnorm(evec,2,2);
    nvec = [evec(:,2), -evec(:,1)] ./ he;

    NEi = numel(tid1);
    K1 = zeros(NEi,2);  K2 = zeros(NEi,2);
    for k = 1:NEi
        tri1 = elem(tid1(k),:);
        tri2 = elem(tid2(k),:);
        uid  = edge(eIdx(k),1);
        vid  = edge(eIdx(k),2);
        K1(k,:) = [find(tri1==uid), find(tri1==vid)];
        K2(k,:) = [find(tri2==uid), find(tri2==vid)];
    end
    mapCode = @(ab) (ab(:,1)==1 & ab(:,2)==2)*1 + ...
                    (ab(:,1)==2 & ab(:,2)==1)*2 + ...
                    (ab(:,1)==2 & ab(:,2)==3)*3 + ...
                    (ab(:,1)==3 & ab(:,2)==2)*4 + ...
                    (ab(:,1)==3 & ab(:,2)==1)*5 + ...
                    (ab(:,1)==1 & ab(:,2)==3)*6;
    code1 = mapCode(K1);  code2 = mapCode(K2);

    D1 = fem.Dlam(:,:,tid1);
    D2 = fem.Dlam(:,:,tid2);

    idx12 = [elem2dof(tid1,:), elem2dof(tid2,:)];
    [Iloc,Jloc] = ndgrid(1:2*nloc, 1:2*nloc);
    Iloc = Iloc(:); Jloc = Jloc(:);
    ii = idx12(:, Iloc).';
    jj = idx12(:, Jloc).';
    ii = ii(:);
    jj = jj(:);
    Pt_all = zeros((2*nloc)^2, NEi);

    [elam, w1d] = fem.quad1d();
    nq = numel(w1d);

    for iq = 1:nq
        t = elam(iq,1); s = elam(iq,2);

        lam12 = mkLam([t s],[1 2]);
        lam21 = mkLam([t s],[2 1]);
        lam23 = mkLam([t s],[2 3]);
        lam32 = mkLam([t s],[3 2]);
        lam31 = mkLam([t s],[3 1]);
        lam13 = mkLam([t s],[1 3]);

        v12 = fem.computeBasisValue_all(1, lam12);
        v21 = fem.computeBasisValue_all(1, lam21);
        v23 = fem.computeBasisValue_all(1, lam23);
        v32 = fem.computeBasisValue_all(1, lam32);
        v31 = fem.computeBasisValue_all(1, lam31);
        v13 = fem.computeBasisValue_all(1, lam13);

        dl12 = fem.computeBasisDlam_all(lam12);
        dl21 = fem.computeBasisDlam_all(lam21);
        dl23 = fem.computeBasisDlam_all(lam23);
        dl32 = fem.computeBasisDlam_all(lam32);
        dl31 = fem.computeBasisDlam_all(lam31);
        dl13 = fem.computeBasisDlam_all(lam13);

        phi   = cat(3, v12, v21, v23, v32, v31, v13);
        dplam = cat(3, dl12,dl21,dl23,dl32,dl31,dl13);

        phiL   = phi(:,:,code1);
        phiR   = phi(:,:,code2);
        dplamL = dplam(:,:,code1);
        dplamR = dplam(:,:,code2);

        gL = pagemtimes(D1, dplamL);
        gR = pagemtimes(D2, dplamR);
        nrow  = permute(nvec,[3 2 1]);
        nphiL = pagemtimes(nrow, gL);
        nphiR = pagemtimes(nrow, gR);

        phi_jump   = cat(2,  phiL, -phiR);
        ngrad_mean = 0.5*cat(2, nphiL,  nphiR);

        wscale = he .* w1d(iq);
        wrow   = permute(wscale, [3 2 1]);

        jt = permute(phi_jump,   [2 1 3]);
        mt = permute(ngrad_mean, [2 1 3]);

        term1 = -beta * pagemtimes(mt, phi_jump);                          % - beta * {n * grad(v)} [u]
        term2 = -pagemtimes(jt, ngrad_mean);                               % - [v] {n * grad(u)}
        term3 = (sigma) * pagemtimes(jt, phi_jump) ./ permute(he,[2 3 1]); % [u] * [v] * (sigma/h)
        Pt_i  = (term1 + term2 + term3) .* wrow;

        Pt_all = Pt_all + reshape(Pt_i, (2*nloc)^2, NEi);
    end

    P = sparse(ii, jj, Pt_all(:), nDof, nDof);
end


%% 逐内边装配代码
% function P = assembleInnerPenalty(fem, node, elem, elem2dof, edge, edge2side, sigma, beta)
%     NE    = size(edge,1);
%     nDof  = max(elem2dof(:));
% 
%     iiC  = cell(NE,1);
%     jjC  = cell(NE,1);
%     ssPC = cell(NE,1);
% 
%     for e = 1 : NE
%         tid1 = edge2side(e, 1);
%         tid2 = edge2side(e, 2);
%         if tid1 == 0 || tid2 == 0, continue; end
%         idx = [elem2dof(tid1,:), elem2dof(tid2,:)];
% 
%         Pt = assembleLocalInnerPenalty(fem, edge(e,:), node(edge(e,:),:),...
%                                        tid1, elem(tid1,:),...
%                                        tid2, elem(tid2,:),...
%                                        sigma, beta);
% 
%         [I,J] = ndgrid(idx,idx);
%         iiC{e}  = I(:);
%         jjC{e}  = J(:);
%         ssPC{e} = Pt(:);
%     end
% 
%     ii  = vertcat(iiC{:});
%     jj  = vertcat(jjC{:});
%     ssP = vertcat(ssPC{:});
%     P = sparse(ii, jj, ssP, nDof, nDof);
% end
% 
% function Pt = assembleLocalInnerPenalty(fem, eid, ep, tid1, elem1, tid2, elem2, sigma, beta)
%     nDof = fem.locDof;
%     Pt = zeros(nDof*2, nDof*2);
%     blank = zeros(1, nDof);
% 
%     u = ep(1,:);
%     v = ep(2,:);
%     e = v - u;
%     he = norm(e);
%     n = [e(2), -e(1)] / he;
% 
%     uid = eid(1);
%     vid = eid(2);
%     K1id = [find(elem1 == uid), find(elem1 == vid)];
%     K2id = [find(elem2 == uid), find(elem2 == vid)];
% 
%     [quadL, w] = fem.quad1d();
%     nq = numel(w);
%     for i = 1 : nq
%         elam = quadL(i,:);
%         lam1 = zeros(1, 3);
%         lam1(K1id) = elam;
%         lam2 = zeros(1, 3);
%         lam2(K2id) = elam;
% 
%         phi_left = [fem.computeBasisValue_all(tid1, lam1), blank];
%         phi_right = [blank, fem.computeBasisValue_all(tid2, lam2)];
% 
%         ngrad_phi_left = [n * fem.computeBasisGrad_all(tid1, lam1), blank];
%         ngrad_phi_right = [blank, n * fem.computeBasisGrad_all(tid2, lam2)];
% 
%         phi_jump = phi_left - phi_right;
%         ngrad_phi_mean = 0.5 * (ngrad_phi_left + ngrad_phi_right);
%         if eid(1) == 1 && eid(2) == 5
%             fuck = w(i) * he * ( ...
%                 - beta * ngrad_phi_mean' * phi_jump ...
%                 - phi_jump' * ngrad_phi_mean ...
%                 + ((sigma/he) * phi_jump') * phi_jump );
%         end
% 
%         Pt = Pt + w(i) * he * ( ...
%                 - beta * ngrad_phi_mean' * phi_jump ...
%                 - phi_jump' * ngrad_phi_mean ...
%                 + ((sigma/he) * phi_jump') * phi_jump );
%     end
% end
% 
