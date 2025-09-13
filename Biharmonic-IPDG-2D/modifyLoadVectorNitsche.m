function Fb = modifyLoadVectorNitsche(fem, node, elem, elem2dof, edge, edge2side, sigma, beta, gradu)
    NE   = size(edge,1);
    nDof = max(elem2dof(:));
    Fb   = zeros(nDof,1);

    for e = 1:NE
        tid1 = edge2side(e,1);
        tid2 = edge2side(e,2);

        if xor(tid1==0, tid2==0)
            idloc = find([tid1, tid2] ~= 0);
            tid   = edge2side(e, idloc);
            idx   = unique(sort(elem2dof(tid,:)));

            Ft = assembleLocalBdryNitsche( ...
                    fem, edge(e,:), node(edge(e,:),:), idx, ...
                    tid, elem(tid,:), elem2dof(tid,:), node(elem(tid,:),:), ...
                    sigma, beta, gradu);

            Fb(idx) = Fb(idx) + Ft;
        end
    end
end

function Ft = assembleLocalBdryNitsche(fem, eid, ep, dof, tid1, elem1, dof1, vtx, sigma, beta, gradu)
    nDof = numel(dof);
    Ft   = zeros(nDof,1);
    [~, locdof1] = ismember(dof1, dof);

    u  = ep(1,:);
    v  = ep(2,:);
    evec = v - u;
    he = norm(evec);

    n = [evec(2), -evec(1)] / he;
    uid = eid(1); vid = eid(2);
    K1id = [find(elem1==uid), find(elem1==vid)];
    wid = 6 - sum(K1id);
    wpt = vtx(wid,:);
    if cross2d(evec, wpt - u) < 0, n = -n; end

    [quadL, wq] = fem.quad1d();
    nq = numel(wq);

    for i = 1:nq
        elam = quadL(i,:);
        lam1 = zeros(1,3);
        lam1(K1id) = elam;

        p = lam1 * vtx;
        dn_phi   = zeros(1, nDof);
        d2n_phi  = zeros(1, nDof);
        dn_phi(locdof1)  = n * fem.computeBasisGrad_all(tid1, lam1);
        d2n_phi(locdof1) = fem.computeBasisDirectedDiff2_all(tid1, lam1, n);

        g1val = dot(n, gradu(p)); % 只使用了 gradu 的外法向分量
        Ft = Ft + wq(i) * he * ( - beta * d2n_phi.' * g1val + (sigma/he) * dn_phi.' * g1val );
    end
end
