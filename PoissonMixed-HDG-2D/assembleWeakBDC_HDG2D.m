function F = assembleWeakBDC_HDG2D(fem1, fem2, node, elem, elem2dof1, elem2dof2, edge, edge2side, alpha, gD)

    NE    = size(edge, 1);
    nDof1 = max(elem2dof1(:));
    nDof2 = max(elem2dof2(:));
    F1 = zeros(nDof1, 1);
    F2 = zeros(nDof2, 1);

    for e = 1:NE
        tid1 = edge2side(e, 1);
        tid2 = edge2side(e, 2);
        if tid1 ~= 0 && tid2 ~= 0
            continue;
        end

        tid   = tid1 + tid2;
        idx1  = elem2dof1(tid,:);
        idx2  = elem2dof2(tid,:);

        vIds  = edge(e,:);
        ep    = node(vIds,:);
        uxy   = ep(1,:); vxy = ep(2,:);
        evec  = vxy - uxy;
        he    = norm(evec);
        n     = [evec(2), -evec(1)] / he;

        tri   = elem(tid,:);
        loc   = [find(tri == vIds(1)), find(tri == vIds(2))];
        wid   = 6 - sum(loc);
        wpt   = node(tri(wid), :);
        if cross2d(evec, wpt - uxy) < 0
            n = -n;
        end

        [quadL, w] = quadpts1(fem1.ord + fem2.ord);
        nq = numel(w);

        for i = 1:nq
            elam = quadL(i,:);
            lam  = zeros(1,3); lam(loc) = elam;

            sigma_in = fem1.computeBasisValue_all(tid, lam);
            nsigma   = n * sigma_in;
            u_in     = fem2.computeBasisValue_all(tid, lam);

            pxy  = lam * node(tri, :);
            gval = gD([pxy(1), pxy(2)]);

            wgt  = w(i) * he;
            F1(idx1) = F1(idx1) + wgt * (nsigma') * gval;
            F2(idx2) = F2(idx2) + wgt * alpha * (u_in') * gval;
        end
    end

    F = [F1; F2];
end
