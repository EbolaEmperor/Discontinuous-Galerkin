function ustar = solveLocalPoisson2D(fem_star, fem1, fem2, node, elem, ...
                                     elem2dof_star, elem2dof1, elem2dof2, ...
                                     sigmah, uh, f)

    NT = size(elem,1);
    nDofStar = max(elem2dof_star(:));
    ustar = zeros(nDofStar,1);

    for t = 1:NT
        vtx = node(elem(t,:), :);
        area = 0.5 * abs(det([vtx(2,:)-vtx(1,:); vtx(3,:)-vtx(1,:)]));

        idxS = elem2dof_star(t,:);
        idx1 = elem2dof1(t,:); 
        idx2 = elem2dof2(t,:);
        ldS  = numel(idxS);

        Kloc = zeros(ldS, ldS);
        Floc = zeros(ldS, 1);
        m_row = zeros(1, ldS);
        m_u   = 0;

        [lam2, w2] = fem_star.quad2d();
        nq2 = numel(w2);
        uh_t = uh(idx2);

        for i = 1:nq2
            lam = lam2(i,:);
            phiS  = fem_star.computeBasisValue_all(t, lam);
            dphiS = fem_star.computeBasisGrad_all(t, lam);

            xq = lam * vtx;
            fq = f(xq);

            wgt = w2(i) * area;
            Kloc = Kloc + (dphiS' * dphiS) * wgt;
            Floc = Floc + (phiS') * fq * wgt;

            m_row = m_row + phiS * wgt;
            phiU  = fem2.computeBasisValue_all(t, lam);
            uval  = phiU * uh_t;
            m_u   = m_u + uval * wgt;
        end

        edgeLoc = [2 3; 3 1; 1 2];
        for ee = 1:3
            a = edgeLoc(ee,1);
            b = edgeLoc(ee,2);
            uxy = vtx(a,:); vxy = vtx(b,:);
            evec = vxy - uxy;
            he = norm(evec);

            n = [evec(2), -evec(1)] / he;
            wid = 6 - (a + b);
            wpt = vtx(wid,:);
            if cross2d(evec, wpt - uxy) < 0
                n = -n;
            end

            [lam1, w1] = fem_star.quad1d();
            nq1 = numel(w1);
            sigma_t = sigmah(idx1);

            for i = 1:nq1
                elam = lam1(i,:);
                lam = zeros(1,3);
                lam([a,b]) = elam;

                sigBasis = fem1.computeBasisValue_all(t, lam);
                nsig_row = n * sigBasis;
                nsig_val = nsig_row * sigma_t;

                phiS = fem_star.computeBasisValue_all(t, lam);

                Floc = Floc + (phiS') * nsig_val * (w1(i) * he);
            end
        end

        Aaug = [Kloc; m_row];
        rhs  = [Floc; m_u];
        ustar(idxS)  = Aaug \ rhs;
    end
end
