function [c, freeDof] = interpStrongBDC(fem, node, elem, elem2dof, domain, u_exact)
    ord = fem.ord;
    locDof = fem.locDof;
    pnt = zeros(locDof, 3);
    ind = 0;
    for i = 0 : ord
        for j = 0 : ord - i
            ind = ind + 1;
            pnt(ind, :) = [1-(i+j)/ord, i/ord, j/ord];
        end
    end

    NT = size(elem, 1);
    nDof = max(elem2dof(:));
    c = zeros(nDof,1);
    isFreeDof = ones(nDof, 1);

    for t = 1 : NT
        vtx = node(elem(t,:),:);
        lagNodes = pnt * vtx;
        isBdNode = domain.findBdryNodes(lagNodes);
        bdIdx = find(isBdNode);
        for id = bdIdx(:)'
            globalID = elem2dof(t,id);
            isFreeDof(globalID) = 0;
            c(globalID) = u_exact(lagNodes(id,:));
        end
    end
    freeDof = find(isFreeDof);
end