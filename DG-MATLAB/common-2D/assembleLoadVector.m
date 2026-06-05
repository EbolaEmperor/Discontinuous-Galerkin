%% 向量化装配右端项，一次性装配所有单元
function F = assembleLoadVector(fem, node, elem, elem2dof, f)
    NT   = size(elem,1);
    nloc = fem.locDof;
    nDof = max(elem2dof(:));

    area = fem.area(:);
    P = permute(reshape(node(elem.', :), 3, NT, 2), [1 3 2]);

    [q_points, w] = fem.quad2d();
    nq = numel(w);

    Phi_ref = zeros(nq, nloc);
    for q = 1:nq
        Phi_ref(q, :) = fem.computeBasisValue_all(1, q_points(q,:));
    end

    PT = pagemtimes(q_points, P);
    PT = permute(PT, [3 2 1]);

    PT_row  = permute(PT, [1 3 2]);
    pts_all = reshape(PT_row, NT*nq, 2);
    fvals   = f(pts_all);
    fvals   = reshape(fvals, NT, nq);

    G = fvals .* (w(:)).';

    B = G * Phi_ref;
    B = B .* area;
    F = accumarray(elem2dof(:), B(:), [nDof, 1]);
end

%% 逐单元装配右端项
% function F = assembleLoadVector(fem, node, elem, elem2dof, f)
% 
% NT = size(elem,1);
% nDof = max(elem2dof(:));
% F = zeros(nDof, 1);
% [q_points, w] = fem.quad2d();
% 
% for t = 1:NT
%     vid = elem(t,:);
%     p = node(vid,:);
%     area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
%     dof = elem2dof(t,:);
% 
%     for q = 1:size(q_points,1)
%         lambda = q_points(q,:);
%         xy = lambda * p;
%         fval = f(xy);
%         phi_vals = fem.computeBasisValue_all(t, lambda);
%         F(dof) = F(dof) + w(q) * area * (fval * phi_vals)';
%     end
% end
% 
% end 