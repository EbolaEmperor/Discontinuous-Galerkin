%% 所有单元一起向量化计算
function errL2 = getL2Err(fem, node, elem, elem2dof, u, u_exact)

    NT   = size(elem,1);
    Uloc = u(elem2dof);
    
    area = fem.area(:);
    Ptmp = node(elem(:), :);
    P    = permute(reshape(Ptmp, NT, 3, 2), [2 3 1]);
    
    [q_points, w] = fem.quad2d();
    nq = numel(w);
    errL2 = 0.0;
    
    for iq = 1:nq
        lam = q_points(iq, :);
    
        phi_ref = fem.computeBasisValue_all(1, lam);
        uh = Uloc * phi_ref.';
    
        pt = pagemtimes(lam, P);
        pts = permute(pt, [3 2 1]);
    
        u_ex  = u_exact(pts);
        wphys = w(iq) * area;
        errL2 = errL2 + sum((uh - u_ex).^2 .* wphys);
    end
    
    errL2 = sqrt(sum(errL2));

end

%% 逐单元计算
% function errL2 = getL2Err(fem, node, elem, elem2dof, u, u_exact)
% 
% NT = size(elem,1);
% errL2 = 0;
% [q_points, w] = fem.quad2d();
% 
% for t = 1:NT
%     p = node(elem(t,:),:);
%     idx = elem2dof(t,:);
%     u_t = u(idx);
% 
%     area = 0.5*abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
% 
%     for i = 1:size(q_points,1)
%         lambda = q_points(i,:);
%         pt = lambda * p;
% 
%         phi_vals = fem.computeBasisValue_all(t, lambda);
%         uh      = u_t(:)' * phi_vals';
%         u_ex = u_exact(pt);
% 
%         errL2 = errL2 + w(i) * area * sum((uh - u_ex).^2);
%     end
% end
% 
% errL2 = sqrt(errL2);
% 
% end 