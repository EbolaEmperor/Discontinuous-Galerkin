%% 所有单元一起向量化计算
function [errH1, errL2, errSemiH1] = getH1Err(fem, node, elem, elem2dof, u, u_exact, grad_u_exact)

    NT   = size(elem,1);
    Uloc = u(elem2dof);
    
    area = fem.area(:);
    Dlam = fem.Dlam;
    
    Ptmp = node(elem(:), :);
    P    = permute(reshape(Ptmp, NT, 3, 2), [2 3 1]);
    
    [q_points, w] = fem.quad2d();
    nq = numel(w);
    
    errL2 = 0.0;
    errSemiH1 = 0.0;
    
    for iq = 1:nq
        lam = q_points(iq, :);
    
        phi_ref   = fem.computeBasisValue_all(1, lam);
        dphi_dlam = fem.computeBasisDlam_all(lam);
    
        dphi_rep = repmat(dphi_dlam, 1, 1, NT);
        grad_phi = pagemtimes(Dlam, dphi_rep);
    
        uh = Uloc * phi_ref.';
        Up = permute(Uloc, [3 2 1]);
        Gp = permute(grad_phi, [2 1 3]);
        uhg3 = pagemtimes(Up, Gp);
        uh_grad = permute(uhg3, [3 2 1]);
    
        pt = pagemtimes(lam, P);
        pts = permute(pt, [3 2 1]);
    
        u_ex      = u_exact(pts);
        u_grad_ex = grad_u_exact(pts);
    
        wphys = w(iq) * area;
        errL2     = errL2     + sum((uh - u_ex).^2 .* wphys);
        errSemiH1 = errSemiH1 + sum(sum((uh_grad - u_grad_ex).^2, 2) .* wphys);
    end
    
    errL2     = sqrt(errL2);
    errSemiH1 = sqrt(errSemiH1);
    errH1     = sqrt(errL2^2 + errSemiH1^2);

end

%% 逐单元计算
% function [errH1, errL2, errSemiH1] = getH1Err(fem, node, elem, elem2dof, u, u_exact, grad_u_exact)
% 
% NT = size(elem,1);
% errL2 = 0;
% errSemiH1 = 0;
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
%         phi_grad = fem.computeBasisGrad_all(t, lambda);
% 
%         uh      = u_t(:)' * phi_vals';
%         uh_grad = u_t(:)' * phi_grad';
% 
%         u_ex = u_exact(pt);
%         u_grad_ex = grad_u_exact(pt);
% 
%         errL2 = errL2 + w(i) * area * sum((uh - u_ex).^2);
%         errSemiH1 = errSemiH1 + w(i) * area * sum((uh_grad - u_grad_ex).^2);
%     end
% end
% 
% errH1 = sqrt(errL2 + errSemiH1);
% errL2 = sqrt(errL2);
% errSemiH1 = sqrt(errSemiH1);
% 
% end 