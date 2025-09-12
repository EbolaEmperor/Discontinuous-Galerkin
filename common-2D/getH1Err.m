function [errH1, errL2, errSemiH1] = getH1Err(fem, node, elem, elem2dof, u, u_exact, grad_u_exact)

NT = size(elem,1);
errL2 = 0;
errSemiH1 = 0;
[q_points, w] = fem.quad2d();

for t = 1:NT
    p = node(elem(t,:),:);
    idx = elem2dof(t,:);
    u_t = u(idx);
    
    area = 0.5*abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
    
    for i = 1:size(q_points,1)
        lambda = q_points(i,:);
        pt = lambda * p;

        phi_vals = fem.computeBasisValue_all(t, lambda);
        phi_grad = fem.computeBasisGrad_all(t, lambda);
        
        uh      = u_t(:)' * phi_vals';
        uh_grad = u_t(:)' * phi_grad';

        u_ex = u_exact(pt);
        u_grad_ex = grad_u_exact(pt);

        errL2 = errL2 + w(i) * area * sum((uh - u_ex).^2);
        errSemiH1 = errSemiH1 + w(i) * area * sum((uh_grad - u_grad_ex).^2);
    end
end

errH1 = sqrt(errL2 + errSemiH1);
errL2 = sqrt(errL2);
errSemiH1 = sqrt(errSemiH1);

end 