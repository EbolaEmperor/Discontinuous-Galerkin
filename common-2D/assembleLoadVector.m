 % 逐单元装配右端项
function F = assembleLoadVector(fem, node, elem, elem2dof, f)

NT = size(elem,1);
nDof = max(elem2dof(:));
F = zeros(nDof, 1);
[q_points, w] = fem.quad2d();

for t = 1:NT
    vid = elem(t,:);
    p = node(vid,:);
    area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
    dof = elem2dof(t,:);

    for q = 1:size(q_points,1)
        lambda = q_points(q,:);
        xy = lambda * p;
        fval = f(xy);
        phi_vals = fem.computeBasisValue_all(t, lambda);
        F(dof) = F(dof) + w(q) * area * (fval * phi_vals)';
    end
end

end 