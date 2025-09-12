function F = assembleLoadVector(fem, node, elem, elem2dof, f)

NT = size(elem,1);
nDof = max(elem2dof(:));
nlDof = size(elem2dof, 2);
Floc = zeros(NT,nlDof);
[q_points, w] = fem.quad2d();

for t = 1:NT
    vid = elem(t,:);
    p = node(vid,:);
    area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));

    for q = 1:size(q_points,1)
        lambda = q_points(q,:);
        xy = lambda * p;
        fval = f(xy);
        phi_vals = fem.computeBasisValue_all(t, lambda);
        Floc(t,:) = Floc(t,:) + w(q) * area * (fval * phi_vals);
    end
end

F = accumarray(elem2dof(:), Floc(:), [nDof, 1]);

end 