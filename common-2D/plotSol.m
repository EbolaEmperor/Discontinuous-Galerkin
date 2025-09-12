function plotSol(fem, node, elem, u, elem2dof)
    hold on;
    axis equal
    NT = size(elem,1);

    n = 5;
    [xGrid, yGrid] = meshgrid(linspace(0,1,n));
    idx = (xGrid + yGrid <= 1 + 1e-12);
    xGrid = xGrid(idx);
    yGrid = yGrid(idx);
    lam = [1-xGrid-yGrid, xGrid, yGrid];
    Nloc = numel(xGrid);

    for t = 1:NT
        idxE = elem2dof(t,:);
        u_t = u(idxE);
        val = zeros(Nloc,1);
        pts = lam * node(elem(t,:), :);
        tri = delaunay(pts(:,1), pts(:,2));
        for i = 1 : Nloc
            phi_vals = fem.computeBasisValue_all(t, lam(i,:));
            val(i) = phi_vals * u_t(:);
        end
        patch('Faces', tri, 'Vertices', pts, ...
          'FaceVertexCData', val, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    end
    
    colorbar
    showmesh(node,elem, 'FaceAlpha', 0);
end