function plotSol(fem, node, elem, u, elem2dof)
    hold on; axis equal
    NT = size(elem,1);

    n = 4;
    [xGrid, yGrid] = meshgrid(linspace(0,1,n));
    idx = (xGrid + yGrid <= 1 + 1e-12);
    xq = xGrid(idx);
    yq = yGrid(idx);
    lam = [1 - xq - yq, xq, yq];
    Nloc = size(lam, 1);

    triRef = delaunay(xq, yq);
    nTri = size(triRef, 1);

    X = node(:,1); Y = node(:,2);
    Ex = X(elem);
    Ey = Y(elem);
    S  = lam.';
    PtsX = Ex * S;
    PtsY = Ey * S;

    locDof = size(elem2dof, 2);
    PhiRef = zeros(Nloc, locDof);
    for i = 1:Nloc
        PhiRef(i,:) = fem.computeBasisValue_all(1, lam(i,:)); 
    end

    U = u(elem2dof.');
    Val = PhiRef * U;
    CData = reshape(Val, [], 1);

    Verts = [reshape(PtsX.', [], 1), reshape(PtsY.', [], 1)];
    Faces = repmat(triRef, NT, 1) + kron((0:NT-1)'*Nloc, ones(nTri,1));

    patch('Faces', Faces, 'Vertices', Verts, ...
          'FaceVertexCData', CData, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    colorbar
    showmesh(node, elem, 'FaceAlpha', 0);
end

%% 逐单元绘制
% function plotSol(fem, node, elem, u, elem2dof)
%     hold on;
%     axis equal
%     NT = size(elem,1);
% 
%     n = 5;
%     [xGrid, yGrid] = meshgrid(linspace(0,1,n));
%     idx = (xGrid + yGrid <= 1 + 1e-12);
%     xGrid = xGrid(idx);
%     yGrid = yGrid(idx);
%     lam = [1-xGrid-yGrid, xGrid, yGrid];
%     Nloc = numel(xGrid);
% 
%     for t = 1:NT
%         idxE = elem2dof(t,:);
%         u_t = u(idxE);
%         val = zeros(Nloc,1);
%         pts = lam * node(elem(t,:), :);
%         tri = delaunay(pts(:,1), pts(:,2));
%         for i = 1 : Nloc
%             phi_vals = fem.computeBasisValue_all(t, lam(i,:));
%             val(i) = phi_vals * u_t(:);
%         end
%         patch('Faces', tri, 'Vertices', pts, ...
%           'FaceVertexCData', val, ...
%           'FaceColor', 'interp', ...
%           'EdgeColor', 'none');
%     end
% 
%     colorbar
%     showmesh(node,elem, 'FaceAlpha', 0);
% end