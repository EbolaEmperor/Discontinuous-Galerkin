function plotSolVec(fem, node, elem, u, elem2edge)

showmesh(node, elem); 
hold on; 
axis equal;

NT = size(elem,1);
cvs = zeros(NT, 4);

for t = 1:NT
    vid = elem(t,:);
    p = node(vid,:);
    center = mean(p,1);

    idxE = elem2edge(t,:);
    u_t = u(idxE);
    
    phi_vals = fem.computeBasisValue_all(t, ones(1,3)/3);
    u_center = phi_vals * u_t(:);
    cvs(t, :) = [center(1), center(2), u_center(1), u_center(2)];
end

quiver(cvs(:,1), cvs(:,2), cvs(:,3), cvs(:,4) , 1, 'b', 'LineWidth', 1);
hold off

end 