function A = assembleMass(fem, node, elem, elem2dof)
    NT = size(elem,1);
    n  = fem.locDof;
    [Ir, Jc] = ndgrid(1:n, 1:n);
    Ir = Ir(:).';
    Jc = Jc(:).';
    
    Gi = zeros(n, n, NT);

    [quadL, w] = fem.quad2d();
    nq = numel(w);
    for i = 1:nq
        lam = quadL(i,:);
        phi = fem.computeBasisValue_all(1, lam);
        Ap = repmat(phi, 1, 1, NT);
        Gi = Gi + pagemtimes(permute(Ap,[2 1 3]), Ap) * w(i);
    end

    Gi = Gi .* reshape(fem.area, 1, 1, NT);
    rowsMat = elem2dof(:, Ir);
    colsMat = elem2dof(:, Jc);
    ii = reshape(rowsMat.', [], 1);
    jj = reshape(colsMat.', [], 1);

    if nargin < 5, nDof = max(elem2dof(:)); end
    A = sparse(ii, jj, Gi(:), nDof, nDof);
end

% function At = assembleLocalMass(fem, p, tid)
%     area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
%     nDof = fem.locDof;
%     At = zeros(nDof, nDof);
%     [quadL, w] = fem.quad2d();
%     nq = numel(w);
% 
%     for i = 1 : nq
%         lam = quadL(i,:);
%         phi = fem.computeBasisValue_all(tid, lam);
%         At = At + phi' * (phi * (w(i) * area));
%     end
% end