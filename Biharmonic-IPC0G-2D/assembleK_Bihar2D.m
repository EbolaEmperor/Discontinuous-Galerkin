function A = assembleK_Bihar2D(fem, node, elem, elem2dof, nDof)
    NT = size(elem,1);
    n  = fem.locDof;

    [Ir, Jc] = ndgrid(1:n, 1:n);
    Ir = Ir(:).';
    Jc = Jc(:).';

    if ~isprop(fem, 'R') || isempty(fem.R)
        R = buildR_from_Dlam(fem.Dlam);    % 3 x 6 x NT
    else
        R = fem.R;                         % 3 x 6 x NT
    end

    Gi = zeros(n, n, NT);

    [quadL, w] = fem.quad2d();
    nq = numel(w);
    for i = 1:nq
        lam = quadL(i,:);
        Hlam = fem.computeBasisHlam_all(lam);
        Hlam_all = repmat(Hlam, 1, 1, NT);
        Hvoigt = pagemtimes(R, Hlam_all);
        Gi = Gi + pagemtimes(permute(Hvoigt, [2 1 3]), Hvoigt) * w(i);
    end

    Gi = Gi .* reshape(fem.area, 1, 1, NT);

    rowsMat = elem2dof(:, Ir);
    colsMat = elem2dof(:, Jc);
    ii = reshape(rowsMat.', [], 1);
    jj = reshape(colsMat.', [], 1);

    if nargin < 5, nDof = max(elem2dof(:)); end
    A = sparse(ii, jj, Gi(:), nDof, nDof);
end

%% 逐单元装配
% function A = assembleK_Bihar2D(fem, node, elem, elem2dof)
%     NT    = size(elem,1);
%     nDof  = max(elem2dof(:));
% 
%     iiC  = cell(NT,1);
%     jjC  = cell(NT,1);
%     ssAC = cell(NT,1);
% 
%     for t = 1 : NT
%         idx = elem2dof(t,:);
%         p = node(elem(t,:),:);
%         At = assembleLocalStiffness(fem, p, t);
% 
%         [I,J] = ndgrid(idx,idx);
%         iiC{t}  = I(:);
%         jjC{t}  = J(:);
%         ssAC{t} = At(:);
%     end
% 
%     ii  = vertcat(iiC{:});
%     jj  = vertcat(jjC{:});
%     ssA = vertcat(ssAC{:});
%     A = sparse(ii,jj,ssA,nDof,nDof);
% end
% 
% function At = assembleLocalStiffness(fem, p, tid)
%     area = 0.5 * abs(det([p(2,:)-p(1,:); p(3,:)-p(1,:)]));
%     nDof = fem.locDof;
%     At = zeros(nDof, nDof);
%     [quadL, w] = fem.quad2d();
%     nq = numel(w);
% 
%     for i = 1 : nq
%         lam = quadL(i,:);
%         hess_phi = fem.computeBasisHessian_all(tid, lam);
%         hess_phi(3,:) = hess_phi(3,:) * sqrt(2);
%         At = At + hess_phi' * (hess_phi * (w(i) * area));
%     end
% end
