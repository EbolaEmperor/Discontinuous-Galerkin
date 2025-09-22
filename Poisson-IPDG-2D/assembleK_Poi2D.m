%% 向量化装配所有单元
function A = assembleK_Poi2D(fem, node, elem, elem2dof)
    NT = size(elem,1);
    E = pagemtimes(permute(fem.Dlam,[2 1 3]), fem.Dlam); 

    n  = fem.locDof;
    [Ir, Jc] = ndgrid(1:n, 1:n);
    Ir = Ir(:);
    Jc = Jc(:);

    offset = repelem((0:NT-1)'*n, n*n, 1);
    ii = repmat(Ir, NT, 1) + offset;
    jj = repmat(Jc, NT, 1) + offset;
    Gi = zeros(n, n, NT);

    [quadL, w] = fem.quad2d();
    nq = numel(w);
    for i = 1 : nq
        lam = quadL(i,:);
        dlamPhi = fem.computeBasisDlam_all(lam);
        Ap = repmat(dlamPhi, 1, 1, NT);
        T  = pagemtimes(E, Ap);
        Gi = Gi + pagemtimes(permute(Ap,[2 1 3]), T) * w(i);
    end

    Gi = Gi .* reshape(fem.area, 1, 1, NT);
    vv = reshape(Gi, [], 1);
    A = sparse(ii, jj, vv, n*NT, n*NT);
end

%% 逐单元装配
% function A = assembleStiffness(fem, node, elem, elem2dof)
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
%         grad_phi = fem.computeBasisGrad_all(tid, lam);
%         At = At + grad_phi' * (grad_phi * (w(i) * area));
%     end
% end