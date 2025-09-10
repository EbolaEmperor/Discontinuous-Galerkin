function [errH1, errL2, errSemiH1] = getH1Err1D(fem, grid, uh, u, du)
    NT   = length(grid) - 1;
    nDof = fem.locDof;

    [quadL, w] = quadpts1(2 * nDof);
    w  = w(:);
    nq = numel(w);

    phi  = zeros(nDof, nq);
    dphi = zeros(nDof, nq);
    for i = 1:nq
        phi(:, i)  = fem.getSpan(quadL(i, :)).';
        dphi(:, i) = fem.diffSpan(quadL(i, :), 1).';
    end

    v1v2 = [grid(1:NT); grid(2:NT+1)];
    h    = diff(grid).';
    pts = quadL * v1v2;

    try
        Uex = u(pts).';
    catch
        Uex = arrayfun(u, pts).';
    end
    
    try
        DUex = du(pts).';
    catch
        DUex = arrayfun(du, pts).';
    end

    uh_mat = reshape(uh, [nDof, NT]);

    Uh_val = (uh_mat.' * phi);
    Uh_dx  = (uh_mat.' * dphi);
    Uh_dx  = Uh_dx ./ h;

    diff0   = Uh_val - Uex;
    diff1   = Uh_dx  - DUex;
    errL2   = sum( (diff0.^2) * w .* h );
    errSemiH1 = sum( (diff1.^2) * w .* h );

    errH1     = sqrt(errL2 + errSemiH1);
    errL2     = sqrt(errL2);
    errSemiH1 = sqrt(errSemiH1);
end

%% 以下为原始的手写代码，向量化代码由 ChatGPT 根据原始代码自动生成
% function [errH1, errL2, errSemiH1] = getH1Err1D(fem, grid, uh, u, du)
% 
%     NT = length(grid) - 1;
%     nDof = fem.locDof;
%     errL2 = 0;
%     errSemiH1 = 0;
% 
%     [quadL, w] = quadpts1(2 * nDof);
%     w = w(:)';
%     nq = numel(w);
%     phi = zeros(nDof, nq);
%     dphi = zeros(nDof, nq);
%     for i = 1 : nq
%         lam = quadL(i, :);
%         phi(:,i) = fem.getSpan(lam)';
%         dphi(:,i) = fem.diffSpan(lam, 1)';
%     end
% 
%     for t = 1 : NT
%         vtx = grid(t:t+1);
%         h = vtx(2) - vtx(1);
%         uh_t = uh((t-1)*nDof + (1:nDof))';
%         u_ex = zeros(1, nq);
%         du_ex = zeros(1, nq);
%         for i = 1 : nq
%             lam = quadL(i, :);
%             pts = dot(lam, vtx);
%             u_ex(i) = u(pts);
%             du_ex(i) = du(pts);
%         end
%         errL2 = errL2 + sum( (uh_t * phi - u_ex).^2 .* w ) * h;
%         errSemiH1 = errSemiH1 + sum( (uh_t * dphi/h - du_ex).^2 .* w ) * h;
%     end
% 
%     errH1 = sqrt(errL2 + errSemiH1);
%     errL2 = sqrt(errL2);
%     errSemiH1 = sqrt(errSemiH1);
% 
% end