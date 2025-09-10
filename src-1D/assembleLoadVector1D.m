function F = assembleLoadVector1D(fem, grid, f)
    NT   = length(grid) - 1;
    nDof = fem.locDof;

    [quadL, w] = quadpts1(2 * (nDof - 1));
    w  = w(:);
    nq = numel(w);

    phi = zeros(nDof, nq);
    for i = 1:nq
        phi(:, i) = fem.getSpan(quadL(i, :)).';
    end

    h = diff(grid);
    v1v2 = [grid(1:NT) ; grid(2:NT+1)];
    pts  = quadL * v1v2;

    try
        fq = f(pts);
    catch
        fq = arrayfun(f, pts);
    end

    weighted = bsxfun(@times, fq, w);
    Fmat = phi * weighted;
    Fmat = bsxfun(@times, Fmat, h);
    F = Fmat(:);
end

%% 以下为原始的手写代码，向量化代码由 ChatGPT 根据原始代码自动生成
% function F = assembleLoadVector1D(fem, grid, f)
%     NT = length(grid) - 1;
%     nDof = fem.locDof;
%     F = zeros(nDof, NT);
% 
%     [quadL, w] = quadpts1(2 * (nDof - 1));
%     w = w(:)';
%     nq = numel(w);
%     phi = zeros(nDof, nq);
%     for i = 1 : nq
%         lam = quadL(i, :);
%         phi(:,i) = fem.getSpan(lam)';
%     end
% 
%     for t = 1 : NT
%         vtx = grid(t:t+1);
%         fq = zeros(1, nq);
%         for i = 1 : nq
%             lam = quadL(i, :);
%             pts = dot(lam, vtx);
%             fq(i) = f(pts);
%         end
%         F(:,t) = sum(phi .* fq .* w, 2) * (vtx(2) - vtx(1));
%     end
% 
%     F = F(:);
% end