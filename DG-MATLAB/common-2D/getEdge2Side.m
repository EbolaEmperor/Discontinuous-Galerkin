%% 向量化处理：寻找每条边两侧单元的编号
function [edge, edge2side] = getEdge2Side(node, elem)
    edge = [elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])];
    edge = sort(edge, 2);
    edge = unique(edge, 'rows');
    NE   = size(edge, 1);
    NT   = size(elem, 1);

    edge2side = zeros(NE, 2);
    nNode   = size(node,1);
    edgeIdS = sparse(edge(:,1), edge(:,2), 1:NE, nNode, nNode);
    tIdx = (1:NT)';

    for loc = [2, 3, 1]
        nxt = mod(loc, 3) + 1;
        thr = 6 - loc - nxt;
        
        u0 = elem(:, loc);
        v0 = elem(:, nxt);
        w  = elem(:, thr);
        u = min(u0, v0);
        v = max(u0, v0);

        eIdx = full(edgeIdS(sub2ind([nNode, nNode], u, v)));  % NT x 1

        pU = node(u, :);
        pV = node(v, :);
        pW = node(w, :);
        a  = pV - pU;
        b  = pW - pU;
        crossz = a(:,1).*b(:,2) - a(:,2).*b(:,1);

        isLeft = crossz > 0;
        Lrows  = eIdx(isLeft);
        Rrows  = eIdx(~isLeft);

        edge2side(Lrows, 1) = tIdx(isLeft);
        edge2side(Rrows, 2) = tIdx(~isLeft);
    end
end

%% 逐单元处理
% function [edge, edge2side] = getEdge2Side(node, elem)
%     edge = [elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])];
%     edge = sort(edge, 2);
%     [edge,~,~] = unique(edge,'rows');
%     NE = size(edge, 1);
%     NT = size(elem, 1);
% 
%     edge2side = zeros(NE, 2);
%     for loc = [2, 3, 1]
%         nxt = mod(loc, 3) + 1;
%         pp = elem(:, [loc, nxt]);
%         w = elem(:, 6-loc-nxt);
%         pp = sort(pp, 2);
%         u = pp(:, 1);
%         v = pp(:, 2);
% 
%         % 这下面我不会写了
%         ind = find(all(edge == [u, v], 2), 1);
%         if cross2d(node(v,:) - node(u,:), node(w,:) - node(u,:)) > 0
%             edge2side(ind, 1) = t;
%         else
%             edge2side(ind, 2) = t;
%         end
% 
%     end
% end
