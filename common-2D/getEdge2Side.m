% edge2side 返回边左右两侧的三角形编号，对于边 [u,v]，定向为 u -> v，
% 此时 edge2side 里的两个三角形编号，第一个是左边的，第二个是右边的。
% 如果有一侧是区域外，则标记为 0

function [edge, edge2side] = getEdge2Side(node, elem)
    edge = [elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])];
    edge = sort(edge, 2);
    [edge,~,~] = unique(edge,'rows');
    NE = size(edge, 1);
    NT = size(elem, 1);

    edge2side = zeros(NE, 2);
    for t = 1 : NT
        for loc = [2, 3, 1]
            nxt = mod(loc, 3) + 1;
            u = elem(t, loc);
            v = elem(t, nxt);
            w = elem(t, 6-loc-nxt);
            if u > v
                tmp = u; u = v; v = tmp;
            end
            ind = find(all(edge == [u, v], 2), 1);
            if cross2d(node(v,:) - node(u,:), node(w,:) - node(u,:)) > 0
                edge2side(ind, 1) = t;
            else
                edge2side(ind, 2) = t;
            end
        end
    end
end