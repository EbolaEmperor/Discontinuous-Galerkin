classdef square
methods

    function [node, elem] = getMesh(obj, h)
        [node,elem] = squaremesh([0 1 0 1], h);
    end

    function isBdEdge = findBdry(obj, node, edge)
        NE = size(edge, 1);
        isHori = @(i) (node(edge(i,1), 2) == node(edge(i,2), 2));
        isVert = @(i) (node(edge(i,1), 1) == node(edge(i,2), 1));
        nodeInBond = @(i) (node(i,1)==0 || node(i,1) == 1) || (node(i,2)==0 || node(i,2)==1);
        isBdEdge = zeros(NE, 1);
        for i = 1 : NE
            if (isHori(i) || isVert(i)) && nodeInBond(edge(i,1)) && nodeInBond(edge(i,2))
                isBdEdge(i) = 1;
            end
        end
    end
    
    function isBdNode = findBdryNodes(obj, node)
        N = size(node, 1);
        tol = 4 * eps;
        nodeInBond = @(node) abs(node(1))<tol || abs(node(1)-1)<tol || abs(node(2))<tol || abs(node(2)-1)<tol;
        isBdNode = zeros(N, 1);
        for i = 1 : N
            if nodeInBond(node(i,:))
                isBdNode(i) = 1;
            end
        end
    end
    
end
end