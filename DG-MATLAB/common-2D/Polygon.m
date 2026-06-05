classdef Polygon
    properties
        vertices
    end
    
    methods
        function obj = Polygon(arg)
            if isscalar(arg)
                % Case 1: arg is the number of vertices for a regular polygon
                n = arg;
                theta = (0 : n - 1) / n * 2 * pi;
                obj.vertices = [cos(theta)', sin(theta)'];
            else
                % Case 2: arg is a matrix of vertices
                vertices = arg;
                if size(vertices, 2) ~= 2
                    error('Vertices must be an N x 2 matrix.');
                end
                
                % Sort vertices in counter-clockwise order
                center = mean(vertices, 1);
                vectors = vertices - center;
                angles = atan2(vectors(:,2), vectors(:,1));
                [~, idx] = sort(angles);
                obj.vertices = vertices(idx, :);
            end
        end

        function [node, elem] = getMesh(obj, h)
            if nargin < 2, h = 1; end
            % Use the centroid as the interior point
            center = mean(obj.vertices, 1);
            node = [obj.vertices; center];
            elem = delaunay(node(:,1), node(:,2));
            elem = fixorder(node,elem);
            k = ceil(log2(1 / h));
            for i = 1 : k
                [node,elem] = uniformrefine(node, elem);
            end
        end

        function isBdEdge = findBdry(obj, node, edge, bdc)
            % 对任意多边形：标记哪些边属于多边形边界的某一条边段
            % bdc 可选：为 M×1 或 1xM 的逻辑向量，用来选择多边形的边（M 为多边形边数）
            
            poly = obj.vertices;
            M = size(poly,1);
            
            if nargin < 4 || isempty(bdc)
                bdc = true(1, M);
            else
                bdc = logical(reshape(bdc, 1, []));
                if numel(bdc) < M
                    % 如果长度不匹配，按可用长度循环或补全为 true
                    bdc = repmat(bdc, 1, ceil(M/numel(bdc)));
                    bdc = bdc(1:M);
                end
            end
            
            tol = 1e-15;
            NE = size(edge,1);
            isBdEdge = false(NE,1);
    
            % 辅助：点是否在线段上
            pointOnSeg = @(p,a,b) (abs(det([b-a; p-a])) <= tol*max(1,norm(b-a))) & (dot(p-a,p-b) <= tol);
    
            % 对每条网格边，判断其两个端点是否同时位于某个被选中的多边形边段上
            for ei = 1:NE
                i1 = edge(ei,1); i2 = edge(ei,2);
                p1 = node(i1,:); p2 = node(i2,:);
                onSome = false;
                for j = 1:M
                    if ~bdc(j), continue; end
                    a = poly(j,:);
                    b = poly(mod(j, M) + 1, :);
                    if pointOnSeg(p1, a, b) && pointOnSeg(p2, a, b)
                        onSome = true; break;
                    end
                end
                isBdEdge(ei) = onSome;
            end
        end
        
        function isBdNode = findBdryNodes(obj, node)
            % 对任意多边形：判断每个点是否在多边形的某一条边段上
            pts = node;
            tol = 1e-15;
            poly = obj.vertices;
            M = size(poly,1);
            N = size(pts,1);
            isBdNode = false(N,1);
    
            for j = 1:M
                a = poly(j,:);
                b = poly(mod(j, M) + 1, :);
                
                v = b - a;
                len_v = norm(v);
                
                % 向量化计算: ap = p - a, bp = p - b
                ap = pts - a;
                bp = pts - b;
                
                % 叉积 (2D) 判断共线: v x ap
                crossVal = v(1) * ap(:,2) - v(2) * ap(:,1);
                
                % 点积判断在线段范围内: ap . bp <= 0
                dotVal = sum(ap .* bp, 2);
                
                onLine = abs(crossVal) <= tol * max(1, len_v);
                onSeg = dotVal <= tol;
                
                isBdNode = isBdNode | (onLine & onSeg);
            end
        end

        function isVertex = findVertex(obj, p)
            isVertex = false(size(p,1), 1);
            for i = 1 : size(obj.vertices, 1)
                err = vecnorm(p - obj.vertices(i,:), inf, 2);
                isVertex = isVertex | (err < 1e-15);
            end
        end
        
    end
end