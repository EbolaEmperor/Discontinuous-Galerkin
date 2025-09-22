classdef Pk

properties
    ord
    locDof
    coef
    node
    Dlam
    area
    locLam
end

methods
    function obj = Pk(ord, node, elem)
        assert(ord >= 1);
        obj.ord = ord;
        obj.locDof = (ord+1) * (ord+2) / 2;
        obj.node = node;
        [obj.Dlam, obj.area] = gradbasis_my(node, elem);
        obj.locLam = zeros(obj.locDof, 3);
        ind = 0;
        for i = 0 : obj.ord
            for j = 0 : obj.ord - i
                ind = ind + 1;
                obj.locLam(ind, :) = [1-(i+j)/obj.ord, i/obj.ord, j/obj.ord];
            end
        end
        obj.coef = obj.initBasis();
    end

    function [elem2dof, nDof] = getDOF(obj, elem)
        NT = size(elem, 1);
        lagNodes = zeros(NT * obj.locDof, 2);
        for t = 1 : NT
            lagNodes((t-1)*obj.locDof + (1:obj.locDof),:) = obj.locLam * obj.node(elem(t,:),:);
        end
        tol = 1e-15;
        % Use approximate unique and find, since the lagrange nodes may not
        % be very exact. (floating error!)
        [lagNodes, ~, ~] = uniquetol(lagNodes, tol, 'ByRows', true);
        nDof = size(lagNodes, 1);
        elem2dof = zeros(NT, obj.locDof);
        for t = 1:NT
            locLagNodes = obj.locLam * obj.node(elem(t,:),:);
            [tf, idx] = ismembertol(locLagNodes, lagNodes, tol, 'ByRows', true, 'DataScale', 1);
            assert(all(tf));
            elem2dof(t, :) = idx(:).';
        end
        fprintf('DoF alligned: NT=%d , nDof=%d\n', NT, nDof);
    end

    function coef = initBasis(obj)
        A = zeros(obj.locDof);
        for i = 1 : obj.locDof
            A(i,:) = obj.getSpan(obj.locLam(i,:));
        end
        coef = A \ eye(obj.locDof);
    end

    function span = getSpan(obj, lam)
        span = polyBasisHomo3D(obj.ord, lam);
    end
    
    function span = gradSpan(obj, lam, dlam)
        span = dlam * polyBasisHomoGrad3D(obj.ord, lam);
    end

    function span = lapSpan(obj, lam, dlam)
        H = polyBasisHomoHess3D(obj.ord, lam);
        S = dlam.' * dlam;
        span = S(1,1).*H(1,:) + S(2,2).*H(2,:) + S(3,3).*H(3,:) ...
             + 2*S(1,2).*H(4,:) + 2*S(1,3).*H(5,:) + 2*S(2,3).*H(6,:);
    end

    function d2 = diff2Span(obj, lam, dlam, dir)
        v = dir(:);
        H = polyBasisHomoHess3D(obj.ord, lam);
        w = dlam.' * v;
        w1 = w(1); w2 = w(2); w3 = w(3);
        d2 = (w1*w1).*H(1,:) + (w2*w2).*H(2,:) + (w3*w3).*H(3,:) ...
           + 2*(w1*w2).*H(4,:) + 2*(w1*w3).*H(5,:) + 2*(w2*w3).*H(6,:);
    end

    function span = hessSpan(obj, lam, dlam)
        H = polyBasisHomoHess3D(obj.ord, lam);
        J = dlam;
        Nspan = size(H,2);
        span = zeros(3, Nspan);
        for k = 1:Nspan
            Hlam = [ H(1,k), H(4,k), H(5,k);
                     H(4,k), H(2,k), H(6,k);
                     H(5,k), H(6,k), H(3,k) ];
            Hxy  = J * Hlam * J.';
            span(1,k) = Hxy(1,1);
            span(2,k) = Hxy(2,2);
            span(3,k) = Hxy(1,2);
        end
    end

    function val = computeBasisValue_all(obj, tid, lam)
        span = obj.getSpan(lam);
        val = span * obj.coef;
    end
    
    function val = computeBasisGrad_all(obj, tid, lam)
        span = obj.gradSpan(lam, obj.Dlam(:,:,tid));
        val = span * obj.coef;
    end

    function span = computeBasisDlam_all(obj, lam)
        % 返回一个 3*n 的结果，第 i 行是 Phi_1,...,Phi_N 对 lam(i) 的偏导
        % 每个三角形算出来都是一样的
        span = polyBasisHomoGrad3D(obj.ord, lam) * obj.coef;
    end

    function val = computeBasisLaplace_all(obj, tid, lam)
        span = obj.lapSpan(lam, obj.Dlam(:,:,tid));
        val = span * obj.coef;
    end

    function val = computeBasisDirectedDiff2_all(obj, tid, lam, dir)
        span = obj.diff2Span(lam, obj.Dlam(:,:,tid), dir);
        val = span * obj.coef;
    end

    function val = computeBasisHessian_all(obj, tid, lam)
        span = obj.hessSpan(lam, obj.Dlam(:,:,tid));
        val = span * obj.coef;
    end

    function [quadL, w] = quad1d(obj)
        [quadL, w] = quadpts1(2*(obj.ord+1));
    end

    function [quadL, w] = quad2d(obj)
        [quadL, w] = quadpts2_my(2*(obj.ord+1));
    end

end

end