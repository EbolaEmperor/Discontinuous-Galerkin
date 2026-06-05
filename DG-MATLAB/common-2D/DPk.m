classdef DPk

properties
    ord
    locDof
    coef
    Dlam
    area
end

methods
    function obj = DPk(ord, node, elem)
        obj.ord = ord;
        obj.locDof = (ord+1) * (ord+2) / 2;
        obj.coef = obj.initBasis();
        [obj.Dlam, obj.area] = gradbasis_my(node, elem);
    end

    function [elem2dof, nDof] = getDOF(obj, elem)
        NT = size(elem, 1);
        elem2dof = 1 : NT * obj.locDof;
        elem2dof = reshape(elem2dof, [obj.locDof, NT])';
        nDof = obj.locDof * NT;
        fprintf('DoF alligned: NT=%d, nDof=%d\n', NT, nDof);
    end

    function coef = initBasis(obj)
        A = zeros(obj.locDof);
        if obj.ord == 0
            pnt = ones(1, 3) / 3;
        else
            pnt = zeros(obj.locDof, 3);
            ind = 0;
            for i = 0 : obj.ord
                for j = 0 : obj.ord - i
                    ind = ind + 1;
                    pnt(ind, :) = [1-(i+j)/obj.ord, i/obj.ord, j/obj.ord];
                end
            end
        end
        for i = 1 : obj.locDof
            A(i,:) = obj.getSpan(pnt(i,:));
        end
        coef = A \ eye(obj.locDof);
    end

    function span = getSpan(obj, lam)
        span = polyBasisHomo3D(obj.ord, lam);
    end
    
    function span = gradSpan(obj, lam, dlam)
        span = dlam * polyBasisHomoGrad3D(obj.ord, lam);
    end

    function coef = getBasis(obj, tid)
        coef = obj.coef;
    end

    function val = computeBasisValue(obj, a, pts)
        N = size(pts,1);
        val = zeros(N,1);
        for i = 1 : N
            span = obj.getSpan(pts(i,:));
            val(i) = span * a(:);
        end
    end
    
    function val = computeBasisGrad(obj, a, pts)
        N = size(pts,1);
        val = zeros(N,2);
        for i = 1 : N
            span = obj.gradSpan(pts(i,:), obj.Dlam(:,:,tid));
            val(i) = span * a(:);
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

    function [quadL, w] = quad1d(obj)
        [quadL, w] = quadpts1(2*(obj.ord+1));
    end

    function [quadL, w] = quad2d(obj)
        [quadL, w] = quadpts2_my(2*(obj.ord+1));
    end

end

end