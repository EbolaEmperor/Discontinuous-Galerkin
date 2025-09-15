classdef vecDPk

properties
    ord
    pW
    locDof
    coef
    Dlam
end

methods
    function obj = vecDPk(ord, node, elem)
        obj.ord = ord;
        obj.pW = ord;
        obj.locDof = (ord+1) * (ord+2) / 2;
        obj.coef = obj.initBasis();
        obj.Dlam = gradbasis_my(node, elem);
    end

    function [elem2dof, nDof] = getDOF(obj, elem)
        NT = size(elem, 1);
        elem2dof = zeros(NT, obj.locDof);
        for t = 1:NT
            elem2dof(t,:) = (t-1)*(obj.locDof) + (1:obj.locDof);
        end
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
        z = zeros(1, length(s));
        span = [span, z; z, span];
    end
    
    function span = divSpan(obj, lam, dlam)
        span = dlam * polyBasisHomoGrad3D(obj.ord, lam);
        span = [span(1,:), span(2,:)];
    end

    function val = computeBasisValue_all(obj, tid, lam)
        span = obj.getSpan(lam);
        val = span * obj.coef;
    end
    
    function val = computeBasisDiv_all(obj, tid, lam)
        span = obj.divSpan(lam, obj.Dlam(:,:,tid));
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