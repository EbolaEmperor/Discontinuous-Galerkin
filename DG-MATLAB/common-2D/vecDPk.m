classdef vecDPk

properties
    ord
    locDof
    Dlam
    area
end

methods
    function obj = vecDPk(ord, node, elem)
        obj.ord = ord;
        obj.locDof = (ord+1) * (ord+2);
        [obj.Dlam, obj.area] = gradbasis_my(node, elem);
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

    function span = getSpan(obj, lam)
        span = polyBasisHomo3D(obj.ord, lam);
        z = zeros(1, length(span));
        span = [span, z; z, span];
    end
    
    function span = divSpan(obj, lam, dlam)
        span = dlam * polyBasisHomoGrad3D(obj.ord, lam);
        span = [span(1,:), span(2,:)];
    end

    function val = computeBasisValue_all(obj, tid, lam)
        val = obj.getSpan(lam);
    end
    
    function val = computeBasisDiv_all(obj, tid, lam)
        val = obj.divSpan(lam, obj.Dlam(:,:,tid));
    end

    function [quadL, w] = quad1d(obj)
        [quadL, w] = quadpts1(2*(obj.ord+1));
    end

    function [quadL, w] = quad2d(obj)
        [quadL, w] = quadpts2_my(2*(obj.ord+1));
    end

end

end