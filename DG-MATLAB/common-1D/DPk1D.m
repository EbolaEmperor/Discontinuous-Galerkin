classdef DPk1D

properties
    ord
    locDof
end

methods
    function obj = DPk1D(k)
        obj.ord = k;
        obj.locDof = k + 1;
    end

    function n = nDof(obj, NT)
        n = obj.locDof * NT;
    end

    function idx = dofMap(obj, t)
        idx = (t-1) * obj.locDof + (1:obj.locDof);
    end

    function span = getSpan(obj, lam)
        span = (lam(:,1) .^ (obj.ord:-1:0)) .* (lam(:,2) .^ (0:obj.ord));
    end

    function span = diffSpan(obj, lam, h)
        span = (-1/h) * (obj.ord:-1:0) .* (lam(:,1) .^ [obj.ord-1:-1:0, 0]) .* (lam(:,2) .^ (0:obj.ord)) + ...
               (1/h) * (0:obj.ord) .* (lam(:,2) .^ [0, 0:obj.ord-1]) .* (lam(:,1) .^ (obj.ord:-1:0));
    end
end

end