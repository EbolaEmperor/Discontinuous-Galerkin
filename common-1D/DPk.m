classdef DPk

properties
    ord
    locDof
end

methods
    function obj = DPk(k)
        obj.ord = k;
        obj.locDof = k + 1;
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