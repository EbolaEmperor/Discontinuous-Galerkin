classdef Pk1D

properties
    ord
    locDof
end

methods
    function obj = Pk1D(k)
        assert(k >= 1);
        obj.ord = k;
        obj.locDof = k + 1;
    end

    function n = nDof(obj, NT)
        n = NT * obj.ord + 1;
    end

    function idx = dofMap(obj, t)
        idx = (t-1) * (obj.locDof-1) + (1:obj.locDof);
    end

    function span = getSpan(obj, lam)
        span = (lam(:,1) .^ (obj.ord:-1:0)) .* (lam(:,2) .^ (0:obj.ord));
    end

    function span = diffSpan(obj, lam, h)
        span = (-1/h) * (obj.ord:-1:0) .* (lam(:,1) .^ [obj.ord-1:-1:0, 0]) .* (lam(:,2) .^ (0:obj.ord)) + ...
               (1/h) * (0:obj.ord) .* (lam(:,2) .^ [0, 0:obj.ord-1]) .* (lam(:,1) .^ (obj.ord:-1:0));
    end

    function span = diff2Span(obj, lam, h)
        nq = size(lam,1);
        if obj.ord <= 1
            span = zeros(nq, obj.locDof);
            return;
        end
    
        a = (obj.ord:-1:0);
        b = (0:obj.ord);
    
        pow_l1_m2 = [obj.ord-2:-1:0, 0, 0];
        pow_l1_m1 = [obj.ord-1:-1:0, 0];
        pow_l2_m1 = [0, 0:obj.ord-1];
        pow_l2_m2 = [0, 0, 0:obj.ord-2];

        term1 = (a .* (a - 1)) .* (lam(:,1) .^ pow_l1_m2) .* (lam(:,2) .^ b);
        term2 = (-2) * (a .* b) .* (lam(:,1) .^ pow_l1_m1) .* (lam(:,2) .^ pow_l2_m1);
        term3 = (b .* (b - 1)) .* (lam(:,1) .^ a) .* (lam(:,2) .^ pow_l2_m2);
    
        span = (term1 + term2 + term3) / (h^2);
    end
end

end