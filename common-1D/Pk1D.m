classdef Pk1D

properties
    ord
    locDof
    coef
end

methods
    function obj = Pk1D(k)
        assert(k >= 1);
        obj.ord = k;
        obj.locDof = k + 1;
        A = zeros(k+1, k+1);
        for i = 1 : k+1
            s = (i-1) / k;
            A(i,:) = ((1-s) .^ (obj.ord:-1:0)) .* (s .^ (0:obj.ord));
        end
        obj.coef = inv(A);
    end

    function span = getSpan(obj, lam)
        span = (lam(:,1) .^ (obj.ord:-1:0)) .* (lam(:,2) .^ (0:obj.ord));
        % make sure that span_i(eta_j) = delta_{ij} for all lagrange nodes eta_j
        span = span * obj.coef;
    end

    function span = diffSpan(obj, lam, h)
        span = (-1/h) * (obj.ord:-1:0) .* (lam(:,1) .^ [obj.ord-1:-1:0, 0]) .* (lam(:,2) .^ (0:obj.ord)) + ...
               (1/h) * (0:obj.ord) .* (lam(:,2) .^ [0, 0:obj.ord-1]) .* (lam(:,1) .^ (obj.ord:-1:0));
        span = span * obj.coef;
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
    
        span = (term1 + term2 + term3) / (h^2) * obj.coef;
    end
end

end