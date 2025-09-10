classdef DPk1D

properties
    ord
    locDof
    coef
end

methods
    function obj = DPk1D(k)
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
end

end