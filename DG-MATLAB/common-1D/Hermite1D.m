classdef Hermite1D

properties
    ord
    locDof
    coef
end

methods
    function obj = Hermite1D(k, h)
        assert(k >= 1);
        obj.ord = k;
        obj.locDof = 2 * (k + 1);
        A = zeros(obj.locDof);
        ff = 2 * k + 1;
        for i = 1 : k+1
            s = (i-1) / k;
            A(2*i-1, :) = ((1-s) .^ (ff:-1:0)) .* (s .^ (0:ff));
            A(2*i, 1:end-1) = -(1/h) * (ff:-1:1) .* ((1-s) .^ (ff-1:-1:0)) .* (s .^ (0:ff-1));
            A(2*i, 2:end) = A(k+1+i,2:end) + (1/h) * (1:ff) .* ((1-s) .^ (ff-1:-1:0)) .* (s .^ (0:ff-1));
        end
        obj.coef = inv(A);
    end

    function n = nDof(obj, NT)
        n = NT * (obj.locDof-2) + 2;
    end

    function idx = dofMap(obj, t)
        idx = (t-1) * (obj.locDof-2) + (1:obj.locDof);
    end

    function span = getSpan(obj, lam)
        ff = 2 * obj.ord + 1;
        span = (lam(:,1) .^ (ff:-1:0)) .* (lam(:,2) .^ (0:ff));
        % make sure that span_i(eta_j) = delta_{ij} for all lagrange nodes eta_j
        span = span * obj.coef;
    end

    function span = diffSpan(obj, lam, h)
        span = zeros(size(lam,1), obj.locDof);
        ff = 2 * obj.ord + 1;
        span(:, 1:ff) = (-1/h) * (ff:-1:1) .* (lam(:,1) .^ (ff-1:-1:0)) .* (lam(:,2) .^ (0:ff-1));
        span(:, 2:ff+1) = span(:,2:ff+1) + (1/h) * (1:ff) .* (lam(:,1) .^ (ff-1:-1:0)) .* (lam(:,2) .^ (0:ff-1));
        span = span * obj.coef;
    end

    function span = diff2Span(obj, lam, h)
        nq = size(lam,1);
        if obj.ord <= 1
            span = zeros(nq, obj.locDof);
            return;
        end
    
        ff = 2 * obj.ord + 1;
        a = (ff:-1:0);
        b = (0:ff);
    
        pow_l1_m2 = [ff-2:-1:0, 0, 0];
        pow_l1_m1 = [ff-1:-1:0, 0];
        pow_l2_m1 = [0, 0:ff-1];
        pow_l2_m2 = [0, 0, 0:ff-2];

        term1 = (a .* (a - 1)) .* (lam(:,1) .^ pow_l1_m2) .* (lam(:,2) .^ b);
        term2 = (-2) * (a .* b) .* (lam(:,1) .^ pow_l1_m1) .* (lam(:,2) .^ pow_l2_m1);
        term3 = (b .* (b - 1)) .* (lam(:,1) .^ a) .* (lam(:,2) .^ pow_l2_m2);
    
        span = (term1 + term2 + term3) / (h^2) * obj.coef;
    end
end

end