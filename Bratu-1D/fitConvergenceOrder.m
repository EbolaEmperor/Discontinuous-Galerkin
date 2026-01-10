function [p, q] = fitConvergenceOrder(e)
%FITCONVERGENCEORDER Estimate convergence order from an error sequence.

    e = e(:);
    n = numel(e);
    if n < 4 || any(~isfinite(e))
        p = NaN;
        q = NaN;
        return;
    end

    tailLen = min(200, n);
    idxTail = (n - tailLen + 1) : n;
    et = e(idxTail);
    it = (idxTail - 1).';
    good = isfinite(et) & (et > 0);
    et = et(good);
    it = it(good);
    if numel(et) < 10
        p = NaN;
        q = NaN;
        return;
    end

    s = polyfit(it, log(et), 1);
    q = exp(s(1));

    e1 = e(1:end-2);
    e2 = e(2:end-1);
    e3 = e(3:end);
    pk = (log(e3) - log(e2)) ./ (log(e2) - log(e1));
    pk = pk(isfinite(pk) & ~isnan(pk));
    if isempty(pk)
        p = NaN;
    else
        p = median(pk(max(1, numel(pk)-50):end));
    end
end

