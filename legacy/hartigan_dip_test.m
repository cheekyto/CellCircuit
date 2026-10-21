function [dip, p_value, boot_dips] = hartigan_dip_test(x, n_boot)
% HARTIGAN_DIP_TEST  Hartigan's Dip Test for unimodality.
%
%   [dip, p_value] = hartigan_dip_test(x)
%   [dip, p_value] = hartigan_dip_test(x, n_boot)
%
%   x       : numeric vector of data
%   n_boot  : number of bootstrap replicates for the p-value
%             (default 2000)
%
%   dip     : Hartigan's dip statistic. Larger values indicate stronger
%             departure from unimodality.
%   p_value : bootstrap p-value against the null of a unimodal
%             (uniform) distribution. Small p-value (e.g. < 0.05)
%             is evidence AGAINST unimodality, i.e. evidence of
%             multimodality (typically bimodality).
%
%   Method: the dip statistic is the smallest possible maximum
%   deviation between the empirical CDF and any unimodal CDF.
%   The optimal unimodal fit consists of the greatest convex minorant
%   (GCM) of the ECDF up to the mode and the least concave majorant
%   (LCM) after the mode. This implementation searches over candidate
%   mode locations and computes the GCM/LCM via convex hulls.
%
%   Note: this is a direct-from-theory implementation, not a port of
%   the original AS 217 Fortran/C code, so exact numerical values may
%   differ slightly from R's diptest package. For high-stakes /
%   publication use, cross-check against a reference implementation.
%
%   Complexity: O(n^2) in the number of data points. Fine for datasets
%   up to a few thousand points; for larger datasets, consider
%   subsampling before testing.

    if nargin < 2 || isempty(n_boot)
        n_boot = 2000;
    end

    x = sort(x(:));
    x = x(~isnan(x));
    n = numel(x);

    if n < 4
        dip = 0;
        p_value = 1;
        boot_dips = [];
        warning('Need at least 4 data points for the dip test.');
        return;
    end

    ecdf = (1:n)' / n;
    dip = local_compute_dip(x, ecdf);

    % Bootstrap null distribution: dip stat is affine-invariant in x,
    % so simulating from Uniform(0,1) is equivalent to using the
    % data's own range.
    boot_dips = zeros(n_boot, 1);
    for b = 1:n_boot
        xs = sort(rand(n, 1));
        boot_dips(b) = local_compute_dip(xs, ecdf);
    end

    p_value = mean(boot_dips >= dip);
end

% ---------------------------------------------------------------
function d = local_compute_dip(x, ecdf)
    n = numel(x);
    best = inf;

    for m = 2:(n-1)
        idxL = 1:m;
        gcm = local_lower_hull_interp(x(idxL), ecdf(idxL));
        devL = max(abs(gcm - ecdf(idxL)));

        idxR = m:n;
        lcm = local_upper_hull_interp(x(idxR), ecdf(idxR));
        devR = max(abs(lcm - ecdf(idxR)));

        dev = max(devL, devR);
        if dev < best
            best = dev;
        end
    end

    d = best;
end

% ---------------------------------------------------------------
function yhat = local_lower_hull_interp(x, y)
% Greatest convex minorant of points (x,y), x sorted ascending.
% Returns the hull's y-value at each input x (piecewise linear).
    n = numel(x);
    hull_idx = zeros(n, 1);
    k = 0;
    for i = 1:n
        while k >= 2 && local_cross(x(hull_idx(k-1)), y(hull_idx(k-1)), ...
                                     x(hull_idx(k)),   y(hull_idx(k)), ...
                                     x(i),             y(i)) <= 0
            k = k - 1;
        end
        k = k + 1;
        hull_idx(k) = i;
    end
    hull_idx = hull_idx(1:k);
    yhat = interp1(x(hull_idx), y(hull_idx), x, 'linear');
end

% ---------------------------------------------------------------
function yhat = local_upper_hull_interp(x, y)
% Least concave majorant = negative of lower hull of negated y.
    yhat = -local_lower_hull_interp(x, -y);
end

% ---------------------------------------------------------------
function c = local_cross(x1, y1, x2, y2, x3, y3)
    c = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
end
