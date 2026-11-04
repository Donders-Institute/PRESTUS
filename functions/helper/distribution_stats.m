function s = distribution_stats(x, edges)
% DISTRIBUTION_STATS  Summary statistics and histogram of a numeric vector
%
% No Statistics Toolbox dependency (percentiles via linear interpolation on
% sorted values). Without edges, the bin count scales with sample size (sqrt
% rule), capped to keep JSON caches small.
%
% Use as:
%   s = distribution_stats(x)
%   s = distribution_stats(x, edges)
%
% Input:
%   x     - numeric array; NaNs are ignored
%   edges - (optional) histogram bin edges, e.g. shared across distributions
%
% Output:
%   s - struct with .n, .mean, .median, .std, .min, .max, .p5, .p25, .p75,
%       .p95 and .histogram (.counts, .edges)
%
% See also: MEDIUM_SETUP, PCT_PREDICT_SKULL_PROPERTIES

    x = double(x(:));
    x = x(~isnan(x));
    s = struct();
    s.n = numel(x);
    if s.n == 0
        s.mean = NaN; s.median = NaN; s.std = NaN; s.min = NaN; s.max = NaN;
        s.p5 = NaN; s.p25 = NaN; s.p75 = NaN; s.p95 = NaN;
        s.histogram = struct('counts', [], 'edges', []);
        return
    end
    s.mean   = mean(x);
    s.median = median(x);
    s.std    = std(x);
    s.min    = min(x);
    s.max    = max(x);
    s.p5  = percentile(x, 5);
    s.p25 = percentile(x, 25);
    s.p75 = percentile(x, 75);
    s.p95 = percentile(x, 95);
    if nargin < 2 || isempty(edges)
        n_bins = min(50, max(10, round(sqrt(s.n))));
        [counts, edges] = histcounts(x, n_bins);
    else
        counts = histcounts(x, edges);
    end
    s.histogram = struct('counts', counts, 'edges', edges);
end

function p = percentile(x, q)
% Percentile via linear interpolation on sorted values.
    x = sort(x(:));
    n = numel(x);
    pos = (q / 100) * (n - 1) + 1;
    lo = min(max(floor(pos), 1), n);
    hi = min(max(ceil(pos), 1), n);
    frac = pos - lo;
    p = x(lo) + frac * (x(hi) - x(lo));
end
