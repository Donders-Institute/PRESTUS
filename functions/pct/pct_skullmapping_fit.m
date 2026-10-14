function [m, b, info] = pct_skullmapping_fit(ute, labels, options)
% PCT_SKULLMAPPING_FIT  Subject-specific UTE->HU line from cortical / trabecular UTE histograms
%
% The fit behind skull_mapping 'kosciessa' (called by PCT_SKULLMAPPING with
% the defaults). Histograms of normalised UTE in charm cortical (7) and
% trabecular (8) bone give three points:
%   lower width crossing of the cortical histogram   -> anchors_hu(1)
%   crossover of the two histograms between peaks     -> anchors_hu(2)
%   upper width crossing of the trabecular histogram  -> anchors_hu(3)
% and a least-squares line through them, HU = m * UTE + b.
%
% Use as:
%   [m, b, info] = pct_skullmapping_fit(ute, labels)
%   [m, b, info] = pct_skullmapping_fit(ute, labels, anchors_hu=[2000 800 300], width_level=0.5)
%
% Input:
%   ute     - normalised UTE (soft-tissue peak = 1), any size
%   labels  - charm final_tissues labels, same size
%   options - name-value:
%     anchors_hu    - HU at the three points (default [1900 700 300])
%     width_level   - fraction of the histogram maximum at which the width
%                     crossings are taken (default 0.25: full width at quarter max)
%     interp_window - UTE range whose histogram bins are replaced by a spline
%                     to suppress the soft-tissue partial-volume bump
%                     (default [0.90 1.10]; [] = no replacement)
%     bin_edges     - histogram bin centres (default 0:0.01:1.5)
%
% Output:
%   m, b - slope and intercept (unrounded; PCT_SKULLMAPPING writes them rounded)
%   info - struct: x, y_cortical, y_trabecular (max-normalised histograms),
%          x_peak_cortical, x_peak_trabecular, x_points (the three points), anchors_hu
%
% See also: PCT_SKULLMAPPING, PCT_CREATE_PSEUDOCT

arguments
    ute
    labels
    options.anchors_hu (1,3) double = [1900 700 300]
    options.width_level (1,1) double {mustBeInRange(options.width_level, 0, 1, 'exclusive')} = 0.25
    options.interp_window double = [0.90 1.10]
    options.bin_edges (1,:) double = 0:0.01:1.5
end

x  = options.bin_edges;
y1 = hist(ute(labels == 7), x); % cortical
y2 = hist(ute(labels == 8), x); % trabecular

% Replace bins in the interpolation window with a spline (partial volume with soft tissue)
if ~isempty(options.interp_window)
    idx = find(x >= options.interp_window(1) & x <= options.interp_window(end));
    keep = setdiff(1:numel(x), idx);
    y1(idx) = interp1(x(keep), y1(keep), x(idx), 'spline');
    y2(idx) = interp1(x(keep), y2(keep), x(idx), 'spline');
end
y1 = y1 ./ max(y1);
y2 = y2 ./ max(y2);

% Peaks (trabecular peak weighted by x)
[pks1, locs1] = findpeaks(y1, x);
[pks2, locs2] = findpeaks(x .* y2, x);
[~, i1] = max(pks1); [~, i2] = max(pks2);
xMax1 = locs1(i1); xMax2 = locs2(i2);

% Lower width crossing of the cortical histogram, around the cortical peak
xLeft1 = width_crossing(x, y1, xMax1, options.width_level, 'left');

% Crossover between the peaks: minimum |moving average (2 bins) of y1 - y2|
range = x >= xMax1 & x <= xMax2;
d = y1(range) - y2(range);
w = 2;
mov = arrayfun(@(i) mean(d(i:i+w-1)), 1:numel(d) - w + 1);
if isempty(mov)
    warn('No crossover point between the cortical and trabecular peaks; using their midpoint. Visual inspection recommended.');
    xCross = mean([xMax1, xMax2]);
else
    [~, iMin] = min(abs(mov));
    xCross = x(iMin + floor(w / 2) + find(range, 1, 'first') - 1);
end

% Upper width crossing of the trabecular histogram, around the trabecular peak
% found above (x-weighted). The global maximum of y2 can be a near-zero mode
% (charm trabecular voxels with air-like UTE), which put this point below the
% crossover and gave a far too shallow line.
xRight2 = width_crossing(x, y2, xMax2, options.width_level, 'right');

X = [xLeft1, xCross, xRight2];
if ~issorted(X, 'strictascend')
    warn(sprintf(['Fit points not in order (cortical %.3f, crossover %.3f, trabecular %.3f); ' ...
        'the line is unreliable. Visual inspection recommended.'], X));
end
c = polyfit(X, options.anchors_hu, 1);
m = c(1); b = c(2);

info = struct('x', x, 'y_cortical', y1, 'y_trabecular', y2, 'x_peak_cortical', xMax1, ...
    'x_peak_trabecular', xMax2, 'x_points', X, 'anchors_hu', options.anchors_hu);
end

function xc = width_crossing(x, y, x_peak, level, side)
% x where y crosses level * y(x_peak) on the given side of the peak at x_peak (linear interpolation)
ip = find(x == x_peak, 1); lv = y(ip) * level;
if strcmp(side, 'left')
    i = find(y(1:ip) < lv, 1, 'last');
    xc = interp1(y(i:i+1), x(i:i+1), lv);
else
    i = find(y(ip:end) < lv, 1, 'first') + ip - 1;
    xc = interp1(y(i-1:i), x(i-1:i), lv);
end
end
