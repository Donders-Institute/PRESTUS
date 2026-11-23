function mismatch = uncertainty_fov_mismatch(t_full, t_fov, temp_0)
% UNCERTAINTY_FOV_MISMATCH  Whole-grid vs. partial-FOV factors of the default medium
%
% With simulation.uncertainty_fov = 'partial', the default medium is simulated
% twice: on the whole grid and in the partial acoustic FOV the bound variants
% use. This function compares the two metrics tables and returns, per metric,
% the factor that maps a partial-FOV value onto the whole grid
% (UNCERTAINTY_FOV_RESCALE applies it to the bound variants).
%
% Metric kinds, by column name:
%   'ratio'  intensity, pressure, MI, temperature rise, CEM43
%            (Isppa*, Ispta*, Ipa*, MI*, Psptp*, Ptp*, riseT*, CEM43*):
%            factor = whole / partial
%   'rise'   absolute temperatures (maxT*, endT*): the factor of the rise above
%            the tissue's baseline, factor = (whole - T0) / (partial - T0)
%   'none'   everything else (focal position and size, volumes, flags): not
%            rescaled; only the difference is reported
% CEM43 is exponential in temperature, so its ratio is a first-order correction.
%
% A factor is NaN (the metric is then left unscaled) when either value is
% missing or not finite, or when the partial-FOV value (or rise) is not positive.
%
% Use as:
%   mismatch = uncertainty_fov_mismatch(t_full, t_fov)
%   mismatch = uncertainty_fov_mismatch(t_full, t_fov, temp_0)
%
% Input:
%   t_full  - metrics table of the default medium on the whole grid
%   t_fov   - metrics table of the default medium in the partial FOV
%   temp_0  - baseline temperatures [deg C]: struct with one field per tissue
%             (thermal.temp_0), or a scalar; tissues it lacks, and metrics
%             without a tissue suffix, use 37 (default: 37)
%
% Output:
%   mismatch - table with one row per numeric metric of both tables:
%              metric, kind, baseline, full, fov, factor, diff_pct
%              (diff_pct = 100 * (fov - full) / full; last table row is used)
%
% See also: UNCERTAINTY_FOV_RESCALE, UNCERTAINTY_PIPELINE, GENERATE_UNCERTAINTY_REPORT

arguments
    t_full table
    t_fov  table
    temp_0 = 37
end

names = intersect(t_full.Properties.VariableNames, t_fov.Properties.VariableNames, 'stable');
keep  = cellfun(@(n) isnumeric(t_full.(n)) && isnumeric(t_fov.(n)) && ...
                     ~isempty(t_full.(n)) && ~isempty(t_fov.(n)), names);
names = names(keep);
n     = numel(names);

kind     = cell(n, 1);
baseline = NaN(n, 1);
full     = NaN(n, 1);
fov      = NaN(n, 1);
factor   = NaN(n, 1);
for m = 1:n
    kind{m} = uncertainty_fov_metric_kind(names{m});
    full(m) = double(t_full.(names{m})(end));
    fov(m)  = double(t_fov.(names{m})(end));
    switch kind{m}
        case 'ratio'
            if isfinite(full(m)) && isfinite(fov(m)) && fov(m) > 0 && full(m) >= 0
                factor(m) = full(m) / fov(m);
            end
        case 'rise'
            baseline(m) = uncertainty_fov_baseline(names{m}, temp_0);
            rise_full = full(m) - baseline(m);
            rise_fov  = fov(m)  - baseline(m);
            if isfinite(rise_full) && isfinite(rise_fov) && rise_fov > 0 && rise_full >= 0
                factor(m) = rise_full / rise_fov;
            end
    end
end
diff_pct = 100 * (fov - full) ./ full;
diff_pct(~isfinite(diff_pct)) = NaN;

mismatch = table(names(:), kind, baseline, full, fov, factor, diff_pct, ...
    'VariableNames', {'metric', 'kind', 'baseline', 'full', 'fov', 'factor', 'diff_pct'});
end
