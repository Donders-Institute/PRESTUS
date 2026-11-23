function kind = uncertainty_fov_metric_kind(name)
% UNCERTAINTY_FOV_METRIC_KIND  How a metric is rescaled from the partial FOV to the whole grid
%
% Use as:
%   kind = uncertainty_fov_metric_kind(name)
%
% Input:
%   name - metrics-table column name
%
% Output:
%   kind - 'ratio' (intensity, pressure, MI, temperature rise, CEM43: multiplied
%          by the factor) | 'rise' (absolute temperature: its rise above the
%          baseline is multiplied) | 'none' (not rescaled)
%
% See also: UNCERTAINTY_FOV_MISMATCH, UNCERTAINTY_FOV_RESCALE

if ~isempty(regexp(name, '^(Isppa|Ispta|Ipa|MI|Psptp|Ptp|riseT|CEM43)($|_)', 'once'))
    kind = 'ratio';
elseif ~isempty(regexp(name, '^(maxT|endT)($|_)', 'once'))
    kind = 'rise';
else
    kind = 'none';
end
end
