function t0 = uncertainty_fov_baseline(name, temp_0)
% UNCERTAINTY_FOV_BASELINE  Baseline temperature of an absolute-temperature metric
%
% Use as:
%   t0 = uncertainty_fov_baseline(name, temp_0)
%
% Input:
%   name   - metrics-table column name, e.g. 'maxT_brain' or 'maxT'
%   temp_0 - baseline temperatures [deg C]: struct with one field per tissue
%            (thermal.temp_0), or a scalar
%
% Output:
%   t0     - temp_0.<tissue> for a column ending in a tissue that temp_0 lists;
%            the scalar temp_0; else 37
%
% See also: UNCERTAINTY_FOV_MISMATCH, UNCERTAINTY_FOV_RESCALE

t0 = 37;
if isnumeric(temp_0) && isscalar(temp_0) && isfinite(temp_0)
    t0 = double(temp_0);
elseif isstruct(temp_0)
    tok = regexp(name, '^[A-Za-z0-9]+_(.+)$', 'tokens', 'once');
    if ~isempty(tok) && isfield(temp_0, tok{1}) && isnumeric(temp_0.(tok{1})) && ...
            isscalar(temp_0.(tok{1})) && isfinite(temp_0.(tok{1}))
        t0 = double(temp_0.(tok{1}));
    end
end
end
