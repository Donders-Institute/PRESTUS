function [values, n_low, n_high] = pct_apply_limits(parameters, quantity, values)
% PCT_APPLY_LIMITS  Regularise CT HU or a CT-derived skull property to configured limits
%
% One scheme for the (pseudo-)CT and every property mapped from it:
% pct.limits.<quantity>_min and _max. A limit is a number, 'water' (the
% configured water value of that property) or [] (none). HU limits act on the
% CT before the skull is defined and the properties are mapped; property
% limits act on the mapped skull values.
%
% Use as:
%   [values, n_low, n_high] = pct_apply_limits(parameters, quantity, values)
%
% Input:
%   parameters - PRESTUS config (pct.limits, medium_properties.water)
%   quantity   - 'hu', 'density', 'sound_speed' or 'alpha_coeff'
%   values     - numeric array
%
% Output:
%   values - the array, clipped to [min, max]
%   n_low  - number of values raised to the minimum
%   n_high - number of values lowered to the maximum
%
% See also: PCT_MAP_PROPERTIES, PCT_HU_REFERENCE

    arguments
        parameters (1,1) struct
        quantity   (1,:) char {mustBeMember(quantity, {'hu', 'density', 'sound_speed', 'alpha_coeff'})}
        values     {mustBeNumeric}
    end
    lo = resolve(parameters, quantity, 'min');
    hi = resolve(parameters, quantity, 'max');
    if ~isempty(lo) && ~isempty(hi) && lo > hi
        error('pct_apply_limits:badLimits', 'pct.limits.%s_min (%g) exceeds pct.limits.%s_max (%g).', ...
            quantity, lo, quantity, hi);
    end
    n_low = 0; n_high = 0;
    if ~isempty(lo)
        below = values < lo;
        n_low = nnz(below);
        values(below) = lo;
    end
    if ~isempty(hi)
        above = values > hi;
        n_high = nnz(above);
        values(above) = hi;
    end
end

function v = resolve(parameters, quantity, side)
    v = get_field(parameters, {'pct', 'limits', [quantity '_' side]}, []);
    if isstring(v), v = char(v); end
    if ischar(v)
        if ~strcmp(v, 'water') || strcmp(quantity, 'hu')
            error('pct_apply_limits:badLimit', ['pct.limits.%s_%s must be a number, [] or ''water'' ' ...
                '(not for hu); got ''%s''.'], quantity, side, v);
        end
        v = parameters.medium_properties.water.(quantity);
    end
    if ~isempty(v) && ~(isscalar(v) && isnumeric(v) && ~isnan(v))
        error('pct_apply_limits:badLimit', 'pct.limits.%s_%s must be a scalar.', quantity, side);
    end
    v = double(v);
end
