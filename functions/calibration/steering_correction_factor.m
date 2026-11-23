function k = steering_correction_factor(sc, targets_mm)
% STEERING_CORRECTION_FACTOR  Measured / modelled focal peak at steering targets
%
% Evaluates a steering correction fitted by FIT_STEERING_CORRECTION (or
% read from the transducer library). A quadratic map is not extrapolated:
% a target outside the measured range is moved onto its edge, coordinate
% by coordinate, with a warning, since a quadratic surface can run away
% outside the data it was fitted to.
%
% Use as:
%   k = steering_correction_factor(sc, targets_mm)
%
% Input:
%   sc         - steering correction struct: model ('scale' | 'quadratic'),
%                scale, and for 'quadratic' coef [10], mean_mm, sd_mm [3],
%                bounds_mm [2x3]
%   targets_mm - [N x 3] targets [mm], transducer frame
%
% Output:
%   k - [N x 1] factor (multiply the modelled drive or peak by k)
%
% See also: FIT_STEERING_CORRECTION, CREATE_MATRIX_KARRAY

arguments
    sc (1,1) struct
    targets_mm {mustBeNumeric}
end

    N = size(targets_mm, 1);
    switch char(sc.model)
        case 'scale'
            k = repmat(sc.scale, N, 1);
        case 'quadratic'
            lo = sc.bounds_mm(1, :); hi = sc.bounds_mm(2, :);
            X = min(max(targets_mm, lo), hi);
            if any(abs(X(:) - targets_mm(:)) > 1e-9)
                warning('steering_correction_factor:outside', ...
                    ['Steering target outside the calibrated range ' ...
                     '(x %.1f..%.1f, y %.1f..%.1f, z %.1f..%.1f mm): the correction at the nearest edge is used.'], ...
                    lo(1), hi(1), lo(2), hi(2), lo(3), hi(3));
            end
            k = steering_quad_terms((X - sc.mean_mm(:)') ./ sc.sd_mm(:)') * sc.coef(:);
        otherwise
            error('steering_correction_factor:model', 'Unknown steering correction model ''%s''.', char(sc.model));
    end
end
