function [flhm_center, flhm_center_index] = get_flhm_center_position(x, y)

% GET_FLHM_CENTER_POSITION Computes the center position of the full-length half-maximum (FLHM).
%
% This function calculates the center position of the FLHM for a given dataset 
% defined by `x` (independent variable) and `y` (dependent variable). The FLHM 
% is defined as the region where `y` values are greater than or equal to half 
% the maximum value of `y`.
%
% Input:
%   x - [1xN] array of independent variable values.
%   y - [1xN] array of dependent variable values corresponding to `x`.
%
% Output:
%   flhm_center        - Scalar value representing the center position of the FLHM.
%   flhm_center_index  - Index of `x` corresponding to the FLHM center.
%
% Notes:
%   - The function assumes that `y` has a single peak and monotonically decreases 
%     on both sides of the peak.
%   - The FLHM is calculated as the midpoint between the first point below half-max 
%     on the left side and the last point below half-max on the right side.

    % Compute half of the maximum range (half-max value)
    halfMax = (min(y) + max(y)) / 2;

    % Find the index of the maximum value in `y`
    [~, max_pos] = max(y);

    % Find where `y` first drops below half-max on the left side of the peak
    index1 = find(y <= halfMax & x < x(max_pos), 1, 'last') + 1;

    % Find where `y` last rises above half-max on the right side of the peak
    index2 = find(y <= halfMax & x > x(max_pos), 1, 'first') - 1;

    % Compute FLHM center as the midpoint between `index1` and `index2`
    flhm_center = (x(index2) - x(index1)) / 2 + x(index1);

    % Find the index in `x` closest to the computed FLHM center
    [~, flhm_center_index] = min(abs(x - flhm_center));
end
