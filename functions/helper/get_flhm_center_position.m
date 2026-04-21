function [flhm_center, flhm_center_index] = get_flhm_center_position(x, y)
% GET_FLHM_CENTER_POSITION  Compute the centre of the full-length half-maximum (FLHM)
%
%   Locates the FLHM region (where y >= half-max) and returns the midpoint
%   between the left and right FLHM crossing positions, along with its index.
%   Assumes y has a single peak with monotonically decreasing flanks.
%
% Use as:
%   [flhm_center, flhm_center_index] = get_flhm_center_position(x, y)
%
% Input:
%   x - [1xN] independent variable (e.g., spatial coordinate) [mm]
%   y - [1xN] dependent variable (e.g., pressure or intensity profile)
%
% Output:
%   flhm_center       - [1x1] coordinate value at the FLHM centre
%   flhm_center_index - [1x1] index into x closest to flhm_center
%
% See also: GET_SLICE_BY_LABEL, MASKED_MAX_3D

arguments
    x (1,:) {mustBeNumeric}
    y (1,:) {mustBeNumeric}
end

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
