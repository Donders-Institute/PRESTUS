function [flhm_center, flhm_center_index] = get_flhm_center_position(x, y)
% GET_FLHM_CENTER_POSITION  Compute the centre of the focal length at half maximum (FLHM)
%
%   The threshold is (min(y) + max(y)) / 2. From the peak, the FLHM region
%   extends on each side to the last sample above the threshold before y
%   first drops to it; returns the midpoint of the two edge samples and the
%   index closest to it. Assumes y has a single peak.
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

    % Threshold halfway between minimum and maximum
    halfMax = (min(y) + max(y)) / 2;

    % Find the index of the maximum value in `y`
    [~, max_pos] = max(y);

    % Left edge: first sample after the last one at or below the threshold left of the peak
    index1 = find(y <= halfMax & x < x(max_pos), 1, 'last') + 1;

    % Right edge: last sample before y first drops to the threshold right of the peak
    index2 = find(y <= halfMax & x > x(max_pos), 1, 'first') - 1;

    % FLHM centre: midpoint between the edge samples
    flhm_center = (x(index2) - x(index1)) / 2 + x(index1);

    % Find the index in `x` closest to the computed FLHM center
    [~, flhm_center_index] = min(abs(x - flhm_center));
end
