function [x, y] = getArc(center, radius, angleStart, angleEnd, npoints)

% GETARC Generates the coordinates of an arc in 2D space.
%
% This function computes the x and y coordinates of points along an arc 
% defined by a center, radius, start angle, and end angle. The arc is 
% sampled uniformly with a specified number of points.
%
% Input:
%   center     - [1x2] array specifying the (x, y) coordinates of the arc's center.
%   radius     - Scalar specifying the radius of the arc.
%   angleStart - Scalar specifying the starting angle of the arc (in radians).
%   angleEnd   - Scalar specifying the ending angle of the arc (in radians).
%   npoints    - Integer specifying the number of points to sample along the arc (default: 200).
%
% Output:
%   x          - [1xnpoints] array containing the x coordinates of points along the arc.
%   y          - [1xnpoints] array containing the y coordinates of points along the arc.

    % Generate angles uniformly between `angleStart` and `angleEnd`
    theta = linspace(angleStart, angleEnd, npoints);

    % Compute x and y coordinates for each angle
    x = radius * cos(theta) + center(1); 
    y = radius * sin(theta) + center(2); 
end
