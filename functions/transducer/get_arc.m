function [x, y] = get_arc(center, radius, angleStart, angleEnd, npoints)

% GET_ARC  Generate x/y coordinates of a 2-D arc
%
% Computes uniformly sampled points along a circular arc defined by a
% centre, radius, start angle, and end angle.
%
% Use as:
%   [x, y] = get_arc(center, radius, angleStart, angleEnd)
%   [x, y] = get_arc(center, radius, angleStart, angleEnd, npoints)
%
% Input:
%   center     - [1x2] (x, y) centre of the arc
%   radius     - arc radius
%   angleStart - start angle [rad]
%   angleEnd   - end angle [rad]
%   npoints    - number of sample points (default: 200)
%
% Output:
%   x - [1xnpoints] x-coordinates along the arc
%   y - [1xnpoints] y-coordinates along the arc
%
% See also: GET_TRANSDUCER_BOX

    arguments
        center (1,2) double
        radius (1,1) double
        angleStart (1,1) double
        angleEnd (1,1) double
        npoints (1,1) double = 200
    end

    % Generate angles uniformly between `angleStart` and `angleEnd`
    theta = linspace(angleStart, angleEnd, npoints);

    % Compute x and y coordinates for each angle
    x = radius * cos(theta) + center(1); 
    y = radius * sin(theta) + center(2); 
end
