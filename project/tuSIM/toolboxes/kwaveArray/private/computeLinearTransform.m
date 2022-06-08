function [R, b] = computeLinearTransform(bowl_pos, focus_pos, radius)
%COMPUTELINEARTRANSFORM Compute a linear transformation.
%
% DESCRIPTION:
%     computeLinearTransform calculates a rotation matrix to tranform the
%     computed bowl (or disc) points to the orientation specified by the
%     bowl and focus positions.
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 5th December 2016
%     last update - 4th February 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2016-2018 Elliott Wise and Bradley Treeby
%
% See also makeCartBowl, makeCartDisc

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
%
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
% more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% compute vector pointing from bowl_pos to focus_pos
beam_vec = focus_pos - bowl_pos;

% normalise to give unit beam vector
beam_vec = beam_vec ./ norm(beam_vec); 

% canonical normalised beam_vec (canonical bowl_pos is [0, 0, 1])
beam_vec0 = [0, 0, -1];                

% find the rotation matrix for the bowl
u = cross(beam_vec0, beam_vec);

% normalise the rotation matrix if not zero
if any(u ~= 0)
    u = u.' ./ norm(u);
end

% find the axis-angle transformation between beam_vec and e1
theta = acos(dot(beam_vec0, beam_vec));
    
% convert axis-angle transformation to a rotation matrix
% https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
A = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
R = cos(theta) * eye(3) + sin(theta) * A + (1 - cos(theta)) * (u * u.');
    
% compute an offset for the bowl, where bowl_centre = move from bowl_pos
% towards focus by radius
if nargin == 3
    b = bowl_pos.' + radius * beam_vec.';
else
    b = 0;
end