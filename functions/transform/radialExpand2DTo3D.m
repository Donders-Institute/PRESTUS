function expanded3D = radialExpand2DTo3D(data2D, method)
% RADIALEXPAND2DTO3D  Radially expand a 2D axisymmetric half-plane into a 3D Cartesian volume
%
% Each axial slice of the input half-plane is interpolated onto a square
% Cartesian grid. The axis of symmetry (r=0) maps to voxel (Nr, Nr) in the
% output, placing one exact on-axis voxel with no half-voxel bias.
%
% Use as:
%   expanded3D = radialExpand2DTo3D(data2D)
%   expanded3D = radialExpand2DTo3D(data2D, method)
%
% Input:
%   data2D - [Nz x Nr] numeric, axisymmetric half-plane
%              Nz: number of axial slices; Nr: radial half-width
%              Column 1 must be the axis of symmetry (r = 0)
%   method - (char, optional) interp1 interpolation method (default: 'linear')
%              Use 'nearest' for integer label/mask fields to prevent
%              fractional blended values at tissue boundaries.
%
% Output:
%   expanded3D - [2*Nr x 2*Nr x Nz] numeric, radially expanded 3D volume
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, GRID_AXISYMMETRY

arguments
    data2D  (:,:) {mustBeNumeric}
    method  (1,:) char {mustBeMember(method, {'linear','nearest','cubic','spline'})} = 'linear'
end

% Extract size
[data_slices, radial_len] = size(data2D);

% Output grid size (square, double radial size)
output_size = 2 * radial_len;

% Create 2D grid. Center is placed at voxel (radial_len, radial_len) so
% that exactly one output voxel sits at r=0 (the axis of symmetry).
% Using radial_len (integer) instead of (output_size+1)/2 avoids the
% half-voxel bias that occurs with even output_size.
[xg, yg] = meshgrid(1:output_size, 1:output_size);
center = radial_len;
rg = sqrt((xg - center).^2 + (yg - center).^2);

% Radial vector for interpolation (0-based)
r_vec = 0:(radial_len - 1);

% Preallocate output
expanded3D = zeros(output_size, output_size, data_slices);

for c = 1:data_slices
    profile = data2D(c, :);
    % Interpolate onto 2D disk, fill outside radius with edge value.
    % 'nearest' is used for label/mask fields to preserve integer values.
    expanded2D = interp1(r_vec, profile, rg, method, profile(end));
    expanded3D(:, :, c) = expanded2D;
end

end