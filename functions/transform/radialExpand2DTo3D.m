function expanded3D = radialExpand2DTo3D(data2D, method, output_size, center)
% RADIALEXPAND2DTO3D  Radially expand a 2D axisymmetric half-plane into a 3D Cartesian volume
%
% Each axial slice of the input half-plane is interpolated onto a square
% Cartesian grid. By default the output is 2*Nr x 2*Nr and the axis lands
% at voxel (Nr, Nr). Pass output_size and center to target a specific
% output grid (e.g. to restore the original bilateral dimensions exactly).
%
% Use as:
%   expanded3D = radialExpand2DTo3D(data2D)
%   expanded3D = radialExpand2DTo3D(data2D, method)
%   expanded3D = radialExpand2DTo3D(data2D, method, output_size, center)
%
% Input:
%   data2D      - [Nz x Nr] numeric, axisymmetric half-plane
%                   Column 1 must be the axis of symmetry (r = 0)
%   method      - (char, optional) interp1 method (default: 'linear')
%                   Use 'nearest' for integer label/mask fields.
%   output_size - (scalar int, optional) side length of the square output
%                   (default: 2*Nr)
%   center      - (scalar int, optional) 1-based column/row index in the
%                   output that corresponds to r = 0 (default: Nr)
%
% Output:
%   expanded3D - [output_size x output_size x Nz] numeric
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, GRID_AXISYMMETRY

arguments
    data2D      (:,:) {mustBeNumeric}
    method      (1,:) char {mustBeMember(method, {'linear','nearest','cubic','spline'})} = 'linear'
    output_size (1,1) {mustBeInteger, mustBeNonnegative} = 0   % 0 = auto
    center      (1,1) {mustBeInteger, mustBeNonnegative} = 0   % 0 = auto
end

% Extract size
[data_slices, radial_len] = size(data2D);

% Resolve defaults
if output_size == 0; output_size = 2 * radial_len; end
if center      == 0; center      = radial_len;      end

% Build the 2D radial-distance grid relative to the axis voxel
[xg, yg] = meshgrid(1:output_size, 1:output_size);
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