function expanded3D = radialExpand2DTo3D(data2D)
% radialExpand2DTo3D - Radially expands 2D axisymmetric data into 3D Cartesian volume
%
% Syntax:
%   expanded3D = radialExpand2DTo3D(data2D)
%
% Input:
%   data2D - [C x Nr] matrix
%               C is number of slices/time points
%               Nr is radial dimension (index 1 = center),
%
% Output:
%   expanded3D - [2*Nr x 2*Nr x C] 3D array, each slice is radial expansion of data2D
%
% Example:
%   expanded3D = radialExpand2DTo3D(sensor_data.p_final);

% Extract size
[data_slices, radial_len] = size(data2D);

% Output grid size (square, double radial size)
output_size = 2 * radial_len;

% Create 2D grid centered at (output_size+1)/2
[xg, yg] = meshgrid(1:output_size, 1:output_size);
center = (output_size + 1) / 2;
rg = sqrt((xg - center).^2 + (yg - center).^2);

% Radial vector for interpolation (0-based)
r_vec = 0:(radial_len - 1);

% Preallocate output
expanded3D = zeros(output_size, output_size, data_slices);

for c = 1:data_slices
    profile = data2D(c, :); % 1 x 35
    % Interpolate onto 2D disk, fill outside radius with edge value
    expanded2D = interp1(r_vec, profile, rg, 'linear', profile(end));
    expanded3D(:, :, c) = expanded2D;
end

end