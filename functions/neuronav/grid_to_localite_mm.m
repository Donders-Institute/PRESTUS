function loc_mm = grid_to_localite_mm(grid_pos, img_header)
% GRID_TO_LOCALITE_MM  Voxel indices to coordinates in the Localite planning image
%
% The heuristic placement writes a Localite planning image with
% CANONICAL_AFFINE_TRANSFORM: a diagonal affine with zero origin, so voxel
% (i,j,k) (1-based) lies at ((i,j,k) - 1) .* voxel size [mm]. This returns
% that position, rounded to 0.01 mm. It assumes that Localite shows the
% world coordinates of the planning image's NIfTI affine.
%
% Use as:
%   loc_mm = grid_to_localite_mm(grid_pos, img_header)
%
% Input:
%   grid_pos   - [Nx3] 1-based voxel indices of the T1 the planning image was made from
%   img_header - niftiinfo header of that T1 (PixelDimensions)
%
% Output:
%   loc_mm     - [Nx3] coordinates in the planning image [mm]
%
% See also: CANONICAL_AFFINE_TRANSFORM, GRID_TO_RAS, TP_SELECT_HEURISTIC_POSITION

arguments
    grid_pos   (:,3) {mustBeNumeric}
    img_header (1,1) struct
end

    loc_mm = round((grid_pos - 1) .* img_header.PixelDimensions(1:3), 2);
end
