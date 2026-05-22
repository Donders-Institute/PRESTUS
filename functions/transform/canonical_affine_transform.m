function nii_out = canonical_affine_transform(nii_in_path, output_path)
% CANONICAL_AFFINE_TRANSFORM  Set NIfTI affine to a diagonal voxel-scaling matrix (LPI origin)
%
% Removes rotations and shear from the NIfTI affine, replacing it with a
% pure diagonal scaling matrix whose diagonal entries equal the voxel
% dimensions in mm. The origin is fixed at [0 0 0] (LPI corner). This
% prevents k-Wave grid-padding issues caused by rotated affines and ensures
% Localite neuronavigation can import the image using simple voxel scaling
% (world_pos = voxel_pos * voxel_size) without a rotation matrix.
% Equivalent to Python: nib.as_closest_canonical() + custom diagonal affine.
%
% Use as:
%   nii_out = canonical_affine_transform(nii_in_path, output_path)
%
% Input:
%   nii_in_path - path to the input T1 NIfTI file (.nii or .nii.gz)
%   output_path - path for the output NIfTI file (uncompressed;
%                 a .gz copy is also written automatically)
%
% Output:
%   nii_out     - modified niftiinfo struct with diagonal affine
%                 (also saved to output_path and output_path.gz)
%
% See also: NIFTIINFO, NIFTIREAD, NIFTIWRITE, RAS_TO_GRID

arguments
    nii_in_path (1,:) char
    output_path (1,:) char
end

%% Load input NIfTI
nii_in = niftiinfo(nii_in_path);
img_data = niftiread(nii_in);

fprintf('Original affine:\n'); disp(nii_in.Transform.T);
fprintf('Shape: [%d %d %d]\n', size(img_data));
fprintf('LOCALITE NOTE: Original affine may contain rotations causing coordinate mismatches\n');

%% Print original affine
nii_original = niftiinfo(nii_in_path);
fprintf('Original affine:\n'); disp(nii_original.Transform.T);

%% Create new diagonal affine matrix for Localite compatibility
% [sx 0  0  0 ]    LOCALITE: world_x = voxel_x * sx
% [0  sy 0  0 ]    LOCALITE: world_y = voxel_y * sy  
% [0  0  sz 0 ]    LOCALITE: world_z = voxel_z * sz
% [0  0  0  1 ]
affine_new = eye(4);
affine_new(1,1) = nii_original.PixelDimensions(1);   % X scaling (mm/voxel)
affine_new(2,2) = nii_original.PixelDimensions(2);   % Y scaling (mm/voxel)  
affine_new(3,3) = nii_original.PixelDimensions(3);   % Z scaling (mm/voxel)

fprintf('New LOCALITE-compatible diagonal affine:\n'); disp(affine_new);
fprintf('Mapping: [voxel_x, voxel_y, voxel_z] → [%.1f*vx, %.1f*vy, %.1f*vz] mm\n', ...
        affine_new(1,1), affine_new(2,2), affine_new(3,3));

%% Create output NIfTI structure
nii_out = nii_in;
nii_out.Transform.T = affine_new;
nii_out.Filename = output_path;
nii_out.ImageSize = size(img_data);  % Preserve dimensions

%% Write output file
niftiwrite(img_data, output_path, nii_out);
gzip(output_path); % also save compressed file
fprintf('Saved LOCALITE-ready T1: %s\n', output_path);
fprintf('Ready for: Localite → k-Wave → PRESTUS simulation pipeline\n');

end
