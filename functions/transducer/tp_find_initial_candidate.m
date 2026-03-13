function [trans_candidate, outer_sphere_3d, parameters] = ...
    tp_find_initial_candidate(img, target, pixel_size, parameters)

% TP_FIND_INITIAL_CANDIDATE Find candidate transducer positions on skull surface
%
% Finds candidate transducer positions by expanding a sphere from target until it intersects 
% skull exterior surface (air-tissue boundary). Generates initial positioning plots and computes
% geometric properties (focus, exit plane, normal vector).
%
% INPUT
%   img           - Segmented head image (nifti final_tissues.nii.gz)
%   target        - 1x3 target coordinates [x,y,z] in voxel space
%   pixel_size    - Scalar voxel size in mm (mean of PixelDimensions)
%   parameters    - Struct with transducer properties (.transducer.curv_radius_mm, etc.)
%
% OUTPUT
%   trans_candidate - Initial random transducer position [x,y,z]
%   outer_sphere_3d - 3D logical mask of valid transducer positions on skull surface
%   parameters      - Parameters with updated min_focal_distance_mm

disp("[TP] Finding candidate position via intersection of skull with expanding sphere ...")

sz = size(img);  % Get image dimensions [Nx Ny Nz]
[t1_x, t1_y, t1_z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));  % Create 3D voxel coordinate grids

% Initialize minimum focal distance (bowl radius in voxels)
if ~isfield(parameters, 'min_focal_distance_mm')
    parameters.min_focal_distance_mm = parameters.expected_focal_distance_bowl;  % Default to bowl geometry
end

outer_sphere_3d = []; outer_sphere = [];  % Pre-allocate for while loop
while numel(find(outer_sphere)) < 1  % Expand until sphere intersects skull exterior (img==0)
    
    % Compute Euclidean distance from target voxel to every voxel in volume
    grid_dist = sqrt((t1_x-target(1)).^2 + (t1_y-target(2)).^2 + (t1_z-target(3)).^2);
    
    % Create binary sphere shell: voxels within ±0.5 voxel of current radius
    dist_sphere = abs(grid_dist - parameters.min_focal_distance_mm/pixel_size) < 0.5;
    
    % 3D: Find skull exterior intersection (sphere shell & outside brain/skull)
    outer_sphere_3d = dist_sphere & (img == 0);
    
    % Extract & visualize central Y-slice through target for debugging
    img_slice = ind2rgb(squeeze(img(:,target(2),:)), viridis(max(img(:))+1));  % Colormap segmentation
    outer_sphere = squeeze(dist_sphere(:,target(2),:)) & squeeze(img(:,target(2),:)) == 0;  % 2D slice intersection
    
    % Overlay white intersection points on slice (for visual verification)
    img_slice(outer_sphere) = 1;  % RGB [1 1 1] = white
    
    % Increment radius by 3mm for next iteration
    parameters.min_focal_distance_mm = parameters.min_focal_distance_mm + 3;
end

% Random transducer candidate from skull surface intersection
outer_idx = find(outer_sphere_3d&t1_y==target(2));
trans_idx = randsample(outer_idx, 1);
trans_pos = [t1_x(trans_idx), t1_y(trans_idx), t1_z(trans_idx)];
trans_xz = trans_pos([1,3]);

trans_candidate.outer_idx = outer_idx;
trans_candidate.trans_idx = trans_idx;
trans_candidate.trans_pos = trans_pos;
trans_candidate.trans_xz = trans_xz;