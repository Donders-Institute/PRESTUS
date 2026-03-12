function [outer_sphere_3d, segm_img_slice, trans_pos, geom_focus_pos, ex_plane_pos, norm_v] = ...
    tp_find_candidate_positions(img, target, pixel_size, parameters, subject_id, target_name)

% TP_FIND_CANDIDATE_POSITIONS Find candidate transducer positions on skull surface
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
%   subject_id    - Scalar subject ID for plotting
%   target_name   - String target name (e.g. 'motor_cortex') for plotting
%
% OUTPUT
%   outer_sphere_3d - 3D logical mask of valid transducer positions on skull surface
%   segm_img_slice  - RGB slice through target(y) with search sphere highlighted
%   trans_pos       - 1x3 initial random transducer position [x,y,z]
%   geom_focus_pos  - 1x3 geometric focus position
%   ex_plane_pos    - 1x3 exit plane position  
%   norm_v          - 1x3 unit normal vector from transducer to target

[t1_x, t1_y, t1_z] = ndgrid(1:size(img));

% Expand sphere until skull intersection found
if ~isfield(parameters, 'min_focal_distance_mm')
    parameters.min_focal_distance_mm = parameters.expected_focal_distance_bowl;
end
outer_sphere_3d = []; outer_sphere = [];
while numel(find(outer_sphere)) < 1
    grid_dist = sqrt((t1_x-target(1)).^2+(t1_y-target(2)).^2+(t1_z-target(3)).^2);
    dist_sphere = abs(grid_dist-parameters.min_focal_distance_mm/pixel_size)<0.5;
    outer_sphere_3d = dist_sphere&img==0;
    segm_img_slice = ind2rgb(squeeze(img(:,target(2),:)), viridis(max(img(:))+1));
    outer_sphere = squeeze(dist_sphere(:,target(2),:))&squeeze(img(:,target(2),:))==0;
    segm_img_slice(outer_sphere) = 1;
    parameters.min_focal_distance_mm = parameters.min_focal_distance_mm + 3;
end

