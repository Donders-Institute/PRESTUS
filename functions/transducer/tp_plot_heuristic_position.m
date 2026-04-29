function tp_plot_heuristic_position(trans_pos, target_pos, img, img_header, parameters, pixel_size, target_name, subject_id)
% TP_PLOT_HEURISTIC_POSITION  Create visualisation of heuristic transducer placement
%
% Rotates the head image to the focal axis, then generates slice-overlay
% plots showing the heuristic transducer and target positions. Saves
% figures to disk.
%
% Use as:
%   tp_plot_heuristic_position(trans_pos, target_pos, img, img_header, ...
%       parameters, pixel_size, target_name, subject_id)
%
% Input:
%   trans_pos   - [1x3] heuristic transducer position in voxel space
%   target_pos  - [1x3] target position in voxel space
%   img         - [Nx x Ny x Nz] segmented head volume
%   img_header  - NIfTI header struct (from niftiinfo)
%   parameters  - (1,1) simulation parameters struct
%   pixel_size  - voxel size [mm]
%   target_name - target label string used in output filename
%   subject_id  - subject identifier used in output filename
%
% See also: TP_SELECT_HEURISTIC_POSITION, TP_PLOT_CANDIDATE_POSITIONS

trans_xyz = gather(trans_pos);
target_xyz = gather(target_pos);
img = gather(img);

sz = size(img);
[t1_x, t1_y, t1_z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
coord_mesh.xyz = [reshape(t1_x,[],1), reshape(t1_y,[],1), reshape(t1_z,[],1)];

%% Align to focal axis 

[img_rotated, trans_xyz_rotated, target_xyz_rotated, transformation_matrix, ~, ~, ~, ~] = ...
    preproc_align_to_focal_axis(...
    img, img_header, round(trans_xyz)', target_xyz', 1, parameters);

%% Visualization

h = figure('units', 'normalized', 'position', [0 0 1 .5]);
subplot(1,3,1);
    colormap([0.3 0.3 0.3; lines(5)]);
    imagesc(squeeze(img_rotated(:,round(trans_xyz_rotated(2)),:)));
    
    % TARGET - Red box
    rectangle('Position', [target_xyz_rotated([3,1]) - 2, 4, 4], ...
                'Curvature', [0,0], ...
                'EdgeColor', 'r', ...
                'LineWidth', 2, 'LineStyle', '-');
    
    % TRANSDUCER (heuristic position) - Blue box  
    rectangle('Position', [trans_xyz_rotated([3,1]) - 2, 4, 4],...
                'Curvature', [0,0],...
                'EdgeColor', 'b', ...
                'LineWidth', 2, 'LineStyle', '-');

    line([trans_xyz_rotated(3) target_xyz_rotated(3)], ...
        [trans_xyz_rotated(1) target_xyz_rotated(1)], 'Color', 'white');

    get_transducer_box(...
        trans_xyz_rotated([1,3]), target_xyz_rotated([1,3]), [], pixel_size, parameters);

subplot(1,3,2);

    show_3d_head(img, ...
        target_xyz, ...
        trans_xyz, ...
        parameters, ...
        pixel_size, ...
        coord_mesh.xyz, ...
        [0 0 0],...
        [0,0],...
        0)

subplot(1,3,3);

    show_3d_head(img, ...
        target_xyz, ...
        trans_xyz, ...
        parameters, ...
        pixel_size, ...
        coord_mesh.xyz, ...
        [0 0 0],...
        [0,0],...
        0)

    view([-175,0])

% save figure in output directory 
output_plot = fullfile(parameters.io.dir_img,...
    sprintf('sub-%03d_heuristic_%s.png', subject_id, target_name));
saveas(h, output_plot, 'png')

%% [Optional] Save copy in localite directory
if isfield(parameters.path, 'localite') && ~isempty(parameters.path.localite)
    saveas(h, fullfile(parameters.path.localite, ...
        sprintf("sub-%03.0f_%s", subject_id, target_name)), 'png');
end

% close figure
close(h);
