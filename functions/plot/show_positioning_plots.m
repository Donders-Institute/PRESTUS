function show_positioning_plots(segmented_img_orig, t1_pixel_size, trans_pos_orig, focus_pos_orig, segmented_img_final, trans_pos_final, focus_pos_final, parameters, output_plot)

% SHOW_POSITIONING_PLOTS  Visualise transducer positioning before and after preprocessing
%
% Displays orthogonal slice overlays and a 3D view of the segmented head
% showing initial and final transducer/focus positions side by side.
% Optionally saves the figure. When multiple pairs are provided, only
% the first row is used.
%
% Use as:
%   show_positioning_plots(segmented_img_orig, t1_pixel_size, trans_pos_orig, ...
%       focus_pos_orig, segmented_img_final, trans_pos_final, focus_pos_final, ...
%       parameters, output_plot)
%
% Input:
%   segmented_img_orig  - [Nx x Ny x Nz] original tissue label volume
%   t1_pixel_size       - T1 voxel size [mm]
%   trans_pos_orig      - [nx3] original transducer position(s) in voxel coordinates
%   focus_pos_orig      - [nx3] original focus position(s) in voxel coordinates
%   segmented_img_final - [Nx x Ny x Nz] final (preprocessed) tissue label volume
%   trans_pos_final     - [nx3] final transducer position(s) in voxel coordinates
%   focus_pos_final     - [nx3] final focus position(s) in voxel coordinates
%   parameters          - (1,1) simulation parameters struct
%   output_plot         - path to save the output figure (empty string to skip)
%
% See also: SHOW_3D_HEAD, PLOT_CORONAL_SLICES

    arguments
        segmented_img_orig  (:,:,:) {mustBeNumeric}
        t1_pixel_size       (1,1)   double
        trans_pos_orig      (:,3)   double
        focus_pos_orig      (:,3)   double
        segmented_img_final (:,:,:) {mustBeNumeric}
        trans_pos_final     (:,3)   double
        focus_pos_final     (:,3)   double
        parameters          (1,1)   struct
        output_plot         (1,:)   char
    end

    %% Determine view settings based on transducer position
    view_pos = [-180, 0];
    slice_cap = [1, 0, 0];
    if trans_pos_final(1,3) > size(segmented_img_final, 3) / 2
        view_pos = [0, 0];
        slice_cap = [-1, 0, 0];
    end

    %% Generate coordinate mesh for visualization
    coord_mesh_xyz = get_xyz_mesh(segmented_img_orig);
    if gpuDeviceCount > 0
        coord_mesh_xyz = gpuArray(coord_mesh_xyz);
    end

    %% Create figure and generate subplots
    h = figure('Position', [200 200 1000 300]);

    % Subplot 1: Original segmentation with initial transducer and focus positions
    subplot(1, 3, 1);
    show_3d_head(segmented_img_orig, focus_pos_orig(1,:), trans_pos_orig(1,:), parameters, ...
                 t1_pixel_size, coord_mesh_xyz, [0 0 0], view_pos, 0);

    % Subplot 2: Original segmentation with slice cap applied
    subplot(1, 3, 2);
    show_3d_head(segmented_img_orig, focus_pos_orig(1,:), trans_pos_orig(1,:), parameters, ...
                 t1_pixel_size, coord_mesh_xyz, slice_cap, view_pos, 0);

    % Subplot 3: Final segmentation with updated transducer and focus positions
    ax3 = subplot(1, 3, 3);
    imagesc(squeeze(segmented_img_final(:, round(trans_pos_final(1,2)), :)));
    
    set(ax3,'dataAspectRatio',[1 1 1]);
    
    % Highlight final focus position
    rectangle('Position', [focus_pos_final(1,[3, 1]) - 2, 4, 4], ...
              'Curvature', [0, 0], 'EdgeColor', 'r', ...
              'LineWidth', 2,'LineStyle', '-');

    % Highlight final transducer position
    rectangle('Position', [trans_pos_final(1,[3, 1]) - 2, 4, 4], ...
              'Curvature', [0, 0], 'EdgeColor', 'b', ...
              'LineWidth', 2,'LineStyle', '-');

    % Draw line connecting transducer and focus positions
    line([trans_pos_final(1,3) focus_pos_final(1,3)], ...
         [trans_pos_final(1,1) focus_pos_final(1,1)], 'Color', 'white');
    
    % Add bounding box for visualization
    get_transducer_box(trans_pos_final(1,[1,3])', focus_pos_final(1,[1,3])', ...
        [], parameters.grid.resolution_mm, parameters);
    
    colormap(ax3,[0.3 0.3 0.3; lines(12)]);

    %% Save plot if output path is provided
    if output_plot ~= ""
        saveas(h ,output_plot ,'png');
        close(h);
    end

end