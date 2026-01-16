function show_positioning_plots(segmented_img_orig, t1_pixel_size, trans_pos_orig, focus_pos_orig, segmented_img_final, trans_pos_final, focus_pos_final, parameters, output_plot)

% SHOW_POSITIONING_PLOTS Visualizes transducer positioning before and after preprocessing.
%
% This function generates visualizations of transducer positioning on segmented brain images. 
% It displays the original segmentation with the initial transducer and focus positions, 
% as well as the final segmentation with updated positions after preprocessing. 
% The function also includes a 3D visualization of the transducer placement.
%
% Input:
%   segmented_img_orig - [Nx x Ny x Nz] matrix representing the original segmented image.
%   t1_pixel_size      - Scalar specifying the pixel size of the T1 image (in mm).
%   trans_pos_orig     - [nx3] array specifying the original transducer position in voxel coordinates.
%   focus_pos_orig     - [nx3] array specifying the original focus position in voxel coordinates.
%   segmented_img_final- [Nx x Ny x Nz] matrix representing the final segmented image after preprocessing.
%   trans_pos_final    - [nx3] array specifying the final transducer position in voxel coordinates.
%   focus_pos_final    - [nx3] array specifying the final focus position in voxel coordinates.
%   parameters         - Struct containing simulation parameters (e.g., grid step size).
%   output_plot        - String specifying the path to save the output plot (if empty, no plot is saved).
%
% Output:
%   None. The function generates visualizations and optionally saves them to a file.
% 
% Notes:
%   - If multiple transducer–focus pairs are provided (n > 1), only the first pair
%     (row 1) is used for visualization in all panels.

    arguments
        segmented_img_orig (:,:,:) 
        t1_pixel_size (1,1)
        trans_pos_orig (:,3)
        focus_pos_orig (:,3)
        segmented_img_final (:,:,:) 
        trans_pos_final (:,3)
        focus_pos_final (:,3)
        parameters struct
        output_plot string
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
    get_transducer_box(trans_pos_final(1,[1,3])', focus_pos_final(1,[1,3])', parameters.grid_step_mm, parameters);
    
    colormap(ax3,[0.3 0.3 0.3; lines(12)]);

    %% Save plot if output path is provided
    if output_plot ~= ""
        saveas(h ,output_plot ,'png');
        close(h);
    end

end