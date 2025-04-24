function [h] = plot_isppa_over_image_2d(...
    Isppa_map, ...
    bg_image, ...
    transducer_bowl, ...
    after_exit_plane_mask, ...
    trans_pos, focus_pos, ...
    max_isppa_pos)

% PLOT_ISPPA_OVER_IMAGE_2D Visualizes ISppa map overlaid on a background image.
%
% This function overlays the spatial peak pulse-average intensity (ISppa) map 
% on a 2D background image (`bg_image`) and highlights key positions such as 
% the transducer position, focus position, and maximum ISppa position. The 
% visualization includes additional masks for regions before and after the 
% transducer's exit plane.
%
% Input:
%   Isppa_map             - [Nx x Ny] matrix representing the ISppa intensity map.
%   bg_image              - [Nx x Ny] matrix representing the background image (e.g., anatomical image).
%   transducer_bowl       - [Nx x Ny] binary mask representing the transducer bowl region.
%   after_exit_plane_mask - [Nx x Ny] binary mask for regions after the transducer's exit plane.
%   trans_pos             - [1x2] array specifying the transducer position in grid coordinates (row, col).
%   focus_pos             - [1x2] array specifying the focus position in grid coordinates (row, col).
%   max_isppa_pos         - [1x2] array specifying the maximum ISppa position in grid coordinates (row, col).
%
% Output:
%   h                     - Handle to the created figure.

    %% Normalize and prepare slices for visualization
    bg_slice = mat2gray(bg_image); % Normalize background image
    transducer_bowl = mat2gray(transducer_bowl); % Normalize transducer bowl mask
    before_exit_plane_mask = ~after_exit_plane_mask; % Compute mask for regions before exit plane

    %% Create figure and overlay images
    h = figure; % Create figure
    ax1 = axes; % First axes for background image
    imagesc(bg_slice + transducer_bowl + 0.2 * before_exit_plane_mask); % Overlay masks on background slice
    colormap(ax1, 'gray'); % Use grayscale colormap for background
    axis image;
    axis off;

    ax2 = axes; % Second axes for ISppa map
    imagesc(ax2, Isppa_map, 'alphadata', (Isppa_map - min(Isppa_map)) / max(Isppa_map(:))); % Overlay ISppa map with transparency
    colormap(ax2, 'viridis'); % Use viridis colormap for ISppa map
    ax2.Visible = 'off'; % Hide second axes visibility

    axis image;
    axis off;

    %% Highlight key positions with rectangles
    rect_size = 2; % Size of rectangles

    % Transducer position (green rectangle)
    rectangle('Position', [trans_pos(2) - rect_size / 2, trans_pos(1) - rect_size / 2, rect_size * 2 + 1, rect_size * 2 + 1], ...
              'EdgeColor', 'g', 'LineWidth', 1, 'LineStyle', '-');

    % Focus position (red rectangle)
    rectangle('Position', [focus_pos(2) - rect_size / 2, focus_pos(1) - rect_size / 2, rect_size * 2 + 1, rect_size * 2 + 1], ...
              'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '-');

    % Maximum ISppa position (blue rectangle)
    rectangle('Position', [max_isppa_pos(2) - rect_size / 2, max_isppa_pos(1) - rect_size / 2, rect_size * 2 + 1, rect_size * 2 + 1], ...
              'EdgeColor', 'b', 'LineWidth', 1, 'LineStyle', '-');

end

