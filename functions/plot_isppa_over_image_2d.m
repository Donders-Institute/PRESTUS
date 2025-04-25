function [h] = plot_isppa_over_image_2d(...
    Isppa_map, ...
    bg_image, ...
    transducer_bowl, ...
    after_exit_plane_mask, ...
    trans_pos, ...
    focus_pos, ...
    max_isppa_pos, ...
    options)

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
%   options         - Struct containing optional visualization settings:
%                     * show_rectangles: Boolean flag to show rectangles for key positions (default: 1).
%                     * rect_size: Size of rectangles for key positions (default: 2).
%                     * isppa_threshold_low/high: Thresholds for alpha scaling of ISppa map.
%                     * isppa_color_range: Range for ISppa map color scaling.
%                     * color_scale: Colormap for ISppa map (default: 'viridis').
%                     * show_colorbar: Boolean flag to display colorbar (default: 1).
% Output:
%   h                     - Handle to the created figure.

    arguments
        Isppa_map (:,:)
        bg_image (:,:)
        transducer_bowl (:,:)
        after_exit_plane_mask (:,:)
        trans_pos (:,2)
        focus_pos (:,2)
        max_isppa_pos (:,2)
        options.show_rectangles = 1
        options.rect_size = 2
        options.isppa_threshold_low (1,1) = min(Isppa_map(:))
        options.isppa_threshold_high (1,1) = min(Isppa_map(:)) + ...
            (max(Isppa_map(:)) - min(Isppa_map(:))) * 0.8
        options.isppa_color_range = []
        options.use_isppa_alpha (1,1) = 1
        options.color_scale = 'viridis'
        options.show_colorbar = 1
    end

    %% Set thresholds and color range for ISppa map
    if options.isppa_threshold_low == options.isppa_threshold_high
        options.isppa_threshold_low = options.isppa_threshold_low - 0.05;
    end
    if isempty(options.isppa_color_range)
        options.isppa_color_range = [max([options.isppa_threshold_low, min(Isppa_map(:)), 0]), max(Isppa_map(:))];
    end

    %% Normalize and prepare slices for visualization
    bg_slice = mat2gray(bg_image); % Normalize background image
    transducer_bowl = mat2gray(transducer_bowl); % Normalize transducer bowl mask
    before_exit_plane_mask = ~after_exit_plane_mask; % Compute mask for regions before exit plane

    %% Create figure and overlay images
    
    h = figure;
    ax1 = axes;
    imagesc(bg_slice + transducer_bowl + 0.2 * before_exit_plane_mask); % Overlay masks on background slice
    colormap(ax1, 'gray');
    axis image;
    axis off;

    ax2 = axes;

    if options.use_isppa_alpha
        isppa_alpha = rescale(Isppa_map, 'InputMin', options.isppa_threshold_low, 'InputMax', options.isppa_threshold_high);
    else
        isppa_alpha = ones(size(Isppa_map));
        isppa_alpha(Isppa_map==min(Isppa_map(:))) = 0;
    end
    
    imagesc(ax2, Isppa_map,'alphadata', isppa_alpha);

    if exist("clim")==2 % renamed in R2022a
        clim(options.isppa_color_range);
    elseif exist("caxis")==2
        caxis(options.isppa_color_range);
    end
    colormap(ax2, options.color_scale);

    ax2.Visible = 'off'; % Hide second axes visibility
    axis image;
    axis off;

    %% Highlight key positions with rectangles

    if options.show_rectangles 
        rect_size = options.rect_size;

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

    linkaxes([ax1, ax2], 'xy'); % Synchronize both axes to avoid shifts
    
    ax2.TightInset;
    ax1.TightInset;
    if options.show_colorbar
       C = colorbar(ax2);
       C.Position(4) = ax1.Position(4)*0.8;
       C.Position(2) = ax1.Position(2)+0.025;
       C.Position(1) = ax1.Position(1)+ax1.Position(3)+0.005;
       C.Color = 'black';
       C.FontSize = 8;
    end

end

