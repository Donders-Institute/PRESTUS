function [h] = plot_overlay_2d(...
    overlay_image, ...
    bg_image, ...
    transducer_bowl, ...
    after_exit_plane_mask, ...
    trans_pos, ...
    focus_pos, ...
    max_intensity_pos, ...
    options)

% PLOT_OVERLAY_2D  Overlay a metric map on a 2D background image
%
% Overlays a 2D metric image on an anatomical background and marks
% transducer, focus, and maximum intensity positions. Supports alpha-scaled
% colouring and optional post-exit-plane masking.
%
% Use as:
%   h = plot_overlay_2d(overlay_image, bg_image, transducer_bowl, ...
%       after_exit_plane_mask, trans_pos, focus_pos, max_intensity_pos)
%   h = plot_overlay_2d(..., options)
%
% Input:
%   overlay_image         - [Nx x Ny] metric map to overlay
%   bg_image              - [Nx x Ny] anatomical background image
%   transducer_bowl       - [Nx x Ny] binary transducer bowl mask
%   after_exit_plane_mask - [Nx x Ny] binary mask for post-exit-plane region
%   trans_pos             - [1x2] transducer position in grid coordinates
%   focus_pos             - [1x2] focus position in grid coordinates
%   max_intensity_pos     - [1x2] maximum intensity position in grid coordinates
%   options               - name-value visualisation settings
%
% Output:
%   h - figure handle
%
% See also: PLOT_OVERLAY, PLOT_T1_WITH_TRANSDUCER

    arguments
        overlay_image (:,:)
        bg_image (:,:)
        transducer_bowl (:,:)
        after_exit_plane_mask (:,:)
        trans_pos (:,2)
        focus_pos (:,2)
        max_intensity_pos (:,2)
        options.show_rectangles = 1
        options.rect_size = 2
        options.overlay_threshold_low (1,1) = min(overlay_image(:))
        options.overlay_threshold_high (1,1) = min(overlay_image(:)) + ...
            (max(overlay_image(:)) - min(overlay_image(:))) * 0.8
        options.overlay_color_range = []
        options.bg_bw_range = []
        options.use_overlay_alpha (1,1) = 1
        options.color_scale = 'viridis'
        options.show_colorbar = 1
    end

    %% Set thresholds and color range for intensity map
    if options.overlay_threshold_low == options.overlay_threshold_high
        options.overlay_threshold_low = options.overlay_threshold_low - 0.05;
    end
    if isempty(options.overlay_color_range)
        options.overlay_color_range = [max([options.overlay_threshold_low, min(overlay_image(:)), 0]), max(overlay_image(:))];
    end
    if isempty(options.bg_bw_range)
        options.bg_bw_range = [0 max(bg_image(:))];
    end

    %% Normalize and prepare slices for visualization
    
    bg_slice = mat2gray(bg_image); % Normalize background image
    transducer_bowl = mat2gray(transducer_bowl); % Normalize transducer bowl mask
    before_exit_plane_mask = ~after_exit_plane_mask; % Compute mask for regions before exit plane

    %% Create figure and overlay images
    
    h = figure;
    ax1 = axes;
    if options.bg_bw_range(2)>0
        imagesc(bg_slice + transducer_bowl + 0.2 * before_exit_plane_mask, options.bg_bw_range); % Overlay masks on background slice
    else
        imagesc(bg_slice + transducer_bowl + 0.2 * before_exit_plane_mask);
    end
    colormap(ax1, 'gray');
    axis image;
    axis off;

    ax2 = axes;

    if options.use_overlay_alpha
        intensity_alpha = rescale(overlay_image, 'InputMin', options.overlay_threshold_low, 'InputMax', options.overlay_threshold_high);
    else
        intensity_alpha = ones(size(overlay_image));
        intensity_alpha(overlay_image==min(overlay_image(:))) = 0;
    end
    
    imagesc(ax2, overlay_image,'alphadata', intensity_alpha);

    if exist("clim")==2 % renamed in R2022a
        clim(options.overlay_color_range);
    elseif exist("caxis")==2
        caxis(options.overlay_color_range);
    end
    % use requested color scale (or revert to MATLAB default)
    try
        colormap(ax2, options.color_scale);
    catch
        colormap(ax2, 'parula');  % Built-in MATLAB default
        warning('Using parula instead of %s', options.color_scale);
    end

    ax2.Visible = 'off'; % Hide second axes visibility
    axis image;
    axis off;

    %% Highlight key positions with rectangles

    if options.show_rectangles 
        rect_size = options.rect_size;

        % if axisymmetric setup (transducer centered at y=1):
        % plot half the horizontal box
        if trans_pos(2) == 1
            rect_size_horizontal = rect_size + 0.5;
        else
            rect_size_horizontal = rect_size * 2 + 1;
        end

        % Transducer position (green rectangle)
        rectangle(...
            'Position', [trans_pos(2) - rect_size / 2, ...
                trans_pos(1) - rect_size / 2, ...
                rect_size_horizontal, ...
                rect_size * 2 + 1], ...
            'EdgeColor', 'g', 'LineWidth', 1, 'LineStyle', '-');

        % Focus position (red rectangle)
        rectangle(...
            'Position', [focus_pos(2) - rect_size / 2, ...
                focus_pos(1) - rect_size / 2, ...
                rect_size_horizontal, ...
                rect_size * 2 + 1], ...
            'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '-');

        % Maximum intensity position (blue rectangle)
        rectangle(...
            'Position', [max_intensity_pos(2) - rect_size / 2, ...
                max_intensity_pos(1) - rect_size / 2, ...
                rect_size_horizontal, ...
                rect_size * 2 + 1], ...
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

