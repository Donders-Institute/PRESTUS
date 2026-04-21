function [bg_slice, transducer_bowl, overlay_image, ax1, ax2, bg_min, bg_max, h] = ...
    plot_overlay(overlay_image, bg_image, transducer_bowl, parameters, slice, trans_pos, focus_pos, max_data_pos, options)

% PLOT_OVERLAY  Overlay a metric map on a 2D slice of a 3D background image
%
% Extracts a specified orthogonal slice from a 3D volume and overlays a
% metric (e.g., acoustic intensity) with alpha-scaled colouring. Marks
% transducer, focus, and maximum intensity positions.
%
% Use as:
%   [bg_slice, transducer_bowl, overlay_image, ax1, ax2, bg_min, bg_max, h] = ...
%       plot_overlay(overlay_image, bg_image, transducer_bowl, parameters, ...
%           slice, trans_pos, focus_pos, max_data_pos)
%   [bg_slice, transducer_bowl, overlay_image, ax1, ax2, bg_min, bg_max, h] = ...
%       plot_overlay(..., options)
%
% Input:
%   overlay_image   - [Nx x Ny x Nz] metric volume to overlay
%   bg_image        - [Nx x Ny x Nz] anatomical background volume
%   transducer_bowl - [Nx x Ny x Nz] binary transducer bowl mask
%   parameters      - (1,1) simulation parameters struct
%   slice           - {axis, index} cell specifying the slice ('x', 'y', or 'z')
%   trans_pos       - [1x3] transducer position in grid coordinates
%   focus_pos       - [1x3] focus position in grid coordinates
%   max_data_pos    - [1x3] maximum intensity position in grid coordinates
%   options         - name-value visualisation settings (show_rectangles,
%                     grid_step, rect_size, overlay thresholds, color_scale, etc.)
%
% Output:
%   bg_slice        - extracted 2D background slice
%   transducer_bowl - extracted 2D transducer bowl slice
%   overlay_image   - extracted 2D overlay slice
%   ax1             - handle to background axes
%   ax2             - handle to overlay axes
%   bg_min          - minimum background intensity
%   bg_max          - maximum background intensity
%   h               - figure handle
%
% See also: PLOT_OVERLAY_2D, PLOT_CORONAL_SLICES

    arguments
        overlay_image (:,:,:)
        bg_image (:,:,:)
        transducer_bowl (:,:,:)
        parameters 
        slice cell
        trans_pos (:,3)
        focus_pos (:,3)
        max_data_pos (1,3)
        options.show_rectangles = 1
        options.grid_step = parameters.grid.resolution_mm
        options.rect_size = 2
        options.overlay_threshold_low (1,1) = min(overlay_image(:))
        options.overlay_threshold_high (1,1) = min(overlay_image(:)) + ...
            (max(overlay_image(:)) - min(overlay_image(:))) * 0.8
        options.overlay_color_range = []
        options.use_overlay_alpha (1,1) = 1
        options.color_scale = 'viridis'
        options.rotation = 0
        options.show_colorbar = 1
        options.ticks = []
        options.tick_labels = []
        options.segmented_img = []
        options.bg_range = []
        options.overlay_segmented = 0  
    end

    %% Validate input positions against image boundaries
    if any(focus_pos > size(bg_image))
        error('Focus point is outside of image boundaries');
    end
    if ~isempty(trans_pos) && any(trans_pos > size(bg_image))
        error('Transducer point is outside of image boundaries');
    end
    if any(max_data_pos > size(bg_image))
        warning('Max intensity point is outside of image boundaries. Regularizing to max. of dimensions');
        max_data_pos = min([max_data_pos; size(bg_image)]);
    end

    %% Set thresholds and color range for intensity map
    if options.overlay_threshold_low == options.overlay_threshold_high
        options.overlay_threshold_low = options.overlay_threshold_low - 0.05;
    end
    if isempty(options.overlay_color_range)
        options.overlay_color_range = [max([options.overlay_threshold_low, min(overlay_image(:)), 0]), max(overlay_image(:))];
    end


    %% Extract transducer info

    tr = parameters.transducer(1);

    %% Extract specified slice from images and masks
    % Determine slicing axis and adjust positions accordingly
    slice_x = 1:size(overlay_image,1);
    slice_y = 1:size(overlay_image,2);
    slice_z = 1:size(overlay_image,3);
   
    if slice{1} == 'x'
        slice_x = slice{2};
        if ~isempty(focus_pos)
            focus_pos = focus_pos(2:3);
        end
        if ~isempty(trans_pos)
            trans_pos = trans_pos(2:3);
            natural_focus = trans_pos + [0, tr.(tr.type).curv_radius_mm / parameters.grid.resolution_mm];
        end
        max_data_pos = max_data_pos(2:3);
    elseif slice{1} == 'y'
        slice_y = slice{2};
        if ~isempty(focus_pos)
            focus_pos = focus_pos([1,3]);
        end
        if ~isempty(trans_pos)
            trans_pos = trans_pos([1,3]);
            natural_focus = trans_pos + [0, tr.(tr.type).curv_radius_mm / parameters.grid.resolution_mm];
        end
        max_data_pos = max_data_pos([1,3]);
    elseif slice{1} == 'z'
        slice_z = slice{2};
        if ~isempty(focus_pos)
            focus_pos = focus_pos(1:2);
        end
        if ~isempty(trans_pos)
            trans_pos = trans_pos([1,2]);
            natural_focus = trans_pos + [0, 0];
        end
        max_data_pos = max_data_pos([1,2]);
    else
        error("slice must be a cell array with the first element of 'x','y', or 'z' and a second element a slice number")
    end
    
    bg_slice = squeeze(bg_image(slice_x, slice_y, slice_z));
    if ~isempty(options.segmented_img)
        segmented_slice = squeeze(options.segmented_img(slice_x, slice_y, slice_z));
    end
    
    if (~isempty(options.segmented_img) && all(options.segmented_img == bg_image,'all')) || (isinteger(bg_image) && max(bg_image(:))< 12)
        bg_slice = mat2gray(bg_slice);
        bg_min =  min(bg_slice(:));
        bg_max =  max(bg_slice(:));

    else
        if ~isempty(options.bg_range)
            bg_min = options.bg_range(1);
            bg_max = options.bg_range(2);
        elseif ~isempty(options.segmented_img)
            %bg_slice(segmented_img_slice==0) = 0;
            bg_slice_brain = bg_slice(segmented_slice>0&segmented_slice<3);
            bg_min =  min(bg_slice_brain(:));
            bg_max =  max(bg_slice_brain(:));
            %bg_slice(bg_slice>0) = rescale(bg_slice(bg_slice>0), 'InputMin', bg_min, 'InputMax', bg_max);
        else
            bg_min = min(bg_slice(:));
            bg_max = max(bg_slice(:));
        end

        bg_slice = mat2gray(bg_slice, double([bg_min, bg_max]));
    end
    if ~isempty(trans_pos)
        transducer_bowl = mat2gray(squeeze(transducer_bowl(slice_x, slice_y, slice_z)));
    end
    overlay_image = gather(squeeze(overlay_image(slice_x, slice_y, slice_z)));

    if options.rotation
        R = [cosd(options.rotation) -sind(options.rotation); sind(options.rotation) cosd(options.rotation)];
        if ~isempty(focus_pos)
            focus_pos = round(R*double(focus_pos'));
        end
        if ~isempty(trans_pos)
            trans_pos = round(R*double(trans_pos'));
			natural_focus = round(R*double(natural_focus'));
        end
        max_data_pos = round(R*double(max_data_pos'));
        overlay_image = imrotate(overlay_image, options.rotation);
        bg_slice = imrotate(bg_slice, options.rotation);
        if options.overlay_segmented
            segmented_slice = imrotate(segmented_slice, options.rotation);
        end
        if ~isempty(trans_pos)
            transducer_bowl = imrotate(transducer_bowl, options.rotation);
        end
        
    end
        
    h = figure;
    ax1 = axes; % the background layer
    bg_img = cat(3, bg_slice, bg_slice, bg_slice);
    overlay_color = [0,0.2,0.7];
    bg_ones = ones(size(bg_img));
    overlay_weight = 0.2;
    overlay_img  = reshape(overlay_color, 1, 1,3).*bg_ones;
    image(bg_img+overlay_weight*overlay_img);
    axis image;
    axis off;
    
    % draw transducer (if grid was properly set up)
    if ~isempty(trans_pos) && parameters.modules.run_grid_setup == 1
        options.grid_step = parameters.grid.resolution_mm;
        plot_transducer_overlay(parameters, trans_pos, focus_pos, natural_focus, max_data_pos, options,  0.3, [0.2 0.6 1])
    end
    
    if options.overlay_segmented
        ax3 = axes;
        imagesc(ax3, segmented_slice,'alphadata', 0.3);
        axis image;
        axis off;
    end
    
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

    ax2.Visible = 'off';
    axis image
    axis off;
    
    if options.show_rectangles 
        rect_size = options.rect_size;
        rectangle('Position', [focus_pos(2)-rect_size/2 focus_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'r', 'LineWidth',1,'LineStyle','-')
        rectangle('Position', [max_data_pos(2)-rect_size/2 max_data_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'b', 'LineWidth',1,'LineStyle','-')
        text(max_data_pos(2), max_data_pos(1)+10, [num2str(round(max_data_pos(2)*parameters.grid.resolution_mm)), 'mm'], 'Color', 'w')
    end

    linkaxes([ax1, ax2], 'xy'); % Synchronize both axes to avoid shifts
    
    ax2.TightInset;
    ax1.TightInset;
    if options.show_colorbar
       C = colorbar(ax2);
       C.Position(4) = ax1.Position(4)*0.935;
       C.Position(2) = ax1.Position(2)+0.025;
       C.Position(1) = ax1.Position(1)+ax1.Position(3)-0.03;
       C.Color = 'white';
       if ~isempty(options.ticks)
           C.Ticks = options.ticks;
       end
       if ~isempty(options.tick_labels)
           C.TickLabels = options.tick_labels;
       end
           
    end
    ax2.Position = ax1.Position;
    
    child = ax1.Children(end);

    if isa(child, 'matlab.graphics.primitive.Image')
        ax2_colour = child.CData(1,1,:);
    elseif isa(child, 'matlab.graphics.chart.primitive.Line')
        ax2_colour = child.Color;
    else
        error('Unknown child type');
    end
    
    if any(ax2_colour > 1)
        ax2_colour = zeros(size(ax2_colour));
    end
    set(gcf,'color',ax2_colour);
end
