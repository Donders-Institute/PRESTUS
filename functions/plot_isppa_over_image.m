function [bg_slice, transducer_bowl, Isppa_map, ax1, ax2, bg_min, bg_max, h] = ...
    plot_isppa_over_image(Isppa_map, bg_image, transducer_bowl, parameters, slice, trans_pos, focus_pos, max_isppa_pos, options)

% PLOT_ISPPA_OVER_IMAGE Visualizes ISppa map overlaid on a 2D slice of a 3D background image.
%
% This function overlays the spatial peak pulse-average intensity (ISppa) map 
% on a specific 2D slice of a 3D background image (`bg_image`). It highlights 
% key positions such as the transducer position, focus position, and maximum ISppa 
% position. The function supports various customization options for visualization.
%
% Input:
%   Isppa_map       - [Nx x Ny x Nz] matrix representing the ISppa intensity map.
%   bg_image        - [Nx x Ny x Nz] matrix representing the 3D background image (e.g., anatomical image).
%   transducer_bowl - [Nx x Ny x Nz] binary mask representing the transducer bowl region.
%   parameters      - Struct containing simulation parameters (e.g., grid step size).
%   slice           - Cell array specifying the slice to visualize:
%                     * First element: axis ('x', 'y', or 'z').
%                     * Second element: slice number along the specified axis.
%   trans_pos       - [1x3] array specifying the transducer position in grid coordinates (row, col, slice).
%   focus_pos       - [1x3] array specifying the focus position in grid coordinates (row, col, slice).
%   max_isppa_pos   - [1x3] array specifying the maximum ISppa position in grid coordinates (row, col, slice).
%   options         - Struct containing optional visualization settings:
%                     * show_rectangles: Boolean flag to show rectangles for key positions (default: 1).
%                     * grid_step: Grid step size in mm (default: from parameters).
%                     * rect_size: Size of rectangles for key positions (default: 2).
%                     * isppa_threshold_low/high: Thresholds for alpha scaling of ISppa map.
%                     * isppa_color_range: Range for ISppa map color scaling.
%                     * color_scale: Colormap for ISppa map (default: 'viridis').
%                     * rotation: Rotation angle for visualization (default: 0).
%                     * show_colorbar: Boolean flag to display colorbar (default: 1).
%
% Output:
%   bg_slice        - Extracted 2D slice from the background image.
%   transducer_bowl - Extracted 2D slice from the transducer bowl mask.
%   Isppa_map       - Extracted 2D slice from the ISppa map.
%   ax1             - Handle to the axes displaying the background image.
%   ax2             - Handle to the axes displaying the ISppa map overlay.
%   bg_min          - Minimum intensity value of the background image slice.
%   bg_max          - Maximum intensity value of the background image slice.
%   h               - Handle to the created figure.

    arguments
        Isppa_map (:,:,:)
        bg_image (:,:,:)
        transducer_bowl (:,:,:)
        parameters 
        slice cell
        trans_pos (:,3)
        focus_pos (:,3)
        max_isppa_pos (1,3)
        options.show_rectangles = 1
        options.grid_step = parameters.grid_step_mm
        options.rect_size = 2
        options.isppa_threshold_low (1,1) = min(Isppa_map(:))
        options.isppa_threshold_high (1,1) = min(Isppa_map(:)) + ...
            (max(Isppa_map(:)) - min(Isppa_map(:))) * 0.8
        options.isppa_color_range = []
        options.use_isppa_alpha (1,1) = 1
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
    if any(max_isppa_pos > size(bg_image))
        error('Max ISPPA point is outside of image boundaries');
    end

    %% Set thresholds and color range for ISppa map
    if options.isppa_threshold_low == options.isppa_threshold_high
        options.isppa_threshold_low = options.isppa_threshold_low - 0.05;
    end
    if isempty(options.isppa_color_range)
        options.isppa_color_range = [max([options.isppa_threshold_low, min(Isppa_map(:)), 0]), max(Isppa_map(:))];
    end

    %% Extract specified slice from images and masks
    % Determine slicing axis and adjust positions accordingly
    slice_x = 1:size(Isppa_map,1);
    slice_y = 1:size(Isppa_map,2);
    slice_z = 1:size(Isppa_map,3);
   
    if strcmp(slice{1}, 'x')
        slice_x = slice{2};
        focus_pos = focus_pos(2:3);
        trans_pos = trans_pos(2:3);
        max_isppa_pos = max_isppa_pos(2:3);
    elseif strcmp(slice{1}, 'y')
        slice_y = slice{2};
        focus_pos = focus_pos([1,3]);
        trans_pos = trans_pos([1,3]);
        max_isppa_pos = max_isppa_pos([1,3]);
    elseif strcmp(slice{1}, 'z')
        slice_z = slice{2};
        focus_pos = focus_pos(1:2);
        trans_pos = trans_pos([1,2]);
        max_isppa_pos = max_isppa_pos([1,2]);
    else
        error("Slice must be a cell array with 'x', 'y', or 'z' as first element and a number as second element");
    end

    % Extract slices from input images and masks
    bg_slice = squeeze(bg_image(slice_x,slice_y,slice_z));
    
    % Normalize background image intensity values
    if isempty(options.bg_range)
        bg_min = min(bg_slice(:));
        bg_max = max(bg_slice(:));
    else
        bg_min = options.bg_range(1);
        bg_max = options.bg_range(2);
    end

    bg_slice = mat2gray(bg_slice,[bg_min,bg_max]);

    % Extract ISppa map and apply rotation if specified
    Isppa_map = squeeze(Isppa_map(slice_x,slice_y,slice_z));

    if options.rotation
        R = [cosd(options.rotation) -sind(options.rotation); sind(options.rotation) cosd(options.rotation)];
        if ~isempty(focus_pos)
        focus_pos = round(R*double(focus_pos'));
        end
        if ~isempty(trans_pos)
            trans_pos = round(R*double(trans_pos'));
        end
        max_isppa_pos = round(R*double(max_isppa_pos'));
        Isppa_map = imrotate(Isppa_map, options.rotation);
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
    
    % draw transducer
    if ~isempty(trans_pos)
        focal_slope = (trans_pos-focus_pos)/norm(trans_pos-focus_pos);
        focal_angle = atan2(focal_slope(2),focal_slope(1));
        
        grid_step = options.grid_step;
        ex_plane_pos = trans_pos - (parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm*[cos(focal_angle); sin(focal_angle)];
        geom_focus_pos = trans_pos - (parameters.transducer.curv_radius_mm)/grid_step*[cos(focal_angle), sin(focal_angle)];
        max_od = max(parameters.transducer.Elements_OD_mm);
        dist_to_ep = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od^2)/grid_step;
        ex_plane_pos_trig = geom_focus_pos + dist_to_ep *[cos(focal_angle), sin(focal_angle)];
        ort_angle = atan(-focal_slope(1)/focal_slope(2));

        arc_halfangle = atan(max_od/2/dist_to_ep/grid_step);

        r = max(parameters.transducer.Elements_OD_mm)/2/grid_step;

        hold on
        lineWidth = 1;

        boxColor = [235, 185, 47]/255*(1-overlay_weight) + overlay_color*overlay_weight;
        LineSmoothing = 'on';

        trans_full_depth = 16/grid_step;
        trans_back = ex_plane_pos_trig+trans_full_depth*focal_slope;
        line([trans_back(2)-r*sin(ort_angle), trans_back(2) + r*sin(ort_angle)], [trans_back(1)-r*cos(ort_angle), trans_back(1) + r*cos(ort_angle)],  'LineWidth', lineWidth, 'Color', boxColor,'LineSmoothing',LineSmoothing )
        line([ex_plane_pos_trig(2), trans_back(2)]-r*sin(ort_angle), [ex_plane_pos_trig(1), trans_back(1)]-r*cos(ort_angle),  'LineWidth', lineWidth,'Color', boxColor,'LineSmoothing',LineSmoothing  )
        line([ex_plane_pos_trig(2), trans_back(2)]+r*sin(ort_angle), [ex_plane_pos_trig(1), trans_back(1)]+r*cos(ort_angle),  'LineWidth', lineWidth,'Color', boxColor,'LineSmoothing',LineSmoothing )
        % exit plane
        line([ex_plane_pos_trig(2), ex_plane_pos_trig(2)]+r*sin(ort_angle), [trans_back(1)-r*cos(ort_angle), trans_back(1) + r*cos(ort_angle)],  'LineWidth', lineWidth,'Color', boxColor,'LineStyle', ':', 'LineSmoothing',LineSmoothing )
        % transducer curvature
        [arc_x, arc_y] = getArc(geom_focus_pos, parameters.transducer.curv_radius_mm/grid_step, focal_angle-arc_halfangle, focal_angle+arc_halfangle );
        plot(arc_y, arc_x, 'Color',boxColor,'LineWidth', lineWidth,'LineSmoothing',LineSmoothing )
    end
    
    if options.overlay_segmented
        ax3 = axes;
        imagesc(ax3, segmented_slice,'alphadata', 0.3);
        axis image;
        axis off;
    end
    
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

    ax2.Visible = 'off';
    axis image
    axis off;
    
    if options.show_rectangles 
        rect_size = options.rect_size;
        if ~isempty(trans_pos)
            rectangle('Position', [trans_pos(2)-1-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', boxColor, 'LineWidth',2,'LineStyle','-')
            text(trans_pos(2)-1, trans_pos(1)+10, [num2str(round((trans_pos(2)-1)*parameters.grid_step_mm))], 'Color', 'w')
            % the following plots the onset of the grid; note: grid is in voxels, not mm
            rectangle('Position', [trans_pos(2)-parameters.transducer.pos_grid(3)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', 'w', 'LineWidth',1,'LineStyle',':')
            rectangle('Position', [parameters.grid_dims(3)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', 'w', 'LineWidth',1,'LineStyle',':')
            % exit plane
            dist_to_exit_plane = parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm;
            % convert to voxels
            dist_to_exit_plane_vox = round(dist_to_exit_plane*(1/parameters.grid_step_mm));
            rectangle('Position', [(trans_pos(2)-1+dist_to_exit_plane_vox)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', boxColor, 'LineWidth',1,'LineStyle',':')
            %rectangle('Position', [(ex_plane_pos(4)-1)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', boxColor, 'LineWidth',1,'LineStyle',':')
            text(trans_pos(2)-1+dist_to_exit_plane_vox, trans_pos(1)+10, [num2str(round((trans_pos(2)-1+dist_to_exit_plane_vox)*parameters.grid_step_mm)), 'mm'], 'Color', 'w')
        end
        rectangle('Position', [focus_pos(2)-rect_size/2 focus_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'r', 'LineWidth',1,'LineStyle','-')
        rectangle('Position', [max_isppa_pos(2)-rect_size/2 max_isppa_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'b', 'LineWidth',1,'LineStyle','-')
        text(max_isppa_pos(2), max_isppa_pos(1)+10, [num2str(round(max_isppa_pos(2)*parameters.grid_step_mm)), 'mm'], 'Color', 'w')
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
    ax2_colour = ax1.Children(length(ax1.Children)).CData(1,1,:);
    if any(ax2_colour > 1)
        ax2_colour = zeros(size(ax2_colour));
    end
    set(gcf,'color',ax2_colour);
end
