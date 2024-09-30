function [bg_slice, transducer_bowl, Isppa_map, ax1, ax2, bg_min, bg_max, h] = ...
    plot_isppa_over_image(Isppa_map, bg_image, transducer_bowl, parameters, slice, trans_pos, focus_pos, max_isppa_pos, options )
    arguments
        Isppa_map (:,:,:)
        bg_image (:,:,:)
        transducer_bowl (:,:,:)
        parameters 
        slice cell
        trans_pos (:,3)
        focus_pos (:,3)
        max_isppa_pos (1,3)
        options.show_rectangles = 1 % show rectangles for transducer, focus, and the real focus
        options.grid_step = parameters.grid_step_mm
        options.rect_size = 2
        options.isppa_threshold_low (1,1) = min(Isppa_map, [], 'all') % threshold for the lower boundary of the alpha values (i.e., alpha would scale from 0 to 1 on the range [isppa_threshold_low, isppa_threshold_high])
        options.isppa_threshold_high (1,1) = min(Isppa_map, [], 'all') + (max(Isppa_map, [], 'all')-min(Isppa_map, [], 'all'))*0.8
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
    
    if any(focus_pos > size(bg_image))
        error('Focus point is outside of image boundaries')
    end
    if ~isempty(trans_pos)
        if any(trans_pos>size(bg_image))
            error('Transducer point is outside of image boundaries')
        end
    end
    if any(max_isppa_pos>size(bg_image))
        error('Max ISPPA point is outside of image boundaries')
    end
    
    if options.isppa_threshold_low == options.isppa_threshold_high
        options.isppa_threshold_low = options.isppa_threshold_low - 0.05;
    end
    if isempty(options.isppa_color_range)
        options.isppa_color_range = [max([options.isppa_threshold_low, min(Isppa_map(:)), 0]), max(Isppa_map(:))];
    end
    
    slice_x = 1:size(Isppa_map,1);
    slice_y = 1:size(Isppa_map,2);
    slice_z = 1:size(Isppa_map,3);
   
    if slice{1} == 'x'
        slice_x = slice{2};
        if ~isempty(focus_pos)
            focus_pos = focus_pos(2:3);
        end
        if ~isempty(trans_pos)
            trans_pos = trans_pos(2:3);
        end
        max_isppa_pos = max_isppa_pos(2:3);
    elseif slice{1} == 'y'
        slice_y = slice{2};
        if ~isempty(focus_pos)
            focus_pos = focus_pos([1,3]);
        end
        if ~isempty(trans_pos)
            trans_pos = trans_pos([1,3]);
        end
        max_isppa_pos = max_isppa_pos([1,3]);
    elseif slice{1} == 'z'
        slice_z = slice{2};
        if ~isempty(focus_pos)
        focus_pos = focus_pos(1:2);
        end
        if ~isempty(trans_pos)
        trans_pos = trans_pos([1,2]);
        end
        max_isppa_pos = max_isppa_pos([1,2]);
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
    %before_exit_plane_mask = ~squeeze(after_exit_plane_mask(slice_x, slice_y, slice_z));
    Isppa_map = gather(squeeze(Isppa_map(slice_x, slice_y, slice_z)));
    %Isppa_map(Isppa_map<=0.5) = 0;

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
        %ex_plane_pos = trans_pos - (parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm*[cos(focal_angle); sin(focal_angle)];
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
            rectangle('Position', [trans_pos(2)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', 'g', 'LineWidth',1,'LineStyle','-')
        end
        rectangle('Position', [focus_pos(2)-rect_size/2 focus_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'r', 'LineWidth',1,'LineStyle','-')
        rectangle('Position', [max_isppa_pos(2)-rect_size/2 max_isppa_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'b', 'LineWidth',1,'LineStyle','-')
    end
    
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
