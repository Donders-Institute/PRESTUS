function [bg_slice, transducer_bowl, Isppa_map] = plot_isppa_over_image(Isppa_map, bg_image, transducer_bowl, parameters, slice, trans_pos, focus_pos, max_isppa_pos, show_rectangles, grid_step, rect_size)
    arguments
        Isppa_map (:,:,:)
        bg_image (:,:,:)
        transducer_bowl (:,:,:)
        parameters 
        slice cell
        trans_pos (1,3)
        focus_pos (1,3)
        max_isppa_pos (1,3)
        show_rectangles = 1 % show rectangles for transducer, focus, and the real focus
        grid_step = parameters.grid_step_mm
        rect_size = 2

    end
    
    slice_x = 1:size(Isppa_map,1);
    slice_y = 1:size(Isppa_map,2);
    slice_z = 1:size(Isppa_map,3);

   
    if slice{1} == 'x'
        slice_x = slice{2};
        focus_pos = focus_pos(2:3);
        trans_pos = trans_pos(2:3);
        max_isppa_pos = max_isppa_pos(2:3);
    elseif slice{1} == 'y'
        slice_y = slice{2};
        focus_pos = focus_pos([1,3]);
        trans_pos = trans_pos([1,3]);
        max_isppa_pos = max_isppa_pos([1,3]);
    elseif slice{1} == 'z'
        slice_z = slice{2};
        focus_pos = focus_pos(1:2);
        trans_pos = trans_pos([1,2]);
        max_isppa_pos = max_isppa_pos([1,2]);
    else
        error("slice must be a cell array with the first element of 'x','y', or 'z' and a second element a slice number")
    end
    
    bg_slice = mat2gray(squeeze(bg_image(slice_x, slice_y, slice_z)));
    transducer_bowl = mat2gray(squeeze(transducer_bowl(slice_x, slice_y, slice_z)));
    %before_exit_plane_mask = ~squeeze(after_exit_plane_mask(slice_x, slice_y, slice_z));
    Isppa_map = squeeze(Isppa_map(slice_x, slice_y, slice_z));

    figure;
    ax1 = axes;
    bg_zeros = zeros(size(bg_slice));
    bg_slice = (bg_slice-min(bg_slice(:)))/ max(bg_slice(:));
    bg_img = cat(3, bg_slice, bg_slice, bg_slice);
    overlay_color = [0,0.2,0.7];
    bg_ones = ones(size(bg_img));
    overlay_weight = 0.2;
    overlay_img  = reshape(overlay_color, 1, 1,3).*bg_ones;
    imagesc(bg_img+overlay_weight*overlay_img);
    %colormap(ax1,'gray');
    axis image;
    axis off;

    focal_slope = (trans_pos-focus_pos)/norm(trans_pos-focus_pos);
    focal_angle = atan2(focal_slope(2),focal_slope(1));
    
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
    %plot(trans_pos(2),trans_pos(1), 'r+', 'MarkerSize', 30, 'LineWidth', 2);

    trans_full_depth = 16/grid_step;
    trans_back = ex_plane_pos_trig+trans_full_depth*focal_slope;
    %line([trans_pos(2), targ_xz(2)], [trans_pos(1), focus_pos(1)],  'LineWidth', lineWidth, 'Color',boxColor )
    line([trans_back(2)-r*sin(ort_angle), trans_back(2) + r*sin(ort_angle)], [trans_back(1)-r*cos(ort_angle), trans_back(1) + r*cos(ort_angle)],  'LineWidth', lineWidth, 'Color', boxColor,'LineSmoothing',LineSmoothing )
    %line([ex_plane_pos_trig(2)-r*sin(ort_angle), ex_plane_pos_trig(2) + r*sin(ort_angle)], [ex_plane_pos_trig(1)-r*cos(ort_angle), ex_plane_pos_trig(1) + r*cos(ort_angle)],  'LineWidth', lineWidth,'Color', 'red')
    line([ex_plane_pos_trig(2), trans_back(2)]-r*sin(ort_angle), [ex_plane_pos_trig(1), trans_back(1)]-r*cos(ort_angle),  'LineWidth', lineWidth,'Color', boxColor,'LineSmoothing',LineSmoothing  )
    line([ex_plane_pos_trig(2), trans_back(2)]+r*sin(ort_angle), [ex_plane_pos_trig(1), trans_back(1)]+r*cos(ort_angle),  'LineWidth', lineWidth,'Color', boxColor,'LineSmoothing',LineSmoothing )
    [arc_x, arc_y] = getArc(geom_focus_pos, parameters.transducer.curv_radius_mm/grid_step, focal_angle-arc_halfangle, focal_angle+arc_halfangle );
    %rectangle('Position', [geom_focus_pos(2)-rect_size/2 geom_focus_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'r', 'LineWidth',1,'LineStyle','-')

    plot(arc_y, arc_x, 'Color',boxColor,'LineWidth', lineWidth,'LineSmoothing',LineSmoothing )

    ax2 = axes;
    Isppa_map(Isppa_map<=0.5) = 0;
    isppa_alpha = (Isppa_map-min(Isppa_map))/max(Isppa_map(:));
    imagesc(ax2,Isppa_map,'alphadata', isppa_alpha);
    colormap(ax2,'viridis');
    %caxis(ax2,[min(Isppa_map(:)) max(Isppa_map(:))]);
    ax2.Visible = 'off';
    %linkprop([ax1 ax2],'Position');
    axis image
    axis off;
    if show_rectangles 
        rectangle('Position', [trans_pos(2)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', 'g', 'LineWidth',1,'LineStyle','-')
        rectangle('Position', [focus_pos(2)-rect_size/2 focus_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'r', 'LineWidth',1,'LineStyle','-')
        rectangle('Position', [max_isppa_pos(2)-rect_size/2 max_isppa_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'b', 'LineWidth',1,'LineStyle','-')
    end
    %colorbar;

end
