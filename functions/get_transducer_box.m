function [transducer_box, ex_plane_pos_trig, geom_focus_pos, dist_to_ep_mm] = get_transducer_box(trans_pos, focus_pos, grid_step, parameters, plot)
    arguments
        trans_pos (1, 2)
        focus_pos (1, 2)
        grid_step (1,1) % grid step in mm
        parameters struct
        plot = 1
    end
    focal_slope = (trans_pos-focus_pos)/norm(trans_pos-focus_pos);
    focal_angle = atan2(focal_slope(2),focal_slope(1));

    geom_focus_pos = trans_pos - (parameters.transducer.curv_radius_mm)/grid_step*[cos(focal_angle), sin(focal_angle)];
    max_od = max(parameters.transducer.Elements_OD_mm);
    dist_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od^2);

    dist_to_ep_grid = dist_to_ep_mm/grid_step;
    ex_plane_pos_trig = geom_focus_pos + dist_to_ep_grid *[cos(focal_angle), sin(focal_angle)];
    ort_angle = atan(-focal_slope(1)/focal_slope(2));

    r = max(parameters.transducer.Elements_OD_mm)/2/grid_step;

    trans_full_depth = 16/grid_step;
    trans_back = ex_plane_pos_trig+trans_full_depth*focal_slope;

    transducer_box = [[trans_back(2)-r*sin(ort_angle), trans_back(1)-r*cos(ort_angle)],...
        [trans_back(2) + r*sin(ort_angle), trans_back(1) + r*cos(ort_angle)],...
        [ex_plane_pos_trig(2) - r*sin(ort_angle), ex_plane_pos_trig(1)-r*cos(ort_angle)],...
        [ex_plane_pos_trig(2) + r*sin(ort_angle), ex_plane_pos_trig(1) + r*cos(ort_angle)]];

    if plot 

        overlay_weight = 0;
        overlay_color = [0,0.2,0.7];
        lineWidth = 1;

        boxColor = [235, 185, 47]/255*(1-overlay_weight) + overlay_color*overlay_weight;
        LineSmoothing = 'on';

        line([trans_back(2)-r*sin(ort_angle), trans_back(2) + r*sin(ort_angle)], [trans_back(1)-r*cos(ort_angle), trans_back(1) + r*cos(ort_angle)],  'LineWidth', lineWidth, 'Color', boxColor,'LineSmoothing',LineSmoothing )

        line([ex_plane_pos_trig(2), trans_back(2)]-r*sin(ort_angle), [ex_plane_pos_trig(1), trans_back(1)]-r*cos(ort_angle),  'LineWidth', lineWidth,'Color', boxColor,'LineSmoothing',LineSmoothing  )
        line([ex_plane_pos_trig(2), trans_back(2)]+r*sin(ort_angle), [ex_plane_pos_trig(1), trans_back(1)]+r*cos(ort_angle),  'LineWidth', lineWidth,'Color', boxColor,'LineSmoothing',LineSmoothing )
    end
end