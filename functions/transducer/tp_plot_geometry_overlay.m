function tp_plot_geometry_overlay(img, target, trans_pos, pixel_size, parameters, subject_id, target_name, t1_x, t1_y, t1_z)
    max_od_mm = max(parameters.transducer.Elements_OD_mm); % Largest element diameter (mm) - defines transducer aperture
    % Sagitta/2: distance from geometric focus to exit plane
    % Derivation: h = R - sqrt(R^2-(D/2)^2) where R=curv_radius, D=max_od_mm
    % Exit plane = h/2 from sphere center = 0.5*sqrt(4R^2-D^2)
    dist_gf_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od_mm^2);
    norm_v = (trans_pos-target)/norm(target-trans_pos); % Unit normal vector pointing from transducer → target
    geom_focus_pos = trans_pos - norm_v*(parameters.transducer.curv_radius_mm)/pixel_size;
    ex_plane_pos = geom_focus_pos+norm_v*dist_gf_to_ep_mm/pixel_size;
    max_od_grid = max_od_mm/pixel_size;
    d = sum(norm_v.*ex_plane_pos);
    orth_plane_disk = abs(t1_x*norm_v(1)+t1_y*norm_v(2)+t1_z*norm_v(3)-d)<0.5 & ...
        sqrt((t1_x-ex_plane_pos(1)).^2+(t1_y-ex_plane_pos(2)).^2+(t1_z-ex_plane_pos(3)).^2) < max_od_grid/2;
    
    skin_only = uint8(img==5); skin_only(orth_plane_disk)=2; skin_only(find(outer_sphere_3d))=3;

    h = figure; imagesc(squeeze(skin_only(:,target(2),:))); hold on;
    trans_xz = trans_pos([1,3]); target_xz = target([1,3]);
    rectangle('Position',[flip(trans_xz)-2, 4, 4],'Curvature',[0,0],'EdgeColor','b','LineWidth',2);
    rectangle('Position',[flip(geom_focus_pos([1,3]))-2, 4, 4],'Curvature',[0,0],'EdgeColor','yellow','LineWidth',2);
    rectangle('Position',[flip(ex_plane_pos([1,3]))-2, 4, 4],'Curvature',[0,0],'EdgeColor','white','LineWidth',2);
    rectangle('Position',[flip(target_xz)-2, 4, 4],'Curvature',[0,0],'EdgeColor','r','LineWidth',2);
    line([trans_xz(2) target_xz(2)], [trans_xz(1) target_xz(1)], 'Color', 'white');
    get_transducer_box(trans_xz, target_xz, pixel_size, parameters);
    saveas(h, fullfile(parameters.output_dir,sprintf('sub-%03d_geometry_%s.png', subject_id, target_name)), 'png'); close(h);
end