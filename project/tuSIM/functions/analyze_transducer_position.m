function [dist_to_target, dist_to_target_mm, prop_intersect, mean_dts, var_dts] = analyze_transducer_position(trans_pos, pos_shift_mm, target, pixel_size, parameters, coord_mesh, full_skull_mask, skin_boundary_coords)

max_od = max(parameters.transducer.Elements_OD_mm);
dist_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od^2);

norm_v = (trans_pos-target)/norm(target-trans_pos);

trans_pos = trans_pos + norm_v*pos_shift_mm/pixel_size;

% dist_to_target = norm(target-trans_pos);
% dist_to_target_mm = dist_to_target * pixel_size;

geom_focus_pos = trans_pos - norm_v*(parameters.transducer.curv_radius_mm)/pixel_size;
ex_plane_pos = geom_focus_pos+norm_v*dist_to_ep_mm/pixel_size;

d = sum(norm_v.*ex_plane_pos);

orth_plane_disk = find(abs(coord_mesh_gpu.x*norm_v(1)+coord_mesh_gpu.y*norm_v(2)+coord_mesh_gpu.z*norm_v(3)-d)<0.5);

orth_plane_disk = orth_plane_disk(sqrt((coord_mesh_gpu.x(orth_plane_disk)-ex_plane_pos(1)).^2+...
    (coord_mesh_gpu.y(orth_plane_disk)-ex_plane_pos(2)).^2+...
    (coord_mesh_gpu.z(orth_plane_disk)-ex_plane_pos(3)).^2) < max_od/2);


prop_intersect = sum(full_skull_mask(orth_plane_disk))/length(orth_plane_disk);
% non_intersecting_vox = orth_plane_disk;
% non_intersecting_vox(full_skull_mask) = 0;

non_intersecting_vox_idx = intersect(orth_plane_disk,find(~full_skull_mask));
non_intersecting_vox_coords = [coord_mesh_gpu.x(non_intersecting_vox_idx), coord_mesh_gpu.y(non_intersecting_vox_idx), coord_mesh_gpu.z(non_intersecting_vox_idx)];
skin_boundary_coords = skin_boundary_coords(pdist2(skin_boundary_coords, ex_plane_pos)<40,:);
dist_to_skin = pdist2(non_intersecting_vox_coords, skin_boundary_coords);
dist_to_skin = min(dist_to_skin, [], 2);
mean_dts = mean(dist_to_skin(:));
var_dts = var(dist_to_skin(:));
dist_to_skin_3d_vol = [];
% dist_to_skin_3d_vol = double(orth_plane_disk);
% dist_to_skin_3d_vol(non_intersecting_vox_idx) = dist_to_skin;
end