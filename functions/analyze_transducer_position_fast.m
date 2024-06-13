function [prop_intersect, mean_dist_skin, var_dist_skin, mean_dist_skull, var_dist_skull] = ...
    analyze_transducer_position_fast(i, norm_v, ex_plane_pos, coord_mesh_gpu,...
    full_skull_mask_idx, skin_boundary_coords, skull_boundary_coords, max_od)

% max_od = max(parameters.transducer.Elements_OD_mm);
% dist_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od^2);
% 
% norm_v = (trans_pos-target)/norm(target-trans_pos);
% 
% trans_pos = trans_pos + norm_v*pos_shift_mm/pixel_size;
% 
% % dist_to_target = norm(target-trans_pos);
% % dist_to_target_mm = dist_to_target * pixel_size;
% 
% geom_focus_pos = trans_pos - norm_v*(parameters.transducer.curv_radius_mm)/pixel_size;
% ex_plane_pos = geom_focus_pos+norm_v*dist_to_ep_mm/pixel_size;

d = sum(norm_v(i,:).*ex_plane_pos(i,:));
% 
orth_plane_disk = find(abs(coord_mesh_gpu.x*norm_v(i,1)+coord_mesh_gpu.y*norm_v(i,2)+coord_mesh_gpu.z*norm_v(i,3)-d)<0.5);

orth_plane_disk = orth_plane_disk(sqrt((coord_mesh_gpu.x(orth_plane_disk)-ex_plane_pos(i,1)).^2+...
    (coord_mesh_gpu.y(orth_plane_disk)-ex_plane_pos(i,2)).^2+...
    (coord_mesh_gpu.z(orth_plane_disk)-ex_plane_pos(i,3)).^2) < (max_od/2+4));
% 
% 
% orth_plane_disk = find(abs(sum(coord_mesh_gpu.xyz.*norm_v(i,:),2)-d)<0.5);
% 
% orth_plane_disk = orth_plane_disk(sqrt(sum((coord_mesh_gpu.xyz(orth_plane_disk,:)-ex_plane_pos(i,:)).^2,2)) < max_od/2);
% 
% 


prop_intersect = length(intersect(full_skull_mask_idx,orth_plane_disk))/length(orth_plane_disk);
% non_intersecting_vox = orth_plane_disk;
% non_intersecting_vox(full_skull_mask) = 0;

non_intersecting_vox_idx = setdiff(orth_plane_disk, full_skull_mask_idx);
if isempty(orth_plane_disk) ||  isempty(non_intersecting_vox_idx) || prop_intersect > 0.3
    if isempty(orth_plane_disk) ||  isempty(non_intersecting_vox_idx)
        prop_intersect = 1;
    end
    mean_dist_skin = gpuArray(0);
    var_dist_skin = gpuArray(0);
    mean_dist_skull = gpuArray(0);
    var_dist_skull = gpuArray(0);
    return;
end

non_intersecting_vox_coords = coord_mesh_gpu.xyz(non_intersecting_vox_idx,:);
skin_boundary_coords = skin_boundary_coords(pdist2(skin_boundary_coords, ex_plane_pos(i,:))<50,:);
skull_boundary_coords = skull_boundary_coords(pdist2(skull_boundary_coords, ex_plane_pos(i,:))<100,:);

dist_to_skin = pdist2(non_intersecting_vox_coords, skin_boundary_coords);
dist_to_skin = min(dist_to_skin, [], 2);

dist_to_skull = pdist2(coord_mesh_gpu.xyz(orth_plane_disk,:), skull_boundary_coords);
dist_to_skull = min(dist_to_skull, [], 2);

mean_dist_skin = mean(dist_to_skin(:));
var_dist_skin = var(dist_to_skin(:));

mean_dist_skull = mean(dist_to_skull(:));
var_dist_skull = var(dist_to_skull(:));
% dist_to_skin_3d_vol = [];
% dist_to_skin_3d_vol = double(orth_plane_disk);
% dist_to_skin_3d_vol(non_intersecting_vox_idx) = dist_to_skin;
end