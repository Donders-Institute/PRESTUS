function [prop_intersect, mean_dist_skin, var_dist_skin, mean_dist_skull, var_dist_skull] = ...
    analyze_transducer_position_fast(i, norm_v, ex_plane_pos, coord_mesh_gpu,...
    full_skull_mask_idx, skin_boundary_coords, skull_boundary_coords, max_od)

d = sum(norm_v(i,:).*ex_plane_pos(i,:));

orth_plane_disk = find(abs(coord_mesh_gpu.x*norm_v(i,1)+coord_mesh_gpu.y*norm_v(i,2)+coord_mesh_gpu.z*norm_v(i,3)-d)<0.5);
orth_plane_disk = orth_plane_disk(sqrt((coord_mesh_gpu.x(orth_plane_disk)-ex_plane_pos(i,1)).^2+...
    (coord_mesh_gpu.y(orth_plane_disk)-ex_plane_pos(i,2)).^2+...
    (coord_mesh_gpu.z(orth_plane_disk)-ex_plane_pos(i,3)).^2) < (max_od/2+4));

prop_intersect = length(intersect(full_skull_mask_idx,orth_plane_disk))/length(orth_plane_disk);

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

end