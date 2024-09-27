function show_3d_head(segmented_img, target_xyz, trans_xyz, parameters, pixel_size, coord_mesh_xyz, crop_at_target, view_angle, open_figure)
arguments
    segmented_img (:,:,:) 
    target_xyz (1,3)
    trans_xyz (1,3)
    parameters struct
    pixel_size (1,1)
    coord_mesh_xyz (:,3) 
    crop_at_target (1,3) = [0,0,0]
    view_angle (1,2) = [0,0]
    open_figure (1,1) = 1
end

% open figure
if open_figure
    figure
end
colormap(gray(80))

%% plot skin surface

segmented_img_to_plot = segmented_img(1:2:end,1:2:end,1:2:end);
% crop if desired
origin_shift = [0 0 0];
if any(crop_at_target)
    original_size = size(segmented_img_to_plot);
    for i = 1:3
        cropping_indices = {':',':',':'};
        if crop_at_target(i) == 1
            cropping_indices(i) = {floor(target_xyz(i)/2):size(segmented_img_to_plot,i)};
        elseif crop_at_target(i) == -1
            cropping_indices(i) = {1:(floor(target_xyz(i)/2)+1)};
            origin_shift(i) = floor(target_xyz(i)/2);
        else
            continue;
        end
        S = struct;
        S.type = '()';
        S.subs = cropping_indices;
        segmented_img_to_plot = subsasgn(segmented_img_to_plot, S, []);
    end
    origin_shift = origin_shift([2,1,3]);
end
if gpuDeviceCount > 0 
    Ds = smooth3(gpuArray(double(segmented_img_to_plot>0)));
else 
    Ds = smooth3(double(segmented_img_to_plot>0));
end
skin_isosurface = isosurface(Ds,0.5);
skin_isosurface.vertices = skin_isosurface.vertices - origin_shift;
hiso = patch(skin_isosurface,...
   'FaceColor',[1,.75,.65],...
   'EdgeColor','none',...
   'facealpha',0.8);
isonormals(Ds,hiso);
hiso.SpecularColorReflectance = 0;
hiso.SpecularExponent = 50;

%% plot GM, WM, CSF, bone (?)

for_caps = segmented_img_to_plot;
for_caps(for_caps>4) = 0; % only GM, WM, CSF, bone
end_slices = isocaps(for_caps*10,4);
end_slices.vertices = end_slices.vertices - origin_shift;
hcap = patch(end_slices,...
   'FaceColor','interp',...
   'EdgeColor','none',...
   'facealpha',0.8);
lighting gouraud
hcap.AmbientStrength = 0.6;

%% plot transducer slice at exit plane
% do not plot a transducer if view is cropped

if ~any(crop_at_target)
    max_od_mm = max(parameters.transducer.Elements_OD_mm);
    max_od_grid = max_od_mm / pixel_size;
    norm_vector = (trans_xyz-target_xyz)/norm(trans_xyz-target_xyz);
    % location of geometric focus
    geom_focus_xyz = trans_xyz - norm_vector*(parameters.transducer.curv_radius_mm)/pixel_size;
    % distance of geometric focus to exit plane
    dist_gf_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od_mm^2);
    % coordinates of exit plane
    ex_plane_xyz = geom_focus_xyz+norm_vector*dist_gf_to_ep_mm/pixel_size;
    % calculate a distance metric
    d = sum(norm_vector.*ex_plane_xyz);
    % identify nearest plane to exit plane
    orth_plane_disk = find(abs(sum(coord_mesh_xyz.*norm_vector,2)-d)<0.5);
    orth_plane_disk = orth_plane_disk(sqrt(sum((coord_mesh_xyz(orth_plane_disk,:)-ex_plane_xyz).^2,2)) < max_od_grid/2);
    orth_plane_3d = zeros(size(segmented_img));
    orth_plane_3d(orth_plane_disk) = 1;
    orth_plane_3d = smooth3(orth_plane_3d);
    orth_plane_3d_surf = isosurface(orth_plane_3d(1:2:end,1:2:end,1:2:end));
    orth_plane_3d_surf.vertices = orth_plane_3d_surf.vertices - origin_shift;
    patch(orth_plane_3d_surf,...
       'FaceColor','blue',...
       'EdgeColor','none');
end

%% plot target sphere

sphere_3d = zeros(size(segmented_img));
sphere_3d(pdist2(coord_mesh_xyz, target_xyz)<3) = 1;
sphere_surf = isosurface(smooth3(sphere_3d(1:2:end,1:2:end,1:2:end)));
sphere_surf.vertices = sphere_surf.vertices - origin_shift;
patch(sphere_surf,'FaceColor','red','EdgeColor','red','facealpha',0.8)

%% adjust viewing angles

if ~isempty(view_angle)
    view(view_angle)
end
xlabel('X')
ylabel('Y')
set(gca, 'YDir','reverse')
zlabel('Z')

if any(crop_at_target)
    original_size = original_size([2,1,3]);
    if any(origin_shift)
        axis(reshape([-origin_shift ; -origin_shift + original_size ], [1 6]));
    end
end

if view_angle(1) < 0
    camlight('right')
else
    camlight('left')
end

end