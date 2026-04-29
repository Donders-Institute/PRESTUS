function tp_plot_geometry_overlay(img, target_pos, trans_pos, pixel_size, parameters, subject_id, target_name, outer_sphere_3d)
% TP_PLOT_GEOMETRY_OVERLAY  Visualise transducer geometry on a skin segmentation slice
%
% Overlays transducer position, geometric focus, exit plane, and ray path
% on the central Y-slice through the target for positioning validation.
% Saves the resulting figure to disk.
%
% Use as:
%   tp_plot_geometry_overlay(img, target_pos, trans_pos, pixel_size, ...
%       parameters, subject_id, target_name, outer_sphere_3d)
%
% Input:
%   img             - [Nx x Ny x Nz] tissue segmentation volume
%   target_pos      - [1x3] target position in voxel space
%   trans_pos       - [1x3] transducer position in voxel space
%   pixel_size      - voxel size [mm]
%   parameters      - (1,1) simulation parameters struct
%   subject_id      - subject identifier used in output filename
%   target_name     - target label string used in output filename
%   outer_sphere_3d - [Nx x Ny x Nz] logical mask of skull surface positions
%
% See also: TP_PLOT_CANDIDATE_POSITIONS, TP_FIND_INITIAL_CANDIDATE

disp("[TP] Generating and visualizing geometric overlay ...")

sz = size(img);  % Get image dimensions [Nx Ny Nz] - overrides input sz
[t1_x, t1_y, t1_z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));  % 3D voxel coordinate grids (X,Y,Z)

% Transducer geometry parameters (physical → voxel space)
max_od_mm = max(parameters.transducer(1).(parameters.transducer(1).type).elem_od_mm); % Largest element diameter (mm) - aperture size
% Sagitta calculation: distance from geometric focus to exit plane
% h = R - sqrt(R^2 - (D/2)^2) where R=curvature radius, D=aperture diameter
% Exit plane lies at h/2 from sphere center along optical axis
dist_gf_to_ep_mm = 0.5 * sqrt(4*parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm^2 - max_od_mm^2);

% Unit normal vector: direction from transducer → target (propagation axis)
norm_v = (trans_pos - target_pos) / norm(target_pos - trans_pos); 

% Geometric focus: sphere center (curvature radius along normal from trans_pos)
geom_focus_pos = trans_pos - norm_v * (parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm / pixel_size);

% Exit plane: halfway between geometric focus and aperture plane
ex_plane_pos = geom_focus_pos + norm_v * (dist_gf_to_ep_mm / pixel_size);

% Aperture disk radius in voxel units
max_od_grid = max_od_mm / pixel_size;

% Define exit plane disk: orthogonal to propagation axis, centered at ex_plane_pos
d = sum(norm_v .* ex_plane_pos);  % Plane equation: n·x = d
orth_plane_disk = abs(t1_x*norm_v(1) + t1_y*norm_v(2) + t1_z*norm_v(3) - d) < 0.5 & ...
                  sqrt((t1_x-ex_plane_pos(1)).^2 + (t1_y-ex_plane_pos(2)).^2 + (t1_z-ex_plane_pos(3)).^2) < max_od_grid/2;

% Create visualization mask: skin(5)=skin, disk=2, sphere_intersection=3
skin_only = uint8(img == 5);  % Extract skin surface only
skin_only(orth_plane_disk) = 2;  % Mark exit plane aperture
skin_only(find(outer_sphere_3d)) = 3;  % Mark skull intersection sphere (from tp_find_candidate_positions)

%% Plot central Y-slice through target

h = figure; 
imagesc(squeeze(skin_only(:, target_pos(2), :))); hold on;  % X-Z slice at target Y
colormap(gray);  % Monochrome for segmentation clarity

% Extract X,Z coordinates for 2D overlay (flip for imagesc convention)
trans_xz = trans_pos([1,3]); target_xz = target_pos([1,3]);

% Draw position markers (4-voxel boxes)
rectangle('Position', [flip(trans_xz)-2, 4, 4], 'Curvature', [0,0], 'EdgeColor', 'b', 'LineWidth', 2);     % Transducer (blue)
rectangle('Position', [flip(geom_focus_pos([1,3]))-2, 4, 4], 'Curvature', [0,0], 'EdgeColor', 'yellow', 'LineWidth', 2);  % Geo focus (yellow)
rectangle('Position', [flip(ex_plane_pos([1,3]))-2, 4, 4], 'Curvature', [0,0], 'EdgeColor', 'white', 'LineWidth', 2);     % Exit plane (white)  
rectangle('Position', [flip(target_xz)-2, 4, 4], 'Curvature', [0,0], 'EdgeColor', 'r', 'LineWidth', 2);     % Target (red)

% Draw propagation path (white line)
line([trans_xz(2) target_xz(2)], [trans_xz(1) target_xz(1)], 'Color', 'white', 'LineWidth', 2);

% Add 3D transducer bowl visualization
get_transducer_box(trans_xz, target_xz, [], pixel_size, parameters);

% Save geometry validation plot
saveas(h, fullfile(parameters.io.figures_preproc_dir, sprintf('sub-%03d_geometry_%s.png', subject_id, target_name)), 'png');
close(h);

end
