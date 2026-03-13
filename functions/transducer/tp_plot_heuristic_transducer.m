function tp_plot_heuristic_transducer(tpos_pars, img, img_info, parameters, pixel_size, target_name, subject_id)
% TP_PLOT_HEURISTIC_TRANSDUCER Create visualization of heuristic transducer placement
%
% INPUT (7 required + 1 optional)
%   tpos_pars   - Table from tp_evaluate_candidate_positions  
%   img         - Segmented head image
%   img_info    - NIfTI header info  
%   parameters  - Parameters struct
%   pixel_size  - Scalar voxel size in mm
%   target_name - String target identifier
%   subject_id  - Scalar subject ID
%   mesh        - [OPTIONAL LAST] Full mesh structure from tp_candidate_mesh()
%
% OUTPUT
%   Saves: sub-XXX_optimal_TARGET.png

%% Select optimal position (multi-criteria heuristic)
tppf = tpos_pars(tpos_pars.prop_intersect<0.05,:);
tppf = tppf(tppf.mean_dist_skull <= quantile(tppf.mean_dist_skull, 0.5) & ...
           tppf.var_dist_skull <= quantile(tppf.var_dist_skull, 0.1),:);
tppf = tppf(tppf.var_dist_skin==min(tppf.var_dist_skin),:);
i = 1;

target_xyz = [tppf.targ_x(i), tppf.targ_y(i), tppf.targ_z(i)];
trans_xyz = [tppf.trans_x(i), tppf.trans_y(i), tppf.trans_z(i)];

%% Recompute mesh variables

mesh = tp_candidate_mesh(img, target_xyz, parameters, pixel_size);

coord_mesh.xyz  = gather(mesh.coord_mesh.xyz);      % struct with field xyz
ex_plane        = gather(mesh.ex_plane);            % [N_cand x 3] exit plane centers
max_od_grid     = gather(mesh.max_od_grid);         % scalar, aperture diameter in voxels
geom_focus      = gather(mesh.geom_focus);          % [N_cand x 3] geometric focus positions

ex_pl_xyz = ex_plane(i,:);
gf_xyz = round(geom_focus(i,:));

%% Align to focal axis 

[rotated_img, trans_xyz, target_xyz, transformation_matrix, ~, ~, ~, ~] = ...
    preproc_align_to_focal_axis(...
    img, img_info, round(trans_xyz)', target_xyz', 1, parameters);

TF = maketform('affine', transformation_matrix);

ex_pl_xyz = round(tformfwd([ex_pl_xyz]', TF));
gf_xyz  = round(tformfwd([gf_xyz]', TF));

%% Visualization

h = figure('units', 'normalized', 'position', [0 0 1 .5]);
    subplot(1,3,1);
        colormap([0.3 0.3 0.3; lines(5)]);
        imagesc(squeeze(rotated_img(:,round(trans_xyz(2)),:)));
        
        % TARGET - Red box
        rectangle('Position', [target_xyz([3,1]) - 2, 4, 4], ...
                  'Curvature', [0,0], ...
                  'EdgeColor', 'r', ...
                  'LineWidth', 2, 'LineStyle', '-');
        
        % TRANSDUCER (heuristic position) - Blue box  
        rectangle('Position', [trans_xyz([3,1]) - 2, 4, 4],...
                  'Curvature', [0,0],...
                  'EdgeColor', 'b', ...
                  'LineWidth', 2, 'LineStyle', '-');
        
        % GEOMETRIC FOCUS (bowl sphere center) - Yellow box
        rectangle('Position', [gf_xyz([3,1]) - 2, 4, 4],...
                  'Curvature', [0,0],...
                  'EdgeColor', 'yellow',... 
                  'LineWidth', 2, 'LineStyle', '-');
        
        % EXIT PLANE (transducer aperture center) - White box
        rectangle('Position', [ex_pl_xyz([3,1]) - 2, 4, 4],...
                  'Curvature', [0,0],...
                  'EdgeColor', 'white', ...
                  'LineWidth', 2, 'LineStyle', '-');

        line([trans_xyz(3) target_xyz(3)], [trans_xyz(1) target_xyz(1)], 'Color', 'white');
        get_transducer_box(trans_xyz([1,3]), target_xyz([1,3]), pixel_size, parameters);

    subplot(1,3,2);

        show_3d_head(img, ...
            target_xyz, ...
            trans_xyz, ...
            parameters, ...
            pixel_size, ...
            coord_mesh.xyz, ...
            [0 0 0],...
            [0,0],...
            0)

    subplot(1,3,3);

        show_3d_head(img, ...
            target_xyz, ...
            trans_xyz, ...
            parameters, ...
            pixel_size, ...
            coord_mesh.xyz, ...
            [0 0 0],...
            [0,0],...
            0)

        view([-175,0])
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_heuristic_%s.png', subject_id, target_name));
    saveas(h, output_plot, 'png')
    close(h);
