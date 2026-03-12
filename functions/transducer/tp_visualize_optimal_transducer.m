function tp_visualize_optimal_transducer(tpos_pars, img_orig, img_info, parameters, pixel_size, target_name, subject_id)
% TP_VISUALIZE_OPTIMAL_TRANSDUCER Create 6-panel visualization of optimal transducer placement
%
% Selects optimal position based on intersection, skull/skin distance criteria. Creates comprehensive
% 6-panel figure showing aligned anatomy, geometric properties, 3D rendering, and validation plots.
%
% INPUT
%   tpos_pars    - Table from tp_evaluate_candidate_positions_gpu()
%   img_orig     - Original segmented head image
%   img_info     - NIfTI header info  
%   parameters   - Parameters struct
%   pixel_size   - Scalar voxel size in mm
%   target_name  - String target identifier
%   subject_id   - Scalar subject ID
%
% OUTPUT
%   Saves: sub-XXX_optimal_TARGET.png (6-panel figure)

% Select optimal position (original multi-step criteria)
tppf = tpos_pars(tpos_pars.prop_intersect<0.05,:);
tppf = tppf(tppf.mean_dist_skull <= quantile(tppf.mean_dist_skull, 0.5) & ...
           tppf.var_dist_skull <= quantile(tppf.var_dist_skull, 0.1),:);
tppf = tppf(tppf.var_dist_skin==min(tppf.var_dist_skin),:);
i = 1; % First optimal candidate
target_xyz = [tppf.targ_x(i), tppf.targ_y(i), tppf.targ_z(i)];
trans_xyz = [tppf.trans_x(i), tppf.trans_y(i), tppf.trans_z(i)];

% Align to focal axis (original preprocessing)
[rotated_img, trans_xyz, target_xyz, transformation_matrix, ~, ~, ~, ~] = ...
    preproc_align_to_focal_axis(img_orig, img_info, round(trans_xyz)', target_xyz', 1, parameters);
TF = maketform('affine', transformation_matrix);

% 6-panel comprehensive visualization
h = figure('units', 'normalized', 'position', [0 0 1 1]);
    subplot(2,3,1);
        colormap([0.3 0.3 0.3; lines(5)])
        imagesc(squeeze(rotated_img(:,round(trans_xyz(2)),:)))
        rectangle('Position',[target_xyz([3,1]) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','r',...
                 'LineWidth',2,'LineStyle','-');
        rectangle('Position',[trans_xyz([3,1]) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','b',...
                 'LineWidth',2,'LineStyle','-');
        rectangle('Position',[gf_xyz([3,1]) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','yellow',...
                 'LineWidth',2,'LineStyle','-');
        rectangle('Position',[ex_pl_xyz([3,1]) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','white',...
                 'LineWidth',2,'LineStyle','-');

        line([trans_xyz(3) target_xyz(3)], [trans_xyz(1) target_xyz(1)], 'Color', 'white')
        get_transducer_box(trans_xyz([1,3]), target_xyz([1,3]), pixel_size, parameters)

    subplot(2,3,2);
        colormap([0.3 0.3 0.3; lines(12)])
        %imagesc(squeeze(grid_dist(:,target(2),:)));
        imagesc(segm_img_slice);
        axis image
        hold on
        target_xz = target([1,3]);

        rectangle('Position',[flip(target_xz) - 2, 4, 4],...
                  'Curvature',[0,0], 'EdgeColor','r',...
                 'LineWidth',2,'LineStyle','-');

        outer_idx = find(outer_sphere_3d&t1_y==target(2));

        trans_idx =  randsample(outer_idx, 1);
        trans_pos = [t1_x(trans_idx), t1_y(trans_idx), t1_z(trans_idx)];
        trans_xz = trans_pos([1,3]);
        rectangle('Position',[flip(trans_xz) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','b',...
                 'LineWidth',2,'LineStyle','-');


        get_transducer_box(trans_xz, target_xz, pixel_size, parameters);

    subplot(2,3,3);

        orth_plane_disk = find(abs(sum(coord_mesh.xyz.*norm_v(i,:),2)-d)<0.5);
        orth_plane_disk = orth_plane_disk(sqrt(sum((coord_mesh.xyz(orth_plane_disk,:)-ex_plane_pos_all(i,:)).^2,2)) < max_od_grid /2);

        orth_plane_3d = zeros(size(segmented_img_orig));
        orth_plane_3d(orth_plane_disk) = 1;
        orth_plane_3d = smooth3(orth_plane_3d);
        skin_idx = find(skin_boundary);
        max_od_grid = max(parameters.transducer.Elements_OD_mm)/pixel_size;

        %show_3d_head(segmented_img_orig, target, shifted_trans_pos_coords(i,:), parameters, pixel_size, coord_mesh.xyz, [-1 1 0],[0,0],0)
        show_3d_head(segmented_img_orig, target, shifted_trans_pos_coords(i,:), parameters, pixel_size, coord_mesh.xyz, [0 0 0],[0,0],0)

        [t1_x_r, t1_y_r, t1_z_r] = ndgrid(1:size(rotated_img, 1),1:size(rotated_img, 2),1:size(rotated_img, 3));
        t1_xyz_r = gpuArray([reshape(t1_x_r,[],1) reshape(t1_y_r,[],1) reshape(t1_z_r,[],1)]);

        sphere_3d =  zeros(size(rotated_img));
        sphere_3d(pdist2(t1_xyz_r, target_xyz)<3) = 1;

    subplot(2,3,4);

        norm_v_r = (trans_xyz-target_xyz)/norm(trans_xyz-target_xyz);
        %norm_v_r = norm_v_r';
        d = sum(norm_v_r.*ex_pl_xyz);
        orth_plane_disk = find(abs(sum(t1_xyz_r.*norm_v_r,2)-d)<0.5);

        orth_plane_disk = orth_plane_disk(sqrt(sum((t1_xyz_r(orth_plane_disk,:)-ex_pl_xyz).^2,2)) < max_od_grid /2);
        orth_plane_3d = zeros(size(rotated_img));
        orth_plane_3d(orth_plane_disk) = 1;
        orth_plane_3d = smooth3(orth_plane_3d);

        segmented_img_to_plot = rotated_img;
        % sM = size(segmented_img_to_plot);
        % N = norm_v(i,:);
        % N_orth = null(N(:).')
        % N_orth = N_orth(:,1)
        % segmented_img_to_plot(:,floor(target_xyz(2)):size(segmented_img_to_plot,2),:) = [];
        %segmented_img_to_plot(:,:,floor(target_xyz(3)):size(segmented_img_to_plot,3)) = [];

        % Dist = reshape(((1:sM(1)) - target(1)) * N_orth(1), [sM(1), 1, 1]) + ...
        %        reshape(((1:sM(2)) - target(2)) * N_orth(2), [1, sM(2), 1]) + ...
        %        reshape(((1:sM(3)) - target(3)) * N_orth(3), [1, 1, sM(3)]);
        % 
        % segmented_img_to_plot(Dist>0) = 0;
        segmented_img_to_plot = segmented_img_to_plot(1:2:end,1:2:end,1:2:end);

        for_caps = segmented_img_to_plot;
        for_caps(for_caps>4) = 0;

        Ds = smooth3(gpuArray(double(segmented_img_to_plot>0)));
        skin_isosurface = isosurface(Ds,0.5);
        %figure
        colormap(gray(80))

        hiso = patch(skin_isosurface,...
           'FaceColor',[1,.75,.65],...
           'EdgeColor','none');
        isonormals(Ds,hiso);

        hcap = patch(isocaps(for_caps*10,4),...
           'FaceColor','interp',...
           'EdgeColor','none');
        lightangle(45,30);
        lighting gouraud
        hcap.AmbientStrength = 0.6;
        hiso.SpecularColorReflectance = 0;
        hiso.SpecularExponent = 50;

        trans_obj = patch(isosurface(orth_plane_3d(1:2:end,1:2:end,1:2:end)),...
           'FaceColor','blue',...
           'EdgeColor','none');

        patch(isosurface(smooth3(sphere_3d(1:2:end,1:2:end,1:2:end))),'FaceColor','red','EdgeColor','red')
        view(90,90)

    subplot(2,3,5)
        colormap([0.3 0.3 0.3; lines(12)])
        %imagesc(squeeze(grid_dist(:,target(2),:)));
        imagesc(segm_img_slice);
        axis image
        hold on
        target_xz = target([1,3]);

        rectangle('Position',[flip(target_xz) - 2, 4, 4],...
                  'Curvature',[0,0], 'EdgeColor','r',...
                 'LineWidth',2,'LineStyle','-');

        outer_idx = find(outer_sphere_3d&t1_y==target(2));

        trans_idx =  randsample(outer_idx, 1);
        trans_pos = [t1_x(trans_idx), t1_y(trans_idx), t1_z(trans_idx)];
        trans_xz = trans_pos([1,3]);
        rectangle('Position',[flip(trans_xz) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','b',...
                 'LineWidth',2,'LineStyle','-');


        get_transducer_box(trans_xz, target_xz, pixel_size, parameters);
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_optimal_%s.png', subject_id, target_name));
    saveas(h, output_plot, 'png')
    close(h);
