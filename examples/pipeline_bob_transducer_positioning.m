clear
%% This is the rough and poorly commented pipeline to get the optimal positions for transducers.
% It assumes that you have already did the segmentation for your subjects in SimNIBS (can be done through the main pipeline)
% The script uses gpuArray so run it on a CUDA-enabled machine. 
% add paths
addpath /project/3015999.02/andche_sandbox/simnibs/simnibs_env/lib/python3.7/site-packages/simnibs/matlab
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/ % change path here

mni_targets = struct('amygdala',[-26 -4 -20], 'ventral_striatum',[-8 6 -10], 'right_mediodorsal_thalamus', [6.5 -18 9], 'right_pulvinar', [15 -27 6.5] );
all_targets = fieldnames(mni_targets);
parameters = load_parameters('sjoerd_config_opt_CTX250-011_64.5mm.yaml');
out_folder = parameters.data_path+'/sim_outputs/';

all_subjs = [3,4,8,9,10];
for subject_id = all_subjs
    fprintf('Current subject: %03i\n', subject_id)

    headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
    filename_segmented_headreco = fullfile(headreco_folder, sprintf('sub-%03d_masks_contr.nii.gz', subject_id));

    segmented_img_orig = niftiread(filename_segmented_headreco);
    segmented_img_head = niftiinfo(filename_segmented_headreco);        
    pixel_size = mean(segmented_img_head.PixelDimensions);
    

    im_center = round(size(segmented_img_orig)/2);

    for target_idx = 1:length(all_targets)
        target_name = all_targets{target_idx};
        fprintf('Current target: %s\n', target_name)
        tpos_output_file = sprintf('%stpars_subj%03i_%s.csv', parameters.data_path, subject_id, target_name);
        if exist(tpos_output_file,'file')
            continue
        end

        simnibs_coords = mni2subject_coords(mni_targets.(target_name), sprintf('%s/m2m_sub-%03i', parameters.data_path, subject_id))
        target = round(transformPointsInverse(segmented_img_head.Transform, simnibs_coords))

        figure
        montage({squeeze(segmented_img_orig(im_center(1),:,:)),...
            squeeze(segmented_img_orig(:,im_center(2),:)),...
            squeeze(segmented_img_orig(:,:,im_center(3)))}, viridis(8), 'Size',[1 3])

        % get list of coordinates at a given distance

        [t1_x, t1_y, t1_z] = ndgrid(1:size(segmented_img_orig, 1),1:size(segmented_img_orig, 2),1:size(segmented_img_orig, 3));
        grid_dist = sqrt((t1_x-target(1)).^2+(t1_y-target(2)).^2+(t1_z-target(3)).^2);
        dist_sphere = abs(grid_dist-parameters.expected_focal_distance_mm/pixel_size)<0.5;
        outer_sphere_3d = dist_sphere&segmented_img_orig==0;

        segm_img_slice = ind2rgb(squeeze(segmented_img_orig(:,target(2),:)), viridis(max(segmented_img_orig(:))+1));
        outer_sphere = squeeze(dist_sphere(:,target(2),:))&squeeze(segmented_img_orig(:,target(2),:))==0;
        segm_img_slice(outer_sphere) = 1;

        figure
        colormap([0.3 0.3 0.3; lines(12)])

        %imagesc(squeeze(grid_dist(:,target(2),:)));
        imagesc(segm_img_slice);
        axis image
        hold on
        target_xz = target([1,3]);

        rectangle('Position',[flip(target_xz) - 2, 4 4],...
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
        max_od_mm = max(parameters.transducer.Elements_OD_mm);

        % normal vector
        dist_gf_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od_mm^2);

        norm_v = (trans_pos-target)/norm(target-trans_pos);
        geom_focus_pos = trans_pos - norm_v*(parameters.transducer.curv_radius_mm)/pixel_size;
        ex_plane_pos = geom_focus_pos+norm_v*dist_gf_to_ep_mm/pixel_size;
        max_od_grid = max(parameters.transducer.Elements_OD_mm)/pixel_size;

        d = sum(norm_v.*ex_plane_pos);
        orth_plane_disk = abs(t1_x*norm_v(1)+t1_y*norm_v(2)+t1_z*norm_v(3)-d)<0.5 & sqrt((t1_x-ex_plane_pos(1)).^2+(t1_y-ex_plane_pos(2)).^2+(t1_z-ex_plane_pos(3)).^2) < max_od_grid/2;

        skin_only = uint8(segmented_img_orig==5);

        skin_only(orth_plane_disk)=2;
        skin_only(outer_sphere_3d) = 3;

        figure
        imagesc(squeeze(skin_only(:,target(2),:)))
        hold on
        trans_xz = trans_pos([1,3]);
        ex_pl_xz = ex_plane_pos([1,3]);
        rectangle('Position',[flip(trans_xz) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','b',...
                 'LineWidth',2,'LineStyle','-');
        rectangle('Position',[flip(geom_focus_pos([1,3])) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','yellow',...
                 'LineWidth',2,'LineStyle','-');
        rectangle('Position',[flip(ex_pl_xz) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','white',...
                 'LineWidth',2,'LineStyle','-');

        rectangle('Position',[flip(target_xz) - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','r',...
                 'LineWidth',2,'LineStyle','-');
        line([trans_xz(2) target_xz(2)], [trans_xz(1) target_xz(1)], 'Color', 'white')
        get_transducer_box(trans_xz, target_xz, pixel_size, parameters)

        all_masks = segmented_img_orig>0;
        all_masks_dilated = imdilate((all_masks), strel('sphere', 1));
        outer_boundary =  all_masks_dilated-all_masks;
        imagesc(squeeze(outer_boundary(:,target(2),:)))
        all_masks_eroded = imerode(all_masks, strel('sphere', 1));
        skin_boundary =  all_masks-all_masks_eroded;
        imagesc(squeeze(skin_boundary(:,target(2),:)))

        skull = segmented_img_orig==4;
        skull_boundary = imfill(segmented_img_orig>0&segmented_img_orig<=4,'holes');
        skull_boundary = imerode(skull_boundary, strel('sphere', 1));
        skull_boundary = skull & ~skull_boundary;
        figure
        imagesc(squeeze(skull_boundary(:,target(2),:)))


        csf_box = regionprops3(segmented_img_orig==3); 
        csf_box = csf_box.BoundingBox; % first 3 numbers are x,y,z of a corner, the other three are dimensions 

        csf_box(1:3) = csf_box(1:3) - 35;
        csf_box(4:6) = csf_box(4:6) + 70;
        csf_box(1:3)
        inside_box = zeros(size(segmented_img_orig));
        inside_box((t1_x>csf_box(1))&(t1_x<(csf_box(1)+csf_box(4)))&...
            (t1_y>csf_box(2))&(t1_y<(csf_box(2)+csf_box(5)))&...
            (t1_z>csf_box(3))&(t1_z<(csf_box(3)+csf_box(6))))=1;
        tmp = segmented_img_orig;
        tmp(~inside_box) = 0;
        imagesc(squeeze(tmp(:,target(2),:)))

        outer_boundary(~inside_box) = 0;
        skin_boundary(~inside_box) = 0;
        imagesc(squeeze(skin_boundary(:,target(2),:)))

        % for each point, put the transducer there, oriented towards the focus, and compute the amount of
        % intersection and the average distance to the scalp

        coord_mesh = struct;
        coord_mesh.x = gpuArray(t1_x);
        coord_mesh.y = gpuArray(t1_y);
        coord_mesh.z = gpuArray(t1_z);

        coord_mesh.xyz = gpuArray([reshape(t1_x,[],1) reshape(t1_y,[],1) reshape(t1_z,[],1)]);

        skin_boundary_coords = gpuArray(coord_mesh.xyz(find(skin_boundary),:));
        skull_boundary_coords = gpuArray(coord_mesh.xyz(find(skull_boundary),:));

        outer_idx = find(outer_boundary);

        coord_rel_to_targ = coord_mesh.xyz - target;
        distances_to_target = sqrt(sum(coord_rel_to_targ.^2,2));
        distances_to_target = distances_to_target(outer_idx);

        close_enough_idx = outer_idx(distances_to_target<(80/pixel_size));
        %close_enough_idx_yplane = intersect(close_enough_idx, find(t1_y==target(2)));

        trans_pos_coords = coord_mesh.xyz(close_enough_idx,:);

        norm_v = gpuArray((trans_pos_coords-target)./repmat(sqrt(sum((trans_pos_coords-target).^2,2)),[1, 3]));

        max_od_mm = max(parameters.transducer.Elements_OD_mm);
        dist_gf_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2-max_od_mm^2);
        dist_tp_to_ep_mm = parameters.transducer.curv_radius_mm - dist_gf_to_ep_mm;
        %parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm
        pos_shift_mm = 5 + dist_tp_to_ep_mm;

        shifted_trans_pos_coords = trans_pos_coords + norm_v*pos_shift_mm/pixel_size;

        geom_focus_pos_all = shifted_trans_pos_coords - norm_v*(parameters.transducer.curv_radius_mm)/pixel_size;
        ex_plane_pos_all = gpuArray(geom_focus_pos_all+norm_v*dist_gf_to_ep_mm/pixel_size);
        all_masks_indx_gpu = gpuArray(find(all_masks>0));

        max_od_mm = max(parameters.transducer.Elements_OD_mm);
        max_od_grid = max_od_mm / pixel_size;

        [prop_intersect, mean_dts, var_dts, mean_dist_skull, var_dist_skull] = arrayfun(@(x) analyze_transducer_position_fast(x, norm_v, ex_plane_pos_all, coord_mesh, all_masks_indx_gpu, skin_boundary_coords, skull_boundary_coords, max_od_grid ), 1:length(close_enough_idx));

        tpos_pars = array2table(gather([close_enough_idx, shifted_trans_pos_coords, repmat(target,[length(close_enough_idx),1]), pdist2( shifted_trans_pos_coords, target), [prop_intersect; mean_dts; var_dts; mean_dist_skull; var_dist_skull]']),...
            'VariableNames',["idx","trans_x","trans_y","trans_z","targ_x","targ_y","targ_z","dist_to_target","prop_intersect","mean_dist_skin","var_dist_skin","mean_dist_skull", "var_dist_skull"]);
        writetable(tpos_pars, tpos_output_file, 'Delimiter',',')
        continue % comment this line for some extra plots of the "best" results 

        sort_var = var_dist_skull;
        sort_var(prop_intersect>0.1) = max(var_dist_skull)+100;
        [v, sort_idx] = sort(sort_var);
        i = sort_idx(1);

        tppf = tpos_pars(tpos_pars.prop_intersect<0.05,:);
        tppf = tppf(tppf.mean_dist_skull <= quantile(tppf.mean_dist_skull, 0.5) & tppf.var_dist_skull <= quantile(tppf.var_dist_skull, 0.1),:);
        tppf = tppf(tppf.var_dist_skin==min(tppf.var_dist_skin),:);

        i = find(close_enough_idx==tppf.idx(1));
        % empty_t1 = zeros(size(outer_boundary));
        % empty_t1(intersect(close_enough_idx, find(t1_y==target(2)))) = 1;
        %empty_t1 = empty_t1 + skin_boundary*2;
        %i = find(close_enough_idx == randsample(intersect(close_enough_idx, find(t1_y==target(2))),1));

        d = sum(norm_v(i,:).*ex_plane_pos_all(i,:));

        % empty_t1(orth_plane_disk) = 1;

        trans_xyz = gather(shifted_trans_pos_coords(i,:));
        ex_pl_xyz = gather(ex_plane_pos_all(i,:));
        gf_xyz = gather(geom_focus_pos_all(i,:));

        [rotated_img, trans_xyz, target_xyz, transformation_matrix, rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
            align_to_focus_axis_and_scale(segmented_img_orig, segmented_img_head, round(trans_xyz)', target', 1, parameters);

        TF = maketform('affine', transformation_matrix);

        ex_pl_xyz = gather(ex_plane_pos_all(i,:));
        gf_xyz = gather(geom_focus_pos_all(i,:));
        ex_pl_xyz = round(tformfwd([ex_pl_xyz]', TF));
        gf_xyz  = round(tformfwd([gf_xyz]', TF));

        figure
        colormap([0.3 0.3 0.3; lines(5)])
        imagesc(squeeze(rotated_img(:,round(trans_xyz(2)),:)))
        rectangle('Position',[target_xyz([3,1])' - 2, 4 4],...
                  'Curvature',[0,0], 'EdgeColor','r',...
                 'LineWidth',2,'LineStyle','-');

        rectangle('Position',[trans_xyz([3,1])' - 2, 4 4],...
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

        orth_plane_disk = find(abs(sum(coord_mesh.xyz.*norm_v(i,:),2)-d)<0.5);
        orth_plane_disk = orth_plane_disk(sqrt(sum((coord_mesh.xyz(orth_plane_disk,:)-ex_plane_pos_all(i,:)).^2,2)) < max_od_grid /2);

        orth_plane_3d = zeros(size(segmented_img_orig));
        orth_plane_3d(orth_plane_disk) = 1;
        orth_plane_3d = smooth3(orth_plane_3d);
        skin_idx = find(skin_boundary);
        max_od_grid = max(parameters.transducer.Elements_OD_mm)/pixel_size;

        show_3d_scalp(segmented_img_orig, target, shifted_trans_pos_coords(i,:), parameters, pixel_size, coord_mesh.xyz, [-1 1 0])

        [t1_x_r, t1_y_r, t1_z_r] = ndgrid(1:size(rotated_img, 1),1:size(rotated_img, 2),1:size(rotated_img, 3));
        t1_xyz_r = gpuArray([reshape(t1_x_r,[],1) reshape(t1_y_r,[],1) reshape(t1_z_r,[],1)]);

        sphere_3d =  zeros(size(rotated_img));
        sphere_3d(pdist2(t1_xyz_r, target_xyz')<3) = 1;



        norm_v_r = (trans_xyz-target_xyz)/norm(trans_xyz-target_xyz);
        norm_v_r = norm_v_r';
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
        figure
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

        % tpos_pars = readtable('tpars_ventral_striatum.csv')
    end
end