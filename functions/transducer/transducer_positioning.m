function transducer_positioning(parameters, pn, subject_id, target_name, mni_targets)
% TRANSDUCER_POSITIONING Heuristic transducer placement for MNI targets

arguments
    parameters struct
    pn struct
    subject_id double
    target_name string
    mni_targets struct
end

%% 1. PATHS & VALIDATION
currentLoc = fileparts(mfilename("fullpath"));
% add functions here to detect path setup function
addpath(genpath(fullfile(currentLoc, '..')));
[parameters] = path_log_setup(parameters, get_prestus_path, subject_id);

%% 2. LOAD SEGMENTATION DATA
headreco_folder = fullfile(pn.seg_path, sprintf('m2m_sub-%03d', subject_id));
filename = fullfile(headreco_folder, 'final_tissues.nii.gz');
seg_img = niftiread(filename);
segmented_img_head = niftiinfo(filename);
pixel_size = mean(segmented_img_head.PixelDimensions);
im_center = round(size(seg_img)/2);

% [DEBUG] plot the segmentation
if parameters.debug
    h = figure; 
    montage({rot90(squeeze(seg_img(im_center(1),:,:))), ...
        rot90(squeeze(seg_img(:,im_center(2),:))), ...
        squeeze(seg_img(:,:,im_center(3)))}, ...
        viridis(8), 'Size', [1 3]);
    saveas(h, fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_segmentation.png', subject_id)), 'png');
    close(h);
end

%% 3. Convert target_vox from MNI to subject space

fprintf('➤ Target: %s\n', target_name);
target_mni = mni_targets.(target_name);

target_vox = map_coordsystems(parameters, target_mni, 'mni', 'grid', segmented_img_head);

%% 5. Find candidate positions

[outer_sphere_3d, segm_img_slice, t1_x, t1_y, t1_z] = ...
    tp_find_candidate_positions(img_orig, target_vox, pixel_size, parameters, subject_id, target_name);

%% 6. Plot initial candidate positions

tp_plot_initial_candidate_positions(segm_img_slice, target_vox, outer_sphere_3d, t1_x, t1_y, t1_z, ...
    pixel_size, parameters, subject_id, target_name);

%% 7. Plot geometry

tp_plot_geometry_overlay(img_orig, target_vox, outer_sphere_3d, ...
    t1_x, t1_y, t1_z, pixel_size, parameters, subject_id, target_name);

%% 8. Candidate EVALUATION

tpos_pars = tp_evaluate_candidate_positions(seg_img, target_vox, parameters, pixel_size);

%% 9. Optimal transducer VISUALIZATION

tp_visualize_optimal_transducer(tpos_pars, seg_img, segmented_img_head, ...
    parameters, pixel_size, target_name, subject_id);

%% 10. Save Results table

tpos_output_file = fullfile(parameters.output_dir, sprintf('tpars_sub-%03i_%s.csv', subject_id, target_name));
writetable(tpos_pars, tpos_output_file, 'Delimiter', ',');
fprintf('Heuristic transducer placement: %s (sub-%03d) → %s\n', target_name, subject_id, tpos_output_file);

end
