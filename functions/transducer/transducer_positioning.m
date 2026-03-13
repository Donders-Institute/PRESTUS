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

m2m_folder = fullfile(pn.seg_path, sprintf('m2m_sub-%03d', subject_id));
filename = fullfile(m2m_folder, 'final_tissues.nii.gz');
img = niftiread(filename);
img_info = niftiinfo(filename);
voxel_size = mean(img_info.PixelDimensions);

% [DEBUG] plot the segmentation
if parameters.debug
    h = figure;
    im_center = round(size(img)/2);
    montage({rot90(squeeze(img(im_center(1),:,:))), ...
        rot90(squeeze(img(:,im_center(2),:))), ...
        squeeze(img(:,:,im_center(3)))}, ...
        viridis(8), 'Size', [1 3]);
    saveas(h, fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_segmentation.png', subject_id)), 'png');
    close(h);
    clear im_center;
end

%% heuristic positioning

% specify output file
tpos_output_file = fullfile(parameters.output_dir, ...
    sprintf('tpars_sub-%03i_%s.csv', subject_id, target_name));

if confirm_overwriting(tpos_output_file, parameters)

    %% Convert target from MNI (mm) to subject grid space (voxels)

    fprintf('➤ Target: %s\n', target_name);
    target_mni = mni_targets.(target_name);
    target_vox = transform_coordinates(...
        parameters, target_mni, 'mni', 'grid', img_info);

    %% Find candidate transducer positions on skull (expanding sphere)

    [trans_candidate, outer_sphere_3d, parameters] = ...
    tp_find_initial_candidate(img, target_vox, voxel_size, parameters);

    %% Plot initial candidate

    tp_plot_candidate_positions(...
        img, target_vox, trans_candidate, ...
        voxel_size, parameters, subject_id, target_name);

    %% Plot geometry

    tp_plot_geometry_overlay(img, target_vox, trans_candidate.trans_pos, ...
        voxel_size, parameters, subject_id, target_name, outer_sphere_3d);

    %% Evaluate candidate according to criteria

    tpos = tp_evaluate_candidate_positions(...
        img, target_vox, parameters, voxel_size);

    %% Save results table

    writetable(tpos, tpos_output_file, 'Delimiter', ',');

else
    disp('Skipping positioning, loading existing output file...')
    tpos = readtable(tpos_output_file, 'Delimiter', ',');
end

%% [Optional] Convert positions to RAS

% TO DO

%% [Optional] Remove ear locations

[tpos] = tp_remove_ear_locations(parameters, tpos);

%% Select heuristic transducer position

[trans_pos, target_pos, ~] = ...
    tp_select_heuristic_position(tpos, voxel_size, subject_id, target_name, parameters);

%% Plot heuristic transducer position

tp_plot_heuristic_position(...
    trans_pos, target_pos, img, img_info, parameters, voxel_size, target_name, subject_id);

fprintf('Heuristic transducer placement: %s (sub-%03d) → %s\n', ...
    target_name, subject_id, tpos_output_file);

%% [Optional] Save T1w image with localite-ready header

if isfield(parameters, 'tp_save_localiteT1') && parameters.tp_save_localiteT1 && ...
    isfield(parameters, 'localite_path') && ~isempty(parameters.localite_path)
    disp("Requested to transform T1 header for localite ...")
    path_t1 = fullfile(m2m_folder, 'T1.nii.gz');
    path_localite_out = fullfile(parameters.localite_path, ...
        ['sub-', sprintf('%03.0f_T1_forneuronav.nii', subject_id)]);
    if confirm_overwriting(path_localite_out, parameters)
        canonical_affine_transform(path_t1, path_localite_out);
        fprintf('Localite-ready T1 (sub-%03d) → %s\n', subject_id, path_localite_out);
    else
        disp('... but localite T1 file already exists ... skipping.')
    end
end

%% Finish

% capture time, RAM, & GB load of pipeline
log_timer('stop','prestus_pipeline')

% end logging
diary('off')

end
