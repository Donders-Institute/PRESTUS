function transducer_positioning(parameters, pn, target_name, mni_targets)
% TRANSDUCER_POSITIONING  Heuristic transducer placement for MNI targets
%
% Loads a subject segmentation, evaluates candidate positions on the skull
% surface, selects the optimal position using heuristic criteria, and
% exports the result for downstream simulation and Localite navigation.
%
% Use as:
%   transducer_positioning(parameters, pn, target_name, mni_targets)
%
% Input:
%   parameters  - (1,1) simulation configuration struct
%   pn          - (1,1) path names struct
%   target_name - target label string
%   mni_targets - (1,1) MNI target coordinates struct
%
% See also: TRANSDUCER_POSITIONING_START, TP_EVALUATE_CANDIDATE_POSITIONS

arguments
    parameters  (1,1) struct
    pn          (1,1) struct
    target_name (1,:) char
    mni_targets (1,1) struct
end

%% 1. PATHS & VALIDATION

currentLoc = fileparts(mfilename("fullpath"));
% add functions here to detect path setup function
addpath(genpath(fullfile(currentLoc, '..')));

[parameters] = path_log_setup(parameters, get_prestus_path);
subject_id = parameters.subject_id;

%% 2. LOAD SEGMENTATION DATA

m2m_folder = fullfile(pn.seg_path, sprintf('m2m_sub-%03d', subject_id));
filename = fullfile(m2m_folder, 'final_tissues.nii.gz');
img = niftiread(filename);
img = gather(img); % Ensure img is on CPU
img_header = niftiinfo(filename);
voxel_size = mean(img_header.PixelDimensions);

% [DEBUG] plot the segmentation
if parameters.simulation.debug
    h = figure;
    im_center = round(size(img)/2);
    montage({rot90(squeeze(img(im_center(1),:,:))), ...
        rot90(squeeze(img(:,im_center(2),:))), ...
        squeeze(img(:,:,im_center(3)))}, ...
        viridis(8), 'Size', [1 3]);
    saveas(h, fullfile(parameters.io.debug_dir_preproc, ...
        sprintf('sub-%03d_segmentation.png', subject_id)), 'png');
    close(h);
    clear im_center;
end

%% heuristic positioning

% specify output file
tpos_output_file = fullfile(parameters.io.output_dir, ...
    sprintf('tpars_sub-%03i_%s.csv', subject_id, target_name));

if confirm_overwriting(tpos_output_file, parameters)

    %% Convert target from MNI (mm) to subject grid space (voxels)

    fprintf('➤ Target: %s\n', target_name);
    target_mni = mni_targets.(target_name);
    target_vox = transform_coordinates(...
        parameters, target_mni, 'mni', 'grid', img_header);

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

%% [Optional] Remove ear locations

[tpos] = tp_remove_ear_locations(parameters, tpos);

%% Select heuristic transducer position

[trans_pos, target_pos, ~] = ...
    tp_select_heuristic_position(tpos, subject_id, target_name, parameters, img_header);

%% Plot heuristic transducer position

tp_plot_heuristic_position(...
    trans_pos, target_pos, img, img_header, parameters, voxel_size, target_name, subject_id);

fprintf('Heuristic transducer placement: %s (sub-%03d) → %s\n', ...
    target_name, subject_id, tpos_output_file);

%% [Optional] Save T1w image with localite-ready header

if isfield(parameters, 'placement') && isfield(parameters.placement, 'localite') && ...
    isfield(parameters.placement.heuristic, 'save_localite_t1') && parameters.placement.heuristic.save_localite_t1 && ...
    isfield(parameters.path, 'localite') && ~isempty(parameters.path.localite)
    disp("Requested to transform T1 header for localite ...")
    path_t1 = fullfile(m2m_folder, 'T1.nii.gz');
    path_localite_out = fullfile(parameters.path.localite, ...
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
