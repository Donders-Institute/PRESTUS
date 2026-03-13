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
img = niftiread(filename);
img_info = niftiinfo(filename);
pixel_size = mean(img_info.PixelDimensions);

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

    %% 3. Convert target from MNI (mm) to subject grid space (voxels)

    fprintf('➤ Target: %s\n', target_name);
    target_mni = mni_targets.(target_name);
    target_vox = map_coordsystems(...
        parameters, target_mni, 'mni', 'grid', img_info);

    %% 5. Find candidate transducer positions on skull (expanding sphere)

    [trans_candidate, outer_sphere_3d, parameters] = ...
    tp_find_initial_candidate(img, target_vox, pixel_size, parameters);

    %% 6. Plot initial candidate

    tp_plot_candidate_positions(...
        target_vox, trans_candidate, ...
        pixel_size, parameters, subject_id, target_name);

    %% 7. Plot geometry

    tp_plot_geometry_overlay(img, target_vox, trans_candidate.trans_pos, ...
        pixel_size, parameters, subject_id, target_name, outer_sphere_3d);

    %% 8. Evaluate candidate according to criteria

    tpos = tp_evaluate_candidate_positions(...
        img, target_vox, parameters, pixel_size);

    %% 9. Save results table

    writetable(tpos, tpos_output_file, 'Delimiter', ',');

else
    disp('Skipping positioning, loading existing output file...')
    tpos = readtable(tpos_output_file, 'Delimiter', ',');
end

%% [Optional] Convert positions to RAS & MNI

% TO DO

%% [Optional] Remove ear locations

function [tpos] = tp_remove_ear_locations(parameters, tpos);

%% Plot heuristic transducer placement

tp_plot_heuristic_transducer(...
    tpos, img, img_info, parameters, pixel_size, target_name, subject_id);

fprintf('Heuristic transducer placement: %s (sub-%03d) → %s\n', ...
    target_name, subject_id, tpos_output_file);

end
