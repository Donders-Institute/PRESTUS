function create_group_MNI_plots(subject_list, parameters, options)
% CREATE_GROUP_MNI_PLOTS  Generate group-level simulation plots in MNI space
%
% Processes per-subject ISPPA, pressure, and heating maps, extracts
% orthogonal slices or peak-intensity slices, and generates composite
% plots with optional ROI/FWHM overlays and statistics. Saves montaged
% images and updates group CSV tables.
%
% Use as:
%   create_group_MNI_plots(subject_list, parameters)
%   create_group_MNI_plots(subject_list, parameters, options)
%
% Supported features:
%   - Wildcard support in filenames via parameters.io.output_affix
%   - Heatmap plotting, FWHM mask logic, ROI overlays and statistics
%
% Input:
%   subject_list - array of subject IDs (numeric)
%   parameters   - (1,1) simulation parameters struct; required fields:
%                    .path.temp_output_dir  — path to subject image outputs
%                    .layers.brain          — label values for brain tissue mask
%                    .layers.water          — label values for water tissue mask
%                    .io.output_affix       — filename suffix/wildcard (e.g. '_ses-*')
%                    .simulation.medium     — used in table filenames
%                    .path.seg_path / .path.data_path — SimNIBS m2m segmentation paths
%                    .thermal.temp_0        — baseline temp; scalar or struct with .skin/.water
%   options      - name-value options:
%                    ROI_MNI_mask (3D logical) — binary ROI volume in MNI space
%                    slice_to_plot (int)        — slice index to display
%                    plot_max_intensity (bool)  — use slice with peak ISPPA (conflicts with slice_to_plot)
%                    slice_label ('x'|'y'|'z') — axis to slice along
%                    rotation (deg)             — image rotation for display
%                    plot_heating (bool)        — include heating map plots
%                    outputs_suffix (string)    — suffix for output files
%                    intensity_thresholds ([lo hi]) — manual ISPPA colorbar range
%                    add_FWHM_boundary (bool)   — outline ISPPA > half-max region
%                    add_ROI_boundary (bool)    — overlay and compute ROI stats
%                    skip_missing (bool)        — skip subjects with missing files
%                    brightness_correction (bool) — normalise brightness across subjects
%                    average_target_brightness (float) — target mean brightness value
%
% Output:
%   - ISPPA and optionally heating plots per subject in MNI space
%   - CSV tables updated with ROI and FWHM statistics
%   - Montaged image per modality across all subjects
%
% See also: COMBINE_PLOTS_BY_SUFFIX

arguments
    subject_list
    parameters   (1,1) struct
    options.ROI_MNI_mask (:,:,:)
    options.slice_to_plot = 0 
    options.plot_max_intensity = 0
    options.slice_label = 'y'
    options.rotation = 90;
    options.plot_heating = 1
    options.outputs_suffix = ''
    options.intensity_thresholds = []
    options.add_FWHM_boundary = 0
    options.add_ROI_boundary = 1
    options.skip_missing = 0
    options.brightness_correction = 0
    options.average_target_brightness = 100
end

assert(xor(options.plot_max_intensity, options.slice_to_plot), ...
    "You should indicate either the slice number to plot ('slice_to_plot') or ask for a max intensity plot ('plot_max_intensity')");
slice_labels = {'x','y','z'};

outputs_path = parameters.temp_output_dir;
bg_range_to_use = [];
intensity_range_to_use = [5, 6];
if isstruct(parameters.thermal.temp_0)
    temp_range_to_use = [parameters.thermal.temp_0.skin, parameters.thermal.temp_0.water + 0.5];
    temp_0_min = parameters.thermal.temp_0.skin;
else
    temp_range_to_use = [parameters.thermal.temp_0, parameters.thermal.temp_0 + 0.5];
    temp_0_min = parameters.thermal.temp_0;
end

full_subject_list = subject_list;

% === FIRST PASS: file existence, background and threshold discovery ===
for subject_i = 1:length(full_subject_list)
    subject_id = full_subject_list(subject_i);

    % -- Setup subject-specific directory and prefix logic --
    subject_dir = sprintf('sub-%03d', subject_id);
    file_base   = sprintf('sub-%03d', subject_id);
    data_dir    = fullfile(outputs_path, subject_dir);

    fprintf('Subject %i, first pass\n', subject_id)

    % -- Determine m2m folder for T1 reference --
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg)
        m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', subject_id));
    else
        m2m_folder = fullfile(parameters.path.anat, sprintf('m2m_sub-%03d', subject_id));
    end

    % -- Construct filename patterns (with possible wildcards) --
    pattern_intensity      = fullfile(data_dir, sprintf('%s_final_intensity_MNI%s.nii.gz',              file_base, parameters.io.output_affix));
    pattern_segmented  = fullfile(data_dir, sprintf('%s_final_medium_masks_MNI%s.nii.gz',        file_base, parameters.io.output_affix));
    pattern_pressure   = fullfile(data_dir, sprintf('%s_final_pressure_MNI%s.nii.gz',            file_base, parameters.io.output_affix));
    pattern_output_tbl = fullfile(data_dir, sprintf('%s_%s_output_table%s.csv',                  file_base, parameters.simulation.medium, parameters.io.output_affix));
    pattern_heating    = fullfile(data_dir, sprintf('%s_final_heating_MNI%s.nii.gz',             file_base, parameters.io.output_affix));

    % -- Find files matching patterns --
    files_intensity     = dir(pattern_intensity);
    files_segmented = dir(pattern_segmented);
    files_pressure  = dir(pattern_pressure);
    files_outtbl    = dir(pattern_output_tbl);
    files_heating   = dir(pattern_heating);

    % -- T1 anatomical reference, not wildcard --
    t1_mni_file = fullfile(m2m_folder, 'toMNI','final_tissues_MNI.nii.gz');
    
    files_to_check = {t1_mni_file};
    if isempty(files_intensity),     files_to_check{end+1} = pattern_intensity; end
    if isempty(files_segmented), files_to_check{end+1} = pattern_segmented; end
    if isempty(files_pressure),  files_to_check{end+1} = pattern_pressure; end
    if isempty(files_outtbl),    files_to_check{end+1} = pattern_output_tbl; end
    if options.plot_heating && isempty(files_heating), files_to_check{end+1} = pattern_heating; end

    if numel(files_to_check) > 1
        if options.skip_missing
            fprintf('Missing files for subject %d; skipping this subject\n', subject_id);
            subject_list(subject_i) = NaN;
            continue;
        else
            error('Missing required files for subject %d: %s; quitting\n', subject_id, strjoin(files_to_check(2:end), ', '));
        end
    end

    % -- SUPPORT MULTIPLE FILES PER SUBJECT: loop over all matches --
    for fidx = 1:numel(files_intensity)
        intensity_map_mni_file      = fullfile(data_dir, files_intensity(fidx).name);
        segmented_image_mni_file= fullfile(data_dir, files_segmented(min(fidx, numel(files_segmented))).name);
        max_pressure_mni_file   = fullfile(data_dir, files_pressure(min(fidx, numel(files_pressure))).name);
        output_pressure_file    = fullfile(data_dir, files_outtbl(min(fidx, numel(files_outtbl))).name);
        if options.plot_heating && ~isempty(files_heating)
            heating_data_mni_file = fullfile(data_dir, files_heating(min(fidx, numel(files_heating))).name);
        end

        % -- Read volumetric NIFTI and CSV data --
        t1_mni             = niftiread(t1_mni_file);
        intensity_map_mni      = niftiread(intensity_map_mni_file);
        segmented_image_mni= niftiread(segmented_image_mni_file);
        max_pressure_map_mni = niftiread(max_pressure_mni_file);

        % -- Logical masks: brain extraction for region calculations --
        brain_ind = parameters.layers.brain;
        new_brain_ind = ones(1, length(brain_ind));
        results_mask_original = changem_vectorized(segmented_image_mni, new_brain_ind, brain_ind);
        results_mask_original(results_mask_original > max(brain_ind)) = 0;
        results_mask = logical(results_mask_original);

        % -- Mask application --
        intensity_map_mni = intensity_map_mni .* results_mask;
        max_pressure_map_mni = max_pressure_map_mni .* results_mask;
        max_pressure = max(max_pressure_map_mni,[],'all');

        % -- Slice selection logic --
        if options.slice_to_plot
            slice_n = options.slice_to_plot;
        else
            [~, I] = max(intensity_map_mni(:));
            [Px, Py, Pz] = ind2sub(size(intensity_map_mni), I);
            max_focus_MNI_grid = [Px, Py, Pz];
            slice_n = max_focus_MNI_grid(strcmp(slice_labels,options.slice_label));
            disp(max_focus_MNI_grid)
        end
        
        t1_slice = get_slice_by_label(t1_mni, options.slice_label, slice_n);
        intensity_slice = get_slice_by_label(intensity_map_mni, options.slice_label, slice_n);
        max_pressure_slice = get_slice_by_label(max_pressure_map_mni, options.slice_label, slice_n);
        bg_min = min(t1_slice,[],'all');
        bg_max = max(t1_slice,[],'all');
        if isempty(bg_range_to_use)
            bg_range_to_use = [bg_min, bg_max];
        else
            if bg_min > bg_range_to_use(1)
                bg_range_to_use(1) = bg_min;
            end
            if bg_max < bg_range_to_use(2)
                bg_range_to_use(2) = bg_max;
            end
        end
        max_isppa = max(intensity_slice(:));
        overlay_threshold_low = min(intensity_slice(max_pressure_slice>=(max_pressure*0.4)));

        if overlay_threshold_low < intensity_range_to_use(1)
            intensity_range_to_use(1) = overlay_threshold_low;
        end
        if max_isppa > intensity_range_to_use(2)
            intensity_range_to_use(2) = max_isppa;
        end

        if isfield(options,'plot_heating') && options.plot_heating == 1 && exist('heating_data_mni_file','var')
            maxT_mni = niftiread(heating_data_mni_file);
            water_ind = parameters.layers.water;
            new_water_ind = zeros(1, length(water_ind));
            heating_mask = logical(changem_vectorized(segmented_image_mni, new_water_ind, water_ind));
            maxT_mni = maxT_mni .* heating_mask;
            maxT_mni(maxT_mni==0) = temp_0_min;

            maxT_slice = get_slice_by_label(maxT_mni, options.slice_label, slice_n);
            maxT_in_slice = max(maxT_slice(:));
            if maxT_in_slice > temp_range_to_use(2)
                temp_range_to_use(2) = maxT_in_slice;
            end
        end
    end
end

% --- Remove unavailable subjects after first pass ---
subject_list(isnan(subject_list)) = [];

% === SECOND PASS: generate images and combine output ===
for subject_i = 1:length(subject_list)
    subject_id = subject_list(subject_i);

    subject_dir = sprintf('sub-%03d', subject_id);
    file_base   = sprintf('sub-%03d', subject_id);
    data_dir    = fullfile(outputs_path, subject_dir);

    fprintf('Subject %i, second pass\n', subject_id)

    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg)
        m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', subject_id));
    else
        m2m_folder = fullfile(parameters.path.anat, sprintf('m2m_sub-%03d', subject_id));
    end

    % -- Patterns for batch processing, as before --
    pattern_intensity         = fullfile(data_dir, sprintf('%s_final_intensity_MNI%s.nii.gz',              file_base, parameters.io.output_affix));
    pattern_segmented     = fullfile(data_dir, sprintf('%s_final_medium_masks_MNI%s.nii.gz',        file_base, parameters.io.output_affix));
    pattern_pressure      = fullfile(data_dir, sprintf('%s_final_pressure_MNI%s.nii.gz',            file_base, parameters.io.output_affix));
    pattern_output_tbl    = fullfile(data_dir, sprintf('%s_%s_output_table%s.csv',                  file_base, parameters.simulation.medium, parameters.io.output_affix));
    pattern_output_tbl_roi= fullfile(data_dir, sprintf('%s_%s_output_table_with_ROI_analysis%s.csv',file_base, parameters.simulation.medium, parameters.io.output_affix));
    pattern_heating       = fullfile(data_dir, sprintf('%s_final_heating_MNI%s.nii.gz',             file_base, parameters.io.output_affix));

    files_intensity     = dir(pattern_intensity);
    files_segmented = dir(pattern_segmented);
    files_pressure  = dir(pattern_pressure);
    files_outtbl    = dir(pattern_output_tbl);
    files_outtbl_roi= dir(pattern_output_tbl_roi);
    files_heating   = dir(pattern_heating);

    if isempty(files_intensity) || isempty(files_segmented) || isempty(files_pressure) || isempty(files_outtbl)
        if options.skip_missing
            fprintf('Missing files for subject %d; skipping this subject\n', subject_id);
            continue;
        else
            error('Missing required files for subject %d; quitting\n', subject_id);
        end
    end

    for fidx = 1:numel(files_intensity)
        intensity_map_mni_file       = fullfile(data_dir, files_intensity(fidx).name);
        segmented_image_mni_file = fullfile(data_dir, files_segmented(min(fidx, numel(files_segmented))).name);
        max_pressure_mni_file    = fullfile(data_dir, files_pressure(min(fidx, numel(files_pressure))).name);
        output_pressure_file     = fullfile(data_dir, files_outtbl(min(fidx, numel(files_outtbl))).name);
        if ~isempty(files_outtbl_roi)
            output_pressure_file_with_ROI_analysis = fullfile(data_dir, files_outtbl_roi(min(fidx, numel(files_outtbl_roi))).name);
        else
            output_pressure_file_with_ROI_analysis = '';
        end
        if options.plot_heating && ~isempty(files_heating)
            heating_data_mni_file = fullfile(data_dir, files_heating(min(fidx, numel(files_heating))).name);
        end

        t1_mni_file = fullfile(m2m_folder, 'toMNI','final_tissues_MNI.nii.gz');
        if ~exist(t1_mni_file, 'file')
            if options.skip_missing
                fprintf('Missing T1 file for subject %d; skipping this subject\n', subject_id);
                continue;
            else
                error('Missing T1 file for subject %d; quitting\n', subject_id);
            end
        end

        % -- Full image/statistics logic from your original code --
        t1_mni = niftiread(t1_mni_file);
        t1_mni_hdr = niftiinfo(t1_mni_file);
        intensity_map_mni = niftiread(intensity_map_mni_file);
        segmented_image_mni = niftiread(segmented_image_mni_file);
        max_pressure_map_mni = niftiread(max_pressure_mni_file);

        if isfield(options,'brightness_correction') && options.brightness_correction == 1
            if isfield(options,'average_target_brightness') && isnumeric(options.average_target_brightness)
                slice_test = get_slice_by_label(t1_mni, options.slice_label, slice_n);
                non_zero_indexes = slice_test ~= 0;
                slice_average = mean(slice_test(non_zero_indexes));
                brightness_correction_factor = options.average_target_brightness / slice_average;
                t1_mni = t1_mni * brightness_correction_factor;
            else
                error('Please fill in a numeric value for "average_target_brightness"');
            end
        end

        brain_ind = parameters.layers.brain;
        new_brain_ind = ones(1, length(brain_ind));
        results_mask_original = changem_vectorized(segmented_image_mni, new_brain_ind, brain_ind);
        results_mask_original(results_mask_original > max(brain_ind)) = 0;
        results_mask = logical(results_mask_original);

        intensity_map_mni = intensity_map_mni .* results_mask;
        max_pressure_map_mni = max_pressure_map_mni .* results_mask;
        max_pressure = max(max_pressure_map_mni,[],'all');

        if isempty(options.intensity_thresholds)
            overlay_threshold_high = min(intensity_map_mni(max_pressure_map_mni >= max_pressure*0.5));
            overlay_threshold_low = min(intensity_map_mni(max_pressure_map_mni >= max_pressure*0.4));
        else
            overlay_threshold_high = options.intensity_thresholds(2);
            overlay_threshold_low = options.intensity_thresholds(1);
        end

        if isfield(options,'plot_heating') && options.plot_heating == 1 && exist('heating_data_mni_file','var')
            maxT_mni = niftiread(heating_data_mni_file);
            water_ind = parameters.layers.water;
            new_water_ind = zeros(1, length(water_ind));
            heating_mask = logical(changem_vectorized(segmented_image_mni, new_water_ind, water_ind));
            maxT_mni = maxT_mni .* heating_mask;
            maxT_mni(maxT_mni==0) = temp_0_min;
        end

        if options.slice_to_plot
            slice_n = options.slice_to_plot;
        else
            [~, I] = max(intensity_map_mni(:));
            [Px, Py, Pz] = ind2sub(size(intensity_map_mni), I);
            max_focus_MNI_grid = [Px, Py, Pz];
            slice_n = max_focus_MNI_grid(strcmp(slice_labels,options.slice_label));
        end

        if subject_id == subject_list(end)
            add_colorbar = 1;
        else
            add_colorbar = 0;
        end

        fwhm_size = sum(max_pressure_map_mni >= max_pressure/2,'all');
        fwhm_mask = max_pressure_map_mni >= max_pressure/2;

        output_table = readtable(output_pressure_file);
        output_table.('fwhm_size_MNI_based_on_pressure') = fwhm_size;
        current_rotation = options.rotation;

        % -- ROI and statistics
        if isfield(options,'ROI_MNI_mask') && options.add_ROI_boundary == 1
            roi_size = sum(options.ROI_MNI_mask,'all');
            avg_intensity_within_roi = mean(intensity_map_mni(logical(options.ROI_MNI_mask)),'all');
            if ~isequal((logical(options.ROI_MNI_mask.*fwhm_mask)), zeros(size(options.ROI_MNI_mask)))
                avg_intensity_within_fwhm_and_roi_overlap = mean(intensity_map_mni(logical(options.ROI_MNI_mask.*fwhm_mask)),'all');
            else
                avg_intensity_within_fwhm_and_roi_overlap = mean(avg_intensity_within_roi, 'all');
            end
            output_table.(sprintf('avg_intensity_within_fwhm_and_roi_overlap%s', options.outputs_suffix)) = avg_intensity_within_fwhm_and_roi_overlap;
            n_voxels_within_roi_above_thresh = sum(options.ROI_MNI_mask & (max_pressure_map_mni >= max_pressure/2),'all');
            [cx,cy,cz] = ndgrid(1:size(options.ROI_MNI_mask,1), 1:size(options.ROI_MNI_mask,2), 1:size(options.ROI_MNI_mask,3));
            w = double(options.ROI_MNI_mask);
            w_sum = sum(w(:));
            weighted_centroid = [sum(cx(:).*w(:)), sum(cy(:).*w(:)), sum(cz(:).*w(:))] / w_sum;
            dist_between_intensity_and_center_of_ROI = norm(max_focus_MNI_grid - weighted_centroid);
            output_table.(sprintf('dist_between_intensity_and_center_of_ROI%s', options.outputs_suffix)) = dist_between_intensity_and_center_of_ROI;
            output_table.(sprintf('avg_intensity_within_roi%s', options.outputs_suffix)) = avg_intensity_within_roi;
            output_table.(sprintf('perc_voxels_within_roi%s', options.outputs_suffix)) = n_voxels_within_roi_above_thresh/roi_size;
            output_table.(sprintf('perc_voxels_within_fwhm%s', options.outputs_suffix)) = n_voxels_within_roi_above_thresh/fwhm_size;
            output_table.(sprintf('roi_size%s', options.outputs_suffix)) = roi_size;
            %writetable(output_table, output_pressure_file_with_ROI_analysis); 
        end

        % -- Plot creation step --
        plot_overlay(...
            intensity_map_mni, ...
            t1_mni, ...
            zeros(size(t1_mni)), ...
            struct(), ...
            {options.slice_label, slice_n}, ...
            uint8.empty([0 3]), ...
            uint8.empty([0 3]), ...
            [0 0 0], ...
            'overlay_threshold_low', overlay_threshold_low, ...
            'overlay_threshold_high', overlay_threshold_high, ...
            'show_rectangles', 0, ...
            'overlay_color_range', intensity_range_to_use, ...
            'grid_step', t1_mni_hdr.PixelDimensions(1), ...
            'overlay_segmented', 0, ...
            'rotation', current_rotation, ...
            'show_colorbar', add_colorbar, ...
            'use_overlay_alpha', 1, ...
            'segmented_img', segmented_image_mni, ...
            'bg_range', bg_range_to_use);

        if isfield(options,'ROI_MNI_mask') && options.add_ROI_boundary == 1
            hold on
            mask_im = get_slice_by_label(options.ROI_MNI_mask, options.slice_label, slice_n);
            visboundaries(imrotate(mask_im, current_rotation))
        end

        if isfield(options,'add_FWHM_boundary') && options.add_FWHM_boundary == 1
            hold on
            mask_im = get_slice_by_label(max_pressure_map_mni >= max_pressure/2, options.slice_label, slice_n);
            visboundaries(imrotate(mask_im, current_rotation), 'Color', 'white','LineStyle', '--','LineWidth',0.5,'EnhanceVisibility',0);
        end

        export_fig(fullfile(data_dir, sprintf('%s_final_intensity_MNI%s%s', file_base, parameters.io.output_affix, options.outputs_suffix)),'-silent','-r320');
        close

        if isfield(options,'plot_heating') && options.plot_heating == 1 && exist('heating_data_mni_file','var')
            plot_overlay(maxT_mni, t1_mni, zeros(size(t1_mni)), struct(), {options.slice_label, slice_n}, ...
                uint8.empty([0 3]), uint8.empty([0 3]), ...
                [0 0 0], 'show_rectangles', 0, 'overlay_color_range', temp_range_to_use, 'grid_step', t1_mni_hdr.PixelDimensions(1),...
                'color_scale', 'inferno', 'rotation', current_rotation, 'show_colorbar', add_colorbar,  'segmented_img', segmented_image_mni, 'bg_range', bg_range_to_use);
            if isfield(options,'ROI_MNI_mask') && options.ROI_MNI_mask == 1
                hold on
                mask_im = get_slice_by_label(options.ROI_MNI_mask, options.slice_label, slice_n);
                visboundaries(imrotate(mask_im, current_rotation))
            end
            export_fig(fullfile(data_dir, sprintf('%s_maxT_MNI%s%s',...
                    file_base, parameters.io.output_affix, options.outputs_suffix)),'-silent','-r320');
            close
        end
    end
end

% --- Assemble output images into montages ---
if isfield(options,'plot_heating') && options.plot_heating == 1
    suffix_list = {sprintf('maxT_MNI%s%s', parameters.io.output_affix, options.outputs_suffix),...
        sprintf('final_intensity_MNI%s%s', parameters.io.output_affix, options.outputs_suffix)};
else
    suffix_list = {sprintf('final_intensity_MNI%s%s', parameters.io.output_affix, options.outputs_suffix)};
end

for suffix_cell = suffix_list
    suffix = suffix_cell{:};
    combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters);
end

fprintf('Completed, see final images in %s\n', convertCharsToStrings(outputs_path))
end
