function create_group_MNI_plots(subject_list, parameters, options)
arguments
    subject_list 
    parameters struct
    options.ROI_MNI_mask (:,:,:)
    options.slice_to_plot = 0 
    options.plot_max_intensity = 0
    options.slice_label = 'y'
    options.rotation = 90;
    options.plot_heating = 1
    options.outputs_suffix = ''
    options.isppa_thresholds = []
    options.add_FWHM_boundary = 0
    options.add_ROI_boundary = 1
    options.skip_missing = 0
end

assert(xor(options.plot_max_intensity,options.slice_to_plot), "You should indicate either the slice number to plot ('slice_to_plot') or ask for a max intensity plot ('plot_max_intensity')")
slice_labels = {'x','y','z'};


outputs_path = parameters.temp_output_dir;
bg_range_to_use = [];
isppa_range_to_use = [5, 6];
temp_range_to_use = [parameters.thermal.temp_0, parameters.thermal.temp_0 + 0.5];    

full_subject_list = subject_list;
% first pass: get background intensity at the target slice & isppa range
% for the colorbar; also checks that the neccessary files exist
for subject_i = 1:length(full_subject_list)
    
    subject_id = full_subject_list(subject_i);
    if parameters.subject_subfolder
        results_prefix = sprintf('sub-%1$03d/sub-%1$03d', subject_id);
    else
        results_prefix = sprintf('sub-%1$03d', subject_id);
    end
    fprintf('Subject %i, first pass\n', subject_id)
    
    headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
    isppa_map_mni_file  = fullfile(outputs_path, sprintf('%s_final_isppa_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    segmented_image_mni_file = fullfile(outputs_path, sprintf('%s_final_medium_masks_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    max_pressure_mni_file = fullfile(outputs_path, sprintf('%s_final_pressure_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    output_pressure_file = fullfile(outputs_path,sprintf('%s_%s_isppa%s.csv', results_prefix, parameters.simulation_medium, parameters.results_filename_affix));
    if strcmp(parameters.segmentation_software, 'headreco')
        t1_mni_file = fullfile(headreco_folder, 'toMNI','T1fs_nu_12DOF_MNI.nii.gz');
    else
        t1_mni_file = fullfile(headreco_folder, 'toMNI','final_tissues_MNI.nii.gz');
    end
    
    files_to_check = {t1_mni_file , segmented_image_mni_file, isppa_map_mni_file, max_pressure_mni_file, output_pressure_file};
    
    if options.plot_heating
        heating_data_mni_file = fullfile(outputs_path,sprintf('%s_final_heating_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
        files_to_check{length(files_to_check)+1} = heating_data_mni_file;    
    end
    all_files_exist = 1;
    for filename = files_to_check
        if ~exist(filename{:}, 'file')
            if options.skip_missing
                sprintf('File does not exist: %s; skipping this subject', filename{:})
                subject_list(subject_i) = [];
                all_files_exist = 0;
                
                break;
            else
                assert(logical(exist(filename{:}, 'file')), sprintf('File does not exist: %s; quitting', filename{:}))
            end
        end
    end
    if ~all_files_exist 
        continue;
    end
    % get T1 in MNI space
    t1_mni = niftiread(t1_mni_file);
    t1_mni_hdr = niftiinfo(t1_mni_file);
    
    Isppa_map_mni = niftiread(isppa_map_mni_file);
    segmented_image_mni = niftiread(segmented_image_mni_file);
    brain_ind = find(strcmp(fieldnames(parameters.layer_labels), 'brain'));
    segmented_mask = segmented_image_mni~=brain_ind;%(segmented_image_mni ==0 | segmented_image_mni ==4 | segmented_image_mni == 5 | segmented_image_mni == 6);
    Isppa_map_mni(segmented_mask) = 0; 
    max_pressure_map_mni = niftiread(max_pressure_mni_file);
    max_pressure_map_mni(segmented_mask) = 0; 
    max_pressure = max(max_pressure_map_mni,[],'all');

    if options.slice_to_plot
        slice_n = options.slice_to_plot;
    else
        [~, I] = max(Isppa_map_mni(:));
        [Px, Py, Pz] = ind2sub(size(Isppa_map_mni), I);
        max_focus_MNI_grid = [Px, Py, Pz];
        slice_n = max_focus_MNI_grid(strcmp(slice_labels,options.slice_label));
        disp(max_focus_MNI_grid)
    end
    t1_slice = get_slice_by_label(t1_mni, options.slice_label, slice_n);
    isppa_slice = get_slice_by_label(Isppa_map_mni, options.slice_label, slice_n);
    max_pressure_slice = get_slice_by_label(max_pressure_map_mni, options.slice_label, slice_n);
    segmented_slice = get_slice_by_label(segmented_image_mni, options.slice_label, slice_n);
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
    max_isppa = max(isppa_slice(:));
    isppa_threshold_low = min(isppa_slice(max_pressure_slice>=(max_pressure*0.4)));

    if isppa_threshold_low<isppa_range_to_use(1)
        isppa_range_to_use(1) = isppa_threshold_low;
    end
    if max_isppa > isppa_range_to_use(2)
        isppa_range_to_use(2) = max_isppa;
    end
    
    if options.plot_heating
        maxT_mni = niftiread(heating_data_mni_file);
        maxT_mni(maxT_mni==0) = parameters.thermal.temp_0;

        maxT_slice = get_slice_by_label(maxT_mni, options.slice_label, slice_n);
        
        maxT_in_slice = max(maxT_slice(:));
        if maxT_in_slice > temp_range_to_use(2)
            temp_range_to_use(2) = maxT_in_slice;
        end
    
    end
end


% second pass, doing the actual images
for subject_i = 1:length(subject_list)
    
    subject_id = subject_list(subject_i);
    if parameters.subject_subfolder
        results_prefix = sprintf('sub-%1$03d/sub-%1$03d', subject_id);
    else
        results_prefix = sprintf('sub-%1$03d', subject_id);
    end
    fprintf('Subject %i, second pass\n', subject_id)
    
    headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
    isppa_map_mni_file  = fullfile(outputs_path, sprintf('%s_final_isppa_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    segmented_image_mni_file = fullfile(outputs_path, sprintf('%s_final_medium_masks_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    max_pressure_mni_file = fullfile(outputs_path, sprintf('%s_final_pressure_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    output_pressure_file = fullfile(outputs_path,sprintf('%s_%s_isppa%s.csv', results_prefix, parameters.simulation_medium, parameters.results_filename_affix));
    if strcmp(parameters.segmentation_software, 'headreco')
        t1_mni_file = fullfile(headreco_folder, 'toMNI','T1fs_nu_12DOF_MNI.nii.gz');
    else
        t1_mni_file = fullfile(headreco_folder, 'toMNI','final_tissues_MNI.nii.gz');
    end

    % get T1 in MNI space
    t1_mni = niftiread(t1_mni_file);
    t1_mni_hdr = niftiinfo(t1_mni_file);
    
    Isppa_map_mni = niftiread(isppa_map_mni_file);
    segmented_image_mni = niftiread(segmented_image_mni_file);
    brain_ind = find(strcmp(fieldnames(parameters.layer_labels), 'brain'));
    segmented_mask = segmented_image_mni~=brain_ind;
    Isppa_map_mni(segmented_mask) = 0; 
    max_pressure_map_mni = niftiread(max_pressure_mni_file);
    max_pressure_map_mni(segmented_mask) = 0; 
    max_pressure = max(max_pressure_map_mni,[],'all');
    
    if isempty(options.isppa_thresholds)
        isppa_threshold_high = min(Isppa_map_mni(max_pressure_map_mni>=max_pressure*0.5));
        isppa_threshold_low = min(Isppa_map_mni(max_pressure_map_mni>=max_pressure*0.4));
    else
        isppa_threshold_high = options.isppa_thresholds(2);
        isppa_threshold_low = options.isppa_thresholds(1);
    end
    if options.plot_heating
        heating_data_mni_file = fullfile(outputs_path,sprintf('%s_final_heating_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
        maxT_mni = niftiread(heating_data_mni_file);
        maxT_mni(maxT_mni==0) = parameters.thermal.temp_0;
    end

    if options.slice_to_plot
        slice_n = options.slice_to_plot;
    else
        [~, I] = max(Isppa_map_mni(:));
        [Px, Py, Pz] = ind2sub(size(Isppa_map_mni), I);
        max_focus_MNI_grid = [Px, Py, Pz];
        slice_n = max_focus_MNI_grid(strcmp(slice_labels,options.slice_label));
    end

    if subject_id == subject_list(end)
        add_colorbar = 1;
    else
        add_colorbar = 0;
    end

    fwhm_size = sum(max_pressure_map_mni >= max_pressure/2,'all');
    fwhm_mask(max_pressure_map_mni ~= 0) = 1;

    curTable = readtable(output_pressure_file);
    curTable.('fwhm_size_MNI_based_on_pressure') = fwhm_size;

    slice = find(strcmp(slice_labels,options.slice_to_plot));
    cur_rotation = options.rotation;
        
    if isfield(options,'ROI_MNI_mask')
        roi_size = sum(options.ROI_MNI_mask,'all');
        avg_isppa_within_roi = mean(Isppa_map_mni(logical(options.ROI_MNI_mask)),'all');
        if ~isequal((logical(options.ROI_MNI_mask.*fwhm_mask)), zeros(size(options.ROI_MNI_mask)))
            avg_isppa_within_fwhm_and_roi_overlap = mean(Isppa_map_mni(logical(options.ROI_MNI_mask.*fwhm_mask)),'all');
        else
            avg_isppa_within_fwhm_and_roi_overlap = mean(avg_isppa_within_roi, 'all');
        end
        curTable.(sprintf('avg_isppa_within_fwhm_and_roi_overlap%s', options.outputs_suffix)) = avg_isppa_within_fwhm_and_roi_overlap;
        n_voxels_within_roi_above_thresh =  sum(options.ROI_MNI_mask & (max_pressure_map_mni >= max_pressure/2),'all');
        props = regionprops(true(size(options.ROI_MNI_mask)), options.ROI_MNI_mask, 'WeightedCentroid');
        dist_between_Isppa_and_center_of_ROI = norm(max_focus_MNI_grid - props.WeightedCentroid);
        curTable.(sprintf('dist_between_Isppa_and_center_of_ROI%s', options.outputs_suffix)) = dist_between_Isppa_and_center_of_ROI;
        curTable.(sprintf('avg_isppa_within_roi%s', options.outputs_suffix)) = avg_isppa_within_roi;
        curTable.(sprintf('perc_voxels_within_roi%s', options.outputs_suffix)) = n_voxels_within_roi_above_thresh/roi_size;
        curTable.(sprintf('perc_voxels_within_fwhm%s', options.outputs_suffix)) = n_voxels_within_roi_above_thresh/fwhm_size;
        curTable.(sprintf('roi_size%s', options.outputs_suffix)) = roi_size;

        writetable(curTable, output_pressure_file);        
    end
 
    plot_isppa_over_image(Isppa_map_mni, t1_mni, zeros(size(t1_mni)), struct(), {options.slice_label, slice_n}, ...
        uint8.empty([0 3]), uint8.empty([0 3]),...
        [0 0 0 ],'isppa_threshold_low',isppa_threshold_low,'isppa_threshold_high',isppa_threshold_high, 'show_rectangles', 0, 'isppa_color_range', isppa_range_to_use, ...
        'grid_step', t1_mni_hdr.PixelDimensions(1), 'overlay_segmented', 0, 'rotation', cur_rotation, 'show_colorbar', add_colorbar, 'use_isppa_alpha', 1, 'segmented_img', segmented_image_mni, 'bg_range', bg_range_to_use);

    if isfield(options,'ROI_MNI_mask') && options.add_ROI_boundary == 1
        hold on
        mask_im = get_slice_by_label(options.ROI_MNI_mask, options.slice_label, slice_n);
        visboundaries(imrotate(mask_im, cur_rotation))
    end

    if isfield(options,'add_FWHM_boundary') && options.add_FWHM_boundary == 1
        hold on
        mask_im = get_slice_by_label(max_pressure_map_mni >= max_pressure/2, options.slice_label, slice_n);
        visboundaries(imrotate(mask_im, cur_rotation), 'Color', 'white','LineStyle', '--','LineWidth',0.5,'EnhanceVisibility',0);
    end

    export_fig(fullfile(outputs_path, sprintf('%s_final_isppa_MNI%s%s',...
            results_prefix, parameters.results_filename_affix, options.outputs_suffix)),'-silent','-r320');
    close

    if options.plot_heating
            
        plot_isppa_over_image(maxT_mni, t1_mni, zeros(size(t1_mni)), struct(), {options.slice_label, slice_n}, ...
            uint8.empty([0 3]), uint8.empty([0 3]),...
            [0 0 0], 'show_rectangles', 0, 'isppa_color_range', temp_range_to_use, 'grid_step', t1_mni_hdr.PixelDimensions(1),...
            'color_scale', 'inferno', 'rotation', cur_rotation, 'show_colorbar', add_colorbar,  'segmented_img', segmented_image_mni, 'bg_range', bg_range_to_use);
            %'ticks', 37:0.05:37.15, 'tick_labels', arrayfun(@(x) sprintf('%.2f', x), 37:0.05:37.15, 'UniformOutput', false));
        if isfield(options,'ROI_MNI_mask')

            hold on
            mask_im = get_slice_by_label(options.ROI_MNI_mask, options.slice_label, slice_n);
            visboundaries(imrotate(mask_im, cur_rotation))
        end
        export_fig(fullfile(outputs_path, sprintf('%s_maxT_MNI%s%s',...
                results_prefix, parameters.results_filename_affix, options.outputs_suffix)),'-silent','-r320');
        close
    end
            
end

% combine individual images in a single image
if options.plot_heating == 1
    suffix_list = {sprintf('maxT_MNI%s%s', parameters.results_filename_affix, options.outputs_suffix),...
        sprintf('final_isppa_MNI%s%s', parameters.results_filename_affix, options.outputs_suffix)
        };
else
    suffix_list = sprintf('final_isppa_MNI%s%s', parameters.results_filename_affix, options.outputs_suffix);
end

for suffix_cell = suffix_list
    suffix = suffix_cell{:};
    combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters);
    
end

fprintf('Completed, see final images in %s\n', convertCharsToStrings(outputs_path))

end