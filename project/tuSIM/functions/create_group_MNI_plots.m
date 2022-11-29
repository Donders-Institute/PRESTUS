function create_group_MNI_plots(outputs_path, subject_list, parameters, options)
arguments
    outputs_path string 
    subject_list 
    parameters struct
    options.ROI_MNI_mask (:,:,:)
    options.slice_to_plot = 0 
    options.plot_max_intensity = 0
    options.slice_label = 'y'
    options.rotation = 90;
    options.plot_heating = 1
    options.outputs_suffix = ''
    options.font_size = 64
    options.isppa_thresholds = []
    options.add_FWHM_boundary = 0
end

assert(xor(options.plot_max_intensity,options.slice_to_plot), "You should indicate either the slice number to plot ('slice_to_plot') or ask for a max intensity plot ('plot_max_intensity')")
slice_labels = {'x','y','z'};



bg_range_to_use = [];
isppa_range_to_use = [5, 6];
temp_range_to_use = [37, 37.5];    
    
% first pass: get background intensity at the target slice & isppa range
% for the colorbar; also checks that the neccessary files exist
for subject_i = 1:length(subject_list)
    
    subject_id = subject_list(subject_i);
    if parameters.subject_subfolder
        results_prefix = sprintf('sub-%1$03d/sub-%1$03d', subject_id);
    else
        results_prefix = sprintf('sub-%1$03d', subject_id);
    end
    fprintf('Subject %i, first pass\n', subject_id)
    
    headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
    isppa_map_mni_file  = fullfile(outputs_path, sprintf('%s_final_isppa_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    segmented_image_mni_file = fullfile(outputs_path, sprintf('%s_segmented_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    max_pressure_mni_file = fullfile(outputs_path, sprintf('%s_final_pressure_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    output_pressure_file = fullfile(outputs_path,sprintf('%s_%s_isppa%s.csv', results_prefix, parameters.simulation_medium, parameters.results_filename_affix));
    
    t1_mni_file = fullfile(headreco_folder, 'toMNI','T1fs_nu_12DOF_MNI.nii.gz');
    files_to_check = {t1_mni_file , segmented_image_mni_file, isppa_map_mni_file, max_pressure_mni_file, output_pressure_file};
    
    if options.plot_heating
        heating_data_mni_file = fullfile(outputs_path,sprintf('%s_heating_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
        files_to_check{length(files_to_check)+1} = heating_data_mni_file;    
    end

    for filename = files_to_check
        assert(logical(exist(filename{:}, 'file')), sprintf('File does not exist: %s; quitting', filename{:}))
    end

    % get T1 in MNI space
    t1_mni = niftiread(t1_mni_file);
    t1_mni_hdr = niftiinfo(t1_mni_file);
    
    Isppa_map_mni = niftiread(isppa_map_mni_file);
    segmented_image_mni = niftiread(segmented_image_mni_file);
    segmented_mask = (segmented_image_mni ==0 | segmented_image_mni ==4 | segmented_image_mni ==5);
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
    segmented_image_mni_file = fullfile(outputs_path, sprintf('%s_segmented_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    max_pressure_mni_file = fullfile(outputs_path, sprintf('%s_final_pressure_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
    t1_mni_file = fullfile(headreco_folder, 'toMNI','T1fs_nu_12DOF_MNI.nii.gz');
    output_pressure_file = fullfile(outputs_path,sprintf('%s_%s_isppa%s.csv', results_prefix, parameters.simulation_medium, parameters.results_filename_affix));

    % get T1 in MNI space
    t1_mni = niftiread(t1_mni_file);
    t1_mni_hdr = niftiinfo(t1_mni_file);
    
    Isppa_map_mni = niftiread(isppa_map_mni_file);
    segmented_image_mni = niftiread(segmented_image_mni_file);
    segmented_mask = (segmented_image_mni ==0 | segmented_image_mni ==4 | segmented_image_mni ==5);
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
        heating_data_mni_file = fullfile(outputs_path,sprintf('%s_heating_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
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
    curTable = readtable(output_pressure_file);
    curTable.('fwhm_size_MNI_based_on_pressure') = fwhm_size;

    slice = find(strcmp(slice_labels,options.slice_to_plot));
    cur_rotation = options.rotation;
        
    if isfield(options,'ROI_MNI_mask')
        roi_size = sum(options.ROI_MNI_mask,'all');
        avg_isppa_within_roi = mean(Isppa_map_mni(options.ROI_MNI_mask),'all');
        n_voxels_within_roi_above_thresh =  sum(options.ROI_MNI_mask & (max_pressure_map_mni >= max_pressure/2),'all');
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

    if isfield(options,'ROI_MNI_mask')
        hold on
        mask_im = get_slice_by_label(options.ROI_MNI_mask, options.slice_label, slice_n);
        visboundaries(imrotate(mask_im, cur_rotation))
    end

    if isfield(options,'add_FWHM_boundary')
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

suffix_list = {sprintf('maxT_MNI%s%s', parameters.results_filename_affix, options.outputs_suffix),...
    sprintf('final_isppa_MNI%s%s', parameters.results_filename_affix, options.outputs_suffix)};

for suffix_cell = suffix_list
    suffix = suffix_cell{:};
    
    if ~exist(sprintf('%stmp/', outputs_path),'dir')
        mkdir(sprintf('%stmp/', outputs_path));
    end
    system(sprintf('cd "%stmp/"; rm *.png', outputs_path));
    if parameters.subject_subfolder
        orig_imsize = size(imread(fullfile(outputs_path, sprintf('sub-%03i/sub-%03i_%s.png', subject_list(1),subject_list(1), suffix))));
        imarray = uint8(zeros(size(imread(fullfile(outputs_path, sprintf('sub-%03i/sub-%03i_%s.png', subject_list(end), subject_list(end), suffix))))));
    else
        orig_imsize = size(imread(fullfile(outputs_path, sprintf('sub-%03i_%s.png', subject_list(1),  suffix))));
        imarray = uint8(zeros(size(imread(fullfile(outputs_path, sprintf('sub-%03i_%s.png', subject_list(end), suffix))))));
    end

    imarray = repmat(imarray, 1, 1, 1, length(subject_list));
    for subject_i = 1:length(subject_list)
        if parameters.subject_subfolder
            old_name = fullfile(outputs_path, sprintf('sub-%03i/sub-%03i_%s.png', subject_list(subject_i), subject_list(subject_i), suffix));
        else
            old_name = fullfile(outputs_path, sprintf('sub-%03i_%s.png', subject_list(subject_i), suffix));
        end
        I = imread(old_name);
        I = padarray(I, max([0 0 0; orig_imsize-size(I)],[], 1), I(1), 'pre');

        I = insertText(I,[0 5],sprintf('P%02i',subject_list(subject_i)),'FontSize', options.font_size, 'BoxOpacity',0,'TextColor','white');
        %imarray(1:size(I, 1), 1:size(I, 2), 1:size(I, 3), subject_i) = I;
        new_name = sprintf('%stmp/%i.png', outputs_path, subject_list(subject_i));
        imwrite(I, new_name);
    end

    system(sprintf('cd "%stmp/"; montage -background black -gravity south  -mode concatenate $(ls -1 *.png | sort -g)  ../all_%s.png', convertCharsToStrings(outputs_path), suffix));
end

fprintf('Completed, see final images in %s\n', convertCharsToStrings(outputs_path))

end