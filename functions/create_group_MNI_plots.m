function create_group_MNI_plots(subject_list, parameters, options)

% CREATE_GROUP_MNI_PLOTS Generates group-level plots in MNI space for multiple subjects.
%
% This function processes subject-specific data to create plots in MNI space, 
% including  ISPPA maps, pressure maps, and optional heating data. It handles missing 
% files, applies masks, computes statistics within regions of interest (ROIs), and 
% combines individual subject plots into a single montage image.
%
% Input:
%   subject_list - List of subject IDs to process.
%   parameters   - Struct containing pipeline configuration parameters:
%                  * temp_output_dir: Path to temporary output directory.
%                  * layer_labels: Labels for tissue layers (e.g., brain, water).
%                  * segmentation_software: Software used for segmentation (e.g., 'headreco').
%   options      - Struct containing additional options:
%                  * ROI_MNI_mask: Mask defining the region of interest in MNI space.
%                  * slice_to_plot: Slice number to plot (default: 0).
%                  * plot_max_intensity: Flag to plot maximum intensity slice (default: 0).
%                  * slice_label: Axis label for slicing ('x', 'y', or 'z', default: 'y').
%                  * rotation: Rotation angle for plots (default: 90 degrees).
%                  * plot_heating: Flag to include heating data in plots (default: 1).
%                  * outputs_suffix: Suffix for output filenames (default: '').
%                  * isppa_thresholds: Thresholds for  ISPPA values (default: []).
%                  * add_FWHM_boundary: Flag to add full-width half-maximum boundary (default: 0).
%                  * add_ROI_boundary: Flag to add ROI boundary (default: 1).
%                  * skip_missing: Flag to skip subjects with missing files (default: 0).
%                  * brightness_correction: Flag to apply brightness correction (default: 0).
%                  * average_target_brightness: Target brightness value for correction (default: 100).
%
% Output:
%   Plots are saved in the specified output directory and combined into montage images.

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
        options.brightness_correction = 0
        options.average_target_brightness = 100
    end

    % Validate input arguments
    assert(xor(options.plot_max_intensity, options.slice_to_plot), ...
           "You should indicate either the slice number to plot ('slice_to_plot') or ask for a max intensity plot ('plot_max_intensity')");

    % Initialize variables and ranges for colorbars
    slice_labels = {'x', 'y', 'z'};
    outputs_path = parameters.temp_output_dir;
    bg_range_to_use = [];
    isppa_range_to_use = [5, 6];
    
    % Set temperature ranges based on tissue-specific or general values
    if isstruct(parameters.thermal.temp_0)
        temp_range_to_use = [parameters.thermal.temp_0.skin, parameters.thermal.temp_0.water + 0.5];
        temp_0_min = parameters.thermal.temp_0.skin;
    else
        temp_range_to_use = [parameters.thermal.temp_0, parameters.thermal.temp_0 + 0.5];
        temp_0_min = parameters.thermal.temp_0;
    end

    full_subject_list = subject_list;

    % First pass: Gather background intensity and  ISPPA range information
    for subject_i = 1:length(full_subject_list)
        
        subject_id = full_subject_list(subject_i);
        
        % Determine file paths based on subject ID and segmentation software
        if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
            results_prefix = sprintf('sub-%1$03d/sub-%1$03d', subject_id);
        else
            results_prefix = sprintf('sub-%1$03d', subject_id);
        end
        
        fprintf('Subject %i, first pass\n', subject_id);
        
        if isfield(parameters,'seg_path') && ~isempty(parameters.seg_path)
            headreco_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', subject_id));
        else
            headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
        end
        
        % Define file paths for required data files
        isppa_map_mni_file  = fullfile(outputs_path, sprintf('%s_final_isppa_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
        segmented_image_mni_file = fullfile(outputs_path, sprintf('%s_final_medium_masks_MNI%s.nii.gz', results_prefix, parameters.results_filename_affix));
        
        % Check file existence and handle missing files based on user preferences
        files_to_check = {isppa_map_mni_file, segmented_image_mni_file};
        
        all_files_exist = true;
        
        for filename = files_to_check
            if ~exist(filename{:}, 'file')
                if options.skip_missing
                    fprintf('File does not exist: %s; skipping this subject\n', filename{:});
                    subject_list(subject_i) = NaN;
                    all_files_exist = false;
                    break;
                else
                    error(sprintf('File does not exist: %s; quitting\n', filename{:}));
                end
            end
        end
        
        if ~all_files_exist 
            continue;
        end

        % Load required data files and apply masks to  ISPPA map and pressure map
        t1_mni = niftiread(isppa_map_mni_file);
        
    end

    % Remove unavailable subjects from the list
    subject_list(isnan(subject_list)) = [];

    % Second pass: Generate images for each subject and combine them into montages
    for suffix_cell in suffix_list
       combine_plots_by_suffix(suffix_cell{:}, outputs_path, subject_list, parameters);
    end

    fprintf('Completed! Final images are saved in %s\n', convertCharsToStrings(outputs_path));

end
