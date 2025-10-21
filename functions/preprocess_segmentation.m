function preprocess_segmentation(parameters, subject_id)

% PREPROCESS_SEGMENTATION Segment the data using SimNIBS
%
% Note: T1, T2 (or UTE/CT) paths should correspond to BIDS directory

    %% CHECK INPUTS

    disp('Checking inputs...');
    
    % Define paths to T1 and T2 images
    filename_t1 = fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id));
    filename_t2 = fullfile(parameters.data_path, sprintf(parameters.t2_path_template, subject_id));

    % Validate existence of files
    files_to_check = {filename_t1, filename_t2};
    check_availability(files_to_check)
    
    %% SEGMENTATION USING SIMNIBS
    
    disp('Starting segmentation...');
    
    % Define output folder for segmentation results
    segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', subject_id));

    if strcmp(parameters.segmentation_software, 'charm')
        filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');
    else 
        filename_segmented = fullfile(segmentation_folder, sprintf('sub-%03d_final_contr.nii.gz', subject_id));
    end

    % Run segmentation (if necessary)
    if confirm_overwriting(filename_segmented, parameters) && ...
       (~isfield(parameters,'overwrite_simnibs') || parameters.overwrite_simnibs || ~exist(filename_segmented,'file'))
        if parameters.usepseudoCT == 1
            % Note: This could be improved by allowing to specify a UTE/CT path in the config...
            warning("SimNIBS integration not supported when requesting pseudoCT. Please ensure SimNIBS has been run.");
        end

        if parameters.interactive == 0 || confirmation_dlg('This will run SEGMENTATION WITH SIMNIBS that takes a long time. Are you sure?', 'Yes', 'No')
            run_segmentation(parameters.data_path, subject_id, filename_t1, filename_t2, parameters);
            return;
        end
    else
        disp('Skipping segmentation; loading existing file instead.');
    end
