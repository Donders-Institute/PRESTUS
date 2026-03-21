function parameters = preproc_segmentation(parameters)

% PREPROC_SEGMENTATION Segment the data using SimNIBS
%
% Note: T1, T2 (or UTE/CT) paths should correspond to BIDS directory

    %% CHECK INPUTS

    disp('Checking inputs...');
    
    % Define paths to T1 and T2 images
    filename_t1 = fullfile(parameters.path.anat, sprintf(parameters.path.t1_pattern, parameters.subject_id));
    filename_t2 = fullfile(parameters.path.anat, sprintf(parameters.path.t2_pattern, parameters.subject_id));

    % Validate existence of files
    files_to_check = {filename_t1, filename_t2};
    check_availability(files_to_check)
    
    %% SEGMENTATION USING SIMNIBS
    
    disp('Starting segmentation...');
    
    % Define output folder for segmentation results
    segmentation_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

    filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');

    % Run segmentation (if necessary)
    if confirm_overwriting(filename_segmented, parameters) && ...
       (~isfield(parameters.io,'overwrite_simnibs') || parameters.io.overwrite_simnibs || ~exist(filename_segmented,'file'))
        if parameters.pct.enabled == 1
            % Note: This could be improved by allowing to specify a UTE/CT path in the config...
            warning("SimNIBS integration not supported when requesting pseudoCT. Please ensure SimNIBS has been run.");
        end

        if parameters.simulation.interactive == 0 || confirmation_dlg('This will run SEGMENTATION WITH SIMNIBS that takes a long time. Are you sure?', 'Yes', 'No')
            segmentation_run(parameters.path.anat, parameters.subject_id, filename_t1, filename_t2, parameters);
            parameters = simnibs_version(segmentation_folder, parameters);
            return;
        end
    else
        disp('Segmentation available...');
        parameters = simnibs_version(segmentation_folder, parameters);
    end
