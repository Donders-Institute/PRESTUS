function parameters = preproc_segmentation(parameters)

% PREPROC_SEGMENTATION  Segment structural MRI using SimNIBS
%
% Locates T1 (and optionally T2/UTE/CT) images following BIDS naming,
% runs SimNIBS charm segmentation, and saves the output. Paths in
% parameters must follow BIDS directory conventions.
%
% Use as:
%   parameters = preproc_segmentation(parameters)
%
% Input:
%   parameters - (1,1) simulation configuration struct with path fields
%
% Output:
%   parameters - updated struct with segmentation path fields
%
% See also: PREPROC_HEAD, SEGMENTATION_RUN

arguments
    parameters (1,1) struct
end

    %% CHECK INPUTS

    disp('Checking inputs...');
    
    % Define paths to T1 and T2 images
    filename_t1 = fullfile(parameters.path.anat, sprintf(parameters.path.t1_pattern, parameters.subject_id));
    if isfield(parameters.path, 't2_pattern') && ~isempty(parameters.path.t2_pattern)
        filename_t2 = fullfile(parameters.path.anat, sprintf(parameters.path.t2_pattern, parameters.subject_id));
    else
        filename_t2 = '';
    end

    % Validate existence of files (T2 is optional)
    files_to_check = {filename_t1};
    if ~isempty(filename_t2)
        files_to_check{end+1} = filename_t2;
    end
    check_availability(files_to_check)
    
    %% SEGMENTATION USING SIMNIBS
    
    disp('Starting segmentation...');
    
    % Define output folder for segmentation results
    segmentation_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

    filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');

    % Run segmentation (if necessary)
    if confirm_overwriting(filename_segmented, parameters) && ...
       (~isfield(parameters.io,'overwrite_simnibs') || parameters.io.overwrite_simnibs || ~exist(filename_segmented,'file'))

        if parameters.simulation.interactive == 0 || confirmation_dlg('This will run SEGMENTATION WITH SIMNIBS that takes a long time. Are you sure?', 'Yes', 'No')
            segmentation_run(parameters.path.anat, parameters.subject_id, filename_t1, filename_t2, parameters);
            parameters = simnibs_version(segmentation_folder, parameters);
            if parameters.pct.enabled == 1
                if isempty(filename_t2)
                    error('pct.enabled = 1 requires a UTE image. Set path.t2_pattern to the UTE file.');
                end
                disp('Starting pseudoCT generation...');
                pct_create_pCT_run(parameters);
            end
            return;
        end
    else
        disp('Segmentation available...');
        parameters = simnibs_version(segmentation_folder, parameters);
        if parameters.pct.enabled == 1
            filename_pseudoCT = fullfile(segmentation_folder, 'pseudoCT.nii.gz');
            if exist(filename_pseudoCT, 'file')
                disp('pseudoCT available...');
            else
                if isempty(filename_t2)
                    error('pct.enabled = 1 requires a UTE image. Set path.t2_pattern to the UTE file.');
                end
                disp('pseudoCT not found — generating from UTE image...');
                pct_create_pCT_run(parameters);
            end
        end
    end

end
