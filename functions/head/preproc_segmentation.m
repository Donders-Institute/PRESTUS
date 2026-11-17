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

    fprintf('[preproc_segmentation] Subject: %03d | T1 pattern: %s | Seg folder: %s\n', ...
        parameters.subject_id, parameters.path.t1_pattern, ...
        fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id)));
    fprintf('[preproc_segmentation] %s\n', 'Checking inputs...');

    % Define paths to T1 and T2 images
    % Patterns may contain sprintf format specifiers (e.g. %03d) or shell-style
    % globs (e.g. *T1w.nii.gz). Expand via dir() after sprintf so both work.
    t1_glob = fullfile(parameters.path.anat, sprintf(parameters.path.t1_pattern, parameters.subject_id));
    d_t1 = dir(t1_glob);
    assert(~isempty(d_t1), 'No T1 file matched pattern: %s', t1_glob);
    filename_t1 = fullfile(d_t1(1).folder, d_t1(1).name);

    % Determine whether segmentation will actually need to run so we can decide
    % whether a missing T2 is a hard error or a soft warning.
    segmentation_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
    filename_segmented  = fullfile(segmentation_folder, 'final_tissues.nii.gz');
    seg_needed = confirm_overwriting(filename_segmented, parameters) && ...
        (~isfield(parameters.io,'overwrite_simnibs') || parameters.io.overwrite_simnibs || ...
         ~exist(filename_segmented,'file'));

    if isfield(parameters.path, 't2_pattern') && ~isempty(parameters.path.t2_pattern)
        t2_glob = fullfile(parameters.path.anat, sprintf(parameters.path.t2_pattern, parameters.subject_id));
        d_t2 = dir(t2_glob);
        if ~isempty(d_t2)
            filename_t2 = fullfile(d_t2(1).folder, d_t2(1).name);
        elseif seg_needed
            % T2 required to run charm — hard error.
            error('preproc_segmentation:t2Missing', ...
                'No T2 file matched pattern: %s\n(Required for SimNIBS segmentation.)', t2_glob);
        else
            % Segmentation already exists; T2 was probably deleted after charm completed.
            warning('preproc_segmentation:t2Missing', ...
                'T2 file not found: %s\nSegmentation exists — continuing without T2.', t2_glob);
            filename_t2 = '';
        end
    else
        filename_t2 = '';
    end

    % Validate existence of input files (T2 is optional).
    files_to_check = {filename_t1};
    if ~isempty(filename_t2)
        files_to_check{end+1} = filename_t2;
    end
    check_availability(files_to_check)

    %% SEGMENTATION USING SIMNIBS

    fprintf('[preproc_segmentation] %s\n', 'Starting segmentation...');

    % Run segmentation (if necessary)
    if seg_needed

        if parameters.simulation.interactive == 0 || confirmation_dlg('This will run SEGMENTATION WITH SIMNIBS that takes a long time. Are you sure?', 'Yes', 'No')
            segmentation_run(parameters.path.anat, parameters.subject_id, filename_t1, filename_t2, parameters);
            % On HPC, segmentation_run submits a job and returns before outputs exist.
            % Return early and let the next pipeline invocation handle post-processing.
            if ~exist(filename_segmented, 'file')
                return;
            end
        else
            return;
        end
    else
        fprintf('[preproc_segmentation] %s\n', 'Segmentation available...');
        % An existing segmentation is reused as is, also when another SimNIBS version or
        % other charm settings are selected now: say so (io.overwrite_simnibs re-runs charm)
        diffs = segmentation_check_existing(segmentation_folder, parameters);
        for i = 1:numel(diffs)
            warning('preproc_segmentation:segmentationDiffers', ...
                '%s %s; it is reused as is (set io.overwrite_simnibs: 1 to segment again).', ...
                segmentation_folder, diffs{i});
        end
    end

    %% SEGMENTATION PIPELINE (charm + pseudoCT share the same pipeline)

    parameters = simnibs_version(segmentation_folder, parameters);


    if parameters.pct.enabled == 1 && ~isempty(pct_ct_file(parameters))
        % A measured CT (path.ct_pattern) is read as is (PCT_LOAD_CT): no UTE, no pseudo-CT
        fprintf('[preproc_segmentation] Measured CT available: %s\n', pct_ct_file(parameters));
    elseif parameters.pct.enabled == 1
        if isempty(filename_t2)
            error(['pct.enabled = 1 requires a CT: set path.ct_pattern to a measured CT, or ' ...
                   'path.t2_pattern to a UTE image to generate a pseudo-CT.']);
        end
        pct_dir = find_pseudoct_dir(parameters);
        if ~isempty(pct_dir)
            fprintf('[preproc_segmentation] %s\n', 'pseudoCT available...');
            check_pct_sidecar(pct_dir, parameters);
        else
            fprintf('[preproc_segmentation] %s\n', 'Starting pseudoCT generation...');
            log_timer('start', 'pct');       % run_stats: stage_pct_s (part of the segmentation stage)
            pct_create_pCT_run(parameters);
            log_timer('stop', 'pct');
        end
    end

    %% MNI SPACE T1
    % charm does not produce a T1-in-MNI image, so create one post-hoc,
    % outside the m2m folder (GET_T1_MNI_FILE)
    if ~strcmp(parameters.simulation.medium, 'phantom')
        path_to_input_img  = fullfile(segmentation_folder, 'T1.nii.gz');
        path_to_output_img = get_t1_mni_file(parameters.path.seg, parameters.subject_id);
        if ~exist(path_to_output_img, 'file')
            out_dir = fileparts(path_to_output_img);
            if ~isfolder(out_dir); mkdir(out_dir); end
            mni_warp_method = 'simnibs';
            if isfield(parameters, 'analysis') && isfield(parameters.analysis, 'mni_warp_method')
                mni_warp_method = parameters.analysis.mni_warp_method;
            end
            convert_final_to_MNI_simnibs(path_to_input_img, segmentation_folder, path_to_output_img, parameters, 'method', mni_warp_method);
        end
    end

    %% DEFACED DERIVATIVE
    % Produces m2m_sub-NNN_defaced/ for sharing; original folder is unchanged.
    % A T1-in-MNI image inside m2m (earlier outputs) is defaced too; the
    % current one lives outside m2m (GET_T1_MNI_FILE) and is not copied.
    % Off by default — enable with parameters.segmentation.deface = 1.
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'deface') && parameters.segmentation.deface
        preproc_deface_m2m(segmentation_folder, parameters);
    end

end
