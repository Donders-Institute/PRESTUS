function [parameters] = path_log_setup(parameters, prestus_path)
% PATH_LOG_SETUP  Add PRESTUS paths, create output directories, and start diary logging
%
% Adds PRESTUS functions and toolbox folders to the MATLAB path (filtering hidden
% directories via safe_addpath), creates the simulation output directory tree
% (output, cache, debug/preproc, debug/medium, debug/source), saves a parameter
% snapshot to cache, opens a diary log file, and writes an initial parameter
% summary. Also verifies k-Wave is on the path and prints version information.
% Relative file and folder settings are made absolute against the PRESTUS root
% (RESOLVE_CONFIG_PATHS); MATLAB's working folder is not changed.
%
% Use as:
%   parameters = path_log_setup(parameters, prestus_path)
%
% Input:
%   parameters   - PRESTUS config; io.dir_output and io.log_file may be absent
%                  on entry and are created by this function
%   prestus_path - absolute path to the PRESTUS root directory
%
% Output:
%   parameters   - input extended with io.dir_output, io.dir_cache,
%                  io.dir_debug (and subdirs), io.log_file, io.filename_table,
%                  io.prestus_version (PRESTUS_VERSION provenance struct)
%
% See also: LOAD_PARAMETERS, PRINT_PARAMETER_SUMMARY, PRESTUS_PIPELINE_START

arguments
    parameters   (1,1) struct
    prestus_path (1,:) char
end

    disp(['Location of PRESTUS: ', prestus_path])

    % Add paths to the 'functions' and toolbox folders.
    % safe_addpath filters out hidden directories (e.g. .claude, .git) to
    % prevent stale file copies from shadowing the current versions.
    functionsLoc = fullfile(prestus_path, 'functions');
    externalLoc = fullfile(prestus_path, 'external');
    allPaths = regexp(path, pathsep, 'Split');

    if ~any(ismember(fullfile(functionsLoc, 'helper'), allPaths))
        safe_addpath(functionsLoc);
        disp(['Adding ', functionsLoc, ' and subfolders']);
    end

    if ~any(ismember(externalLoc, allPaths))
        safe_addpath(externalLoc);
        disp(['Adding ', externalLoc, ' and subfolders']);
    end

    % If there are paths to be added, add them; this is mostly for batch runs
    if isfield(parameters.startup, 'paths_to_add') && ~isempty(parameters.startup.paths_to_add)
        for nPaths = 1:length(parameters.startup.paths_to_add)
            addpath(parameters.startup.paths_to_add{nPaths})
            disp(['Adding ', parameters.startup.paths_to_add{nPaths}]);
        end
    end

    % If the path and subpaths need to be added, use this instead
    if isfield(parameters.startup, 'subpaths_to_add') && ~isempty(parameters.startup.subpaths_to_add)
        for nPaths = 1:length(parameters.startup.subpaths_to_add)
            addpath(genpath(parameters.startup.subpaths_to_add{nPaths}))
            disp(['Adding ', parameters.startup.subpaths_to_add{nPaths}, ' and subfolders']);
        end
    end
    
    % Print current PRESTUS version (Git-based) and keep it with the run:
    % io.prestus_version goes into the parameter snapshots (.mat, .yaml),
    % from which result summaries read which code produced the outputs.
    prestus_version(prestus_path);
    [~, ~, parameters.io.prestus_version] = prestus_version(prestus_path);

    % Verify that kwave is added and print k-wave information
    if ~exist('makeBowl','file')
        error('kwave not added');
    else
        kwave_version(getkWavePath);
    end

    % return to PRESTUS path
    cd(prestus_path);

    subject_id = parameters.subject_id;

    % [SIMULATION OUTPUT] Make subfolder (if enabled) and check if directory exists
    if isfield(parameters.path, 'sim') && ~isempty(parameters.path.sim) && ...
            (isstring(parameters.path.sim) || ischar(parameters.path.sim)) && ...
            ~any(strcmp(parameters.path.sim, {"", ''}))
        parameters.io.dir_output = get_output_dir(parameters);
        if ~isfolder(parameters.io.dir_output); mkdir(parameters.io.dir_output); end
    end

    % [LOCALITE EXPORT] Per-subject folder for files meant for import into
    % Localite (heuristic placement: position table, plot, Localite-ready T1).
    % path.localite_export; path.localite is the deprecated fallback. The
    % path.* values are left untouched (path.localite may also be the
    % deprecated alias of path.localite_raw, see neuronav_resolve_paths).
    parameters.io.dir_localite_export = '';
    export_root = '';
    for key = {'localite_export', 'localite'}
        if isfield(parameters.path, key{1}) && ...
                (isstring(parameters.path.(key{1})) || ischar(parameters.path.(key{1}))) && ...
                strlength(string(parameters.path.(key{1}))) > 0
            export_root = char(parameters.path.(key{1}));
            break;
        end
    end
    if ~isempty(export_root)
        parameters.io.dir_localite_export = fullfile(export_root, sprintf('sub-%03d', subject_id));
        if ~isfolder(parameters.io.dir_localite_export); mkdir(parameters.io.dir_localite_export); end
    end

    % Output subdirectories
    % ├── nii/     — NIfTIs (T1w and MNI, differentiated by filename entity)
    % ├── img/     — all PNG visualizations (acoustic, thermal, preprocessing)
    % ├── log/     — diary log files
    % ├── cache/   — regenerable intermediates (checkpoints, matrices)
    % └── debug/   — diagnostic artefacts written only when simulation.debug=1
    %     ├── manifest.jsonl — one line per debug output (DEBUG_REGISTER)
    %     ├── preproc/   — head preprocessing (rotation, cropping, skull visualizations)
    %     ├── medium/    — medium mapping (grid-space property matrices, pCT)
    %     ├── source/    — source/transducer setup (element distribution plots)
    %     └── babelbrain/— BabelBrain FOV/domain-coverage diagnostics
    % Reports (.html) and tables (.csv) are written directly to dir_output.
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        out = char(parameters.io.dir_output);   % all io.dir_* are char (see GET_OUTPUT_DIR)
        parameters.io.dir_output = out;

        parameters.io.dir_nii              = fullfile(out, 'nii');
        % nii_T1w_dir and nii_MNI_dir both point to nii/; space is encoded in filenames
        parameters.io.dir_nii_T1w          = fullfile(out, 'nii');
        parameters.io.dir_nii_MNI          = fullfile(out, 'nii');
        parameters.io.dir_img              = fullfile(out, 'img');
        parameters.io.dir_tabular          = out;
        parameters.io.dir_reports          = out;
        parameters.io.dir_logs             = fullfile(out, 'log');
        parameters.io.dir_cache            = fullfile(out, 'cache');
        parameters.io.dir_babelbrain       = fullfile(out, 'babelbrain');
        parameters.io.dir_debug            = fullfile(out, 'debug');
        parameters.io.dir_debug_preproc    = fullfile(out, 'debug', 'preproc');
        parameters.io.dir_debug_medium     = fullfile(out, 'debug', 'medium');
        parameters.io.dir_debug_source     = fullfile(out, 'debug', 'source');
        parameters.io.dir_debug_babelbrain = fullfile(out, 'debug', 'babelbrain');

        nii_enabled = ~(isfield(parameters, 'modules') && ...
                        isfield(parameters.modules, 'run_nifti_creation') && ...
                        parameters.modules.run_nifti_creation == 0);
        always_create = {parameters.io.dir_img, parameters.io.dir_logs, parameters.io.dir_cache};
        if nii_enabled; always_create{end+1} = parameters.io.dir_nii; end
        for d = always_create
            if ~isfolder(d{1}); mkdir(d{1}); end
        end

        if debug_on(parameters)
            % Folders of the stages selected by simulation.debug_stages
            % (DEBUG_STAGES errors on an unknown stage name)
            stages = debug_stages(parameters);
            dirs = [{parameters.io.dir_debug}, ...
                cellfun(@(s) parameters.io.(['dir_debug_' s]), stages, 'UniformOutput', false)];
            for d = dirs
                if ~isfolder(d{1}); mkdir(d{1}); end
            end
        end
    end

    % create segmentation folder if it doesn't exist
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg) && ~isfolder(parameters.path.seg)
        mkdir(parameters.path.seg);
    end

    % Resolve the subject-specific pCT output directory into parameters.io.dir_pct
    % (GET_PCT_DIR): {path.pct}/sub-NNN, or {parent of path.seg}/pct/sub-NNN when
    % path.pct is empty. Never the SimNIBS m2m folder.
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg)
        parameters.io.dir_pct = get_pct_dir(parameters);
        % Only create the directory when pCT is actually enabled.
        pct_enabled = isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && ...
                      parameters.pct.enabled;
        measured_ct = strlength(string(get_field(parameters, {'path', 'ct_pattern'}, ''))) > 0;
        if pct_enabled && ~measured_ct && ~isfolder(parameters.io.dir_pct)
            mkdir(parameters.io.dir_pct);
        end
    end
    
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        % Save parameter snapshot to cache (regenerable, not a primary output)
        filename_parameters = fullfile(parameters.io.dir_cache, ...
            sprintf('sub-%03d_%s_parameters_%s%s.mat', ...
            subject_id, parameters.simulation.medium, ...
            string(datetime('now'), 'yyMMdd_HHmm'), parameters.io.output_affix));
        save(filename_parameters, 'parameters');
        clear filename_parameters;

        % Write human-readable resolved-parameters YAML to log folder so that
        % results are reproducible even if the transducer library is updated.
        try
            filename_params_yaml = fullfile(parameters.io.dir_logs, ...
                sprintf('sub-%03d_%s_parameters%s.yaml', ...
                subject_id, parameters.simulation.medium, parameters.io.output_affix));
            yaml.dumpFile(filename_params_yaml, parameters);
        catch yaml_err
            warning('prestus:yamlSnapshot', ...
                'Could not write parameter snapshot YAML: %s', yaml_err.message);
        end

        % Create a log.
        % Use a pre-assigned path if one was set -- by uncertainty_pipeline
        % (so every stage log has a known, deterministic path) or at HPC
        % submission (so hpc_job_report can parse it); otherwise fall back to
        % a timestamp-based name.
        if isfield(parameters.io, 'log_file') && ~isempty(parameters.io.log_file)
            filename_log = parameters.io.log_file;
        else
            filename_log = default_log_file(parameters, parameters.io.dir_logs);
        end
        parameters.io.log_file = filename_log;
        diary(filename_log);
        v = parameters.io.prestus_version;
        fprintf('PRESTUS %s (commit %s, %s; uncommitted changes: %s)\n', ...
            v.version, v.commit, v.commit_date, v.dirty);
    end
    for k = 1:numel(resolved_paths)
        fprintf('[path_log_setup] relative setting taken from the PRESTUS root: %s\n', resolved_paths{k});
    end

    % summarize parameters
    print_parameter_summary(parameters)

    % Define the filename of the summary table
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        parameters.io.filename_table = ...
            fullfile(parameters.io.dir_tabular, sprintf('sub-%03d_%s%s.csv', ...
            subject_id, parameters.simulation.medium, parameters.io.output_affix));
    end

    % suppress unnecessary warnings from export_fig when running without OpenGL
    warning('off','MATLAB:prnRenderer:opengl');
    
    % display GPU information (if requested)
    if strcmp(parameters.simulation.code_type, 'cpp_gpu') || strcmp(parameters.simulation.code_type, 'matlab_gpu')
        fprintf('========================================\n');
        fprintf('GPU INFO \n');
        fprintf('========================================\n\n');
        gpuDevice()
        fprintf('========================================\n\n');
    end

    % set initial time, RAM, GB state
    log_timer('start','prestus_pipeline', parameters.io.dir_output);