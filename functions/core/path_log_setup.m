function [parameters] = path_log_setup(parameters, prestus_path)
% PATH_LOG_SETUP  Add PRESTUS paths, create output directories, and start diary logging
%
% Adds PRESTUS functions and toolbox folders to the MATLAB path (filtering hidden
% directories via safe_addpath), creates the simulation output directory tree
% (output, cache, debug/preproc, debug/medium, debug/source), saves a parameter
% snapshot to cache, opens a diary log file, and writes an initial parameter
% summary. Also verifies k-Wave is on the path and prints version information.
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
%                  io.dir_debug (and subdirs), io.log_file, io.filename_table
%
% See also: LOAD_PARAMETERS, PRINT_PARAMETER_SUMMARY, PRESTUS_PIPELINE_START

arguments
    parameters   (1,1) struct
    prestus_path (1,:) char
end

    % currentLoc | Root of PRESTUS' 'functions' folder

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
            disp(['Adding ', parameters.startup.subpaths_to_add{nPaths}, 'and subfolders']);
        end
    end
    
    % Print current PRESTUS version (Git-based)
    prestus_version(prestus_path);

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

    % [LOCALITE OUTPUT] Make subfolder (if enabled) and check if directory exists
    if isfield(parameters.path, 'localite') && ...
            (isstring(parameters.path.localite) || ischar(parameters.path.localite)) && ...
            ~any(strcmp(parameters.path.localite, {"", ''}))
        parameters.path.localite = fullfile(parameters.path.localite, sprintf('sub-%03d', subject_id));
        if ~exist(parameters.path.localite); mkdir(parameters.path.localite); end
    end

    % Output subdirectories
    % ├── nii/     — NIfTIs (T1w and MNI, differentiated by filename entity)
    % ├── img/     — all PNG visualizations (acoustic, thermal, preprocessing)
    % ├── log/     — diary log files
    % ├── cache/   — regenerable intermediates (checkpoints, matrices)
    % └── debug/   — diagnostic artefacts written only when debug=1
    %     ├── preproc/  — head preprocessing (rotation, cropping, skull visualizations)
    %     ├── medium/   — medium mapping (grid-space property matrices, pCT)
    %     └── source/   — source/transducer setup (element distribution plots)
    % Reports (.html) and tables (.csv) are written directly to dir_output.
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        out = parameters.io.dir_output;

        parameters.io.dir_nii              = fullfile(out, 'nii');
        % nii_T1w_dir and nii_MNI_dir both point to nii/; space is encoded in filenames
        parameters.io.dir_nii_T1w          = fullfile(out, 'nii');
        parameters.io.dir_nii_MNI          = fullfile(out, 'nii');
        parameters.io.dir_img              = fullfile(out, 'img');
        parameters.io.dir_tabular          = out;
        parameters.io.dir_reports          = out;
        parameters.io.dir_logs             = fullfile(out, 'log');
        parameters.io.dir_cache            = fullfile(out, 'cache');
        parameters.io.dir_debug            = fullfile(out, 'debug');
        parameters.io.dir_debug_preproc    = fullfile(out, 'debug', 'preproc');
        parameters.io.dir_debug_medium     = fullfile(out, 'debug', 'medium');
        parameters.io.dir_debug_source     = fullfile(out, 'debug', 'source');

        for d = {parameters.io.dir_nii, parameters.io.dir_img, ...
                 parameters.io.dir_logs, parameters.io.dir_cache}
            if ~isfolder(d{1}); mkdir(d{1}); end
        end

        if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'debug') && ...
                parameters.simulation.debug == 1
            for d = {parameters.io.dir_debug, parameters.io.dir_debug_preproc, ...
                     parameters.io.dir_debug_medium, parameters.io.dir_debug_source}
                if ~isfolder(d{1}); mkdir(d{1}); end
            end
        end
    end

    % create segmentation folder if it doesn't exist
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg) && ~isfolder(parameters.path.seg)
        mkdir(parameters.path.seg);
    end

    % Resolve the subject-specific pCT output directory into parameters.io.dir_pct.
    % path.pct is a user-configured base directory (preserved as-is).
    %   empty    → io.dir_pct = {path.seg}/m2m_sub-NNN/  (uses existing SimNIBS m2m folder)
    %   non-empty → io.dir_pct = {path.pct}/sub-NNN/     (dedicated dir, PRESTUS naming)
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg)
        if ~isfield(parameters.path, 'pct') || isempty(parameters.path.pct)
            parameters.io.dir_pct = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', subject_id));
        else
            parameters.io.dir_pct = fullfile(parameters.path.pct, sprintf('sub-%03d', subject_id));
        end
        % Only create the directory when pCT is actually enabled; an unconditional
        % mkdir here (after cd(prestus_path) above) would create sub-NNN in the
        % PRESTUS root if path.pct is a relative path.
        pct_enabled = isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && ...
                      parameters.pct.enabled;
        if pct_enabled && ~isfolder(parameters.io.dir_pct)
            mkdir(parameters.io.dir_pct);
        end
    end
    
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        % Save parameter snapshot to cache (regenerable, not a primary output)
        filename_parameters = fullfile(parameters.io.dir_cache, ...
            sprintf('sub-%03d_%s%s_parameters_%s.mat', ...
            subject_id, parameters.simulation.medium, parameters.io.output_affix, ...
            string(datetime('now'), 'yyMMdd_HHmm')));
        save(filename_parameters, 'parameters');
        clear filename_parameters;

        % Create a log.
        % Use a pre-assigned path if uncertainty_pipeline set one (so all
        % five stage logs have known, deterministic paths); otherwise fall
        % back to a timestamp-based name.
        if isfield(parameters.io, 'log_file') && ~isempty(parameters.io.log_file)
            filename_log = parameters.io.log_file;
        else
            filename_log = fullfile(parameters.io.dir_logs, ...
                sprintf('sub-%03d_%s%s_%s.txt', ...
                subject_id, parameters.simulation.medium, parameters.io.output_affix, ...
                string(datetime('now'), 'yyMMdd_HHmm')));
        end
        parameters.io.log_file = filename_log;
        diary(filename_log);
    end

    % summarize parameters
    print_parameter_summary(parameters)

    % Define the filename of the summary table
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        parameters.io.filename_table = ...
            fullfile(parameters.io.dir_tabular, sprintf('sub-%03d_%s%s.csv', ...
            subject_id, parameters.simulation.medium, parameters.io.output_affix));
    end

    % suppress unneccessary warnings from export_fig when running without OpenGL
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