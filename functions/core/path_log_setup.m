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
%   parameters   - PRESTUS config; io.output_dir and io.log_file may be absent
%                  on entry and are created by this function
%   prestus_path - absolute path to the PRESTUS root directory
%
% Output:
%   parameters   - input extended with io.output_dir, io.cache_dir,
%                  io.debug_dir (and subdirs), io.log_file, io.filename_output_table
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
    toolboxesLoc = fullfile(prestus_path, 'toolboxes');
    allPaths = regexp(path, pathsep, 'Split');

    if ~any(ismember(fullfile(functionsLoc, 'helper'), allPaths))
        safe_addpath(functionsLoc);
        disp(['Adding ', functionsLoc, ' and subfolders']);
    end

    if ~any(ismember(toolboxesLoc, allPaths))
        safe_addpath(toolboxesLoc);
        disp(['Adding ', toolboxesLoc, ' and subfolders']);
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
        parameters.io.output_dir = get_output_dir(parameters);
        if ~isfolder(parameters.io.output_dir); mkdir(parameters.io.output_dir); end
    end

    % [LOCALITE OUTPUT] Make subfolder (if enabled) and check if directory exists
    if isfield(parameters.path, 'localite') && ...
            (isstring(parameters.path.localite) || ischar(parameters.path.localite)) && ...
            ~any(strcmp(parameters.path.localite, {"", ''}))
        parameters.path.localite = fullfile(parameters.path.localite, sprintf('sub-%03d', subject_id));
        if ~exist(parameters.path.localite); mkdir(parameters.path.localite); end
    end

    % Output subdirectories
    % ├── cache/         — regenerable intermediates (checkpoints, matrices, T1-space property maps)
    % └── debug/         — diagnostic artefacts written only when debug=1
    %     ├── preproc/   — head preprocessing (rotation, cropping, skull visualizations)
    %     ├── medium/    — medium mapping (grid-space property matrices, pCT)
    %     └── source/    — source/transducer setup (element distribution plots)
    if isfield(parameters.io, 'output_dir') && ~isempty(parameters.io.output_dir)
        out = parameters.io.output_dir;

        parameters.io.cache_dir        = fullfile(out, 'cache');
        parameters.io.debug_dir        = fullfile(out, 'debug');
        parameters.io.debug_dir_preproc = fullfile(out, 'debug', 'preproc');
        parameters.io.debug_dir_medium  = fullfile(out, 'debug', 'medium');
        parameters.io.debug_dir_source  = fullfile(out, 'debug', 'source');

        for d = {parameters.io.cache_dir, parameters.io.debug_dir, ...
                 parameters.io.debug_dir_preproc, parameters.io.debug_dir_medium, ...
                 parameters.io.debug_dir_source}
            if ~isfolder(d{1}); mkdir(d{1}); end
        end
    end

    % create segmentation folder if it doesn't exist
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg) && ~isfolder(parameters.path.seg)
        mkdir(parameters.path.seg);
    end

    % Resolve the subject-specific pCT output directory into parameters.io.pct_dir.
    % path.pct is a user-configured base directory (preserved as-is).
    %   empty    → io.pct_dir = {path.seg}/m2m_sub-NNN/  (uses existing SimNIBS m2m folder)
    %   non-empty → io.pct_dir = {path.pct}/sub-NNN/     (dedicated dir, PRESTUS naming)
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg)
        if ~isfield(parameters.path, 'pct') || isempty(parameters.path.pct)
            parameters.io.pct_dir = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', subject_id));
        else
            parameters.io.pct_dir = fullfile(parameters.path.pct, sprintf('sub-%03d', subject_id));
        end
        % Only create the directory when pCT is actually enabled; an unconditional
        % mkdir here (after cd(prestus_path) above) would create sub-NNN in the
        % PRESTUS root if path.pct is a relative path.
        pct_enabled = isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && ...
                      parameters.pct.enabled;
        if pct_enabled && ~isfolder(parameters.io.pct_dir)
            mkdir(parameters.io.pct_dir);
        end
    end
    
    if isfield(parameters.io, 'output_dir') && ~isempty(parameters.io.output_dir)
        % Save parameter snapshot to cache (regenerable, not a primary output)
        filename_parameters = fullfile(parameters.io.cache_dir, ...
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
            filename_log = fullfile(parameters.io.output_dir, ...
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
    if isfield(parameters.io, 'output_dir') && ~isempty(parameters.io.output_dir)
        parameters.io.filename_output_table = ...
            fullfile(parameters.io.output_dir,sprintf('sub-%03d_%s_output_table%s.csv', ...
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
    log_timer('start','prestus_pipeline', parameters.io.output_dir);