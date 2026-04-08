function [parameters] = path_log_setup(parameters, prestus_path)

    % currentLoc | Root of PRESTUS' 'functions' folder

    disp(['Location of PRESTUS: ', prestus_path])

    % Add paths to the 'functions' and toolbox folders
    functionsLoc = fullfile(prestus_path, 'functions');
    toolboxesLoc = fullfile(prestus_path, 'toolboxes');
    allPaths = regexp(path,pathsep,'Split');

    % Add 'functions' folders
    if ~any(ismember(fullfile(functionsLoc, 'helper'),allPaths))
        addpath(genpath(functionsLoc));
        disp(['Adding ', functionsLoc, ' and subfolders']);
    else
    end

    % Add 'toolboxes' 
    if ~any(ismember(toolboxesLoc,allPaths))
        addpath(genpath(toolboxesLoc));
        disp(['Adding ', toolboxesLoc, ' and subfolders']);
    else
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
    if isfield(parameters.path, 'sim') && isstring(parameters.path.sim) && ~strcmp(parameters.path.sim, "")
        if isfield(parameters.path, 'subject_subfolder') && parameters.path.subject_subfolder == 1
            parameters.io.output_dir = fullfile(parameters.path.sim, sprintf('sub-%03d', subject_id));
        else
            parameters.io.output_dir = parameters.path.sim;
        end
        if ~isfolder(parameters.io.output_dir); mkdir(parameters.io.output_dir); end
    end

    % [LOCALITE OUTPUT] Make subfolder (if enabled) and check if directory exists
    if isfield(parameters.path, 'localite') && isstring(parameters.path.sim) && ~strcmp(parameters.path.localite, "")
        if isfield(parameters.path, 'subject_subfolder') && parameters.path.subject_subfolder == 1
            parameters.path.localite = fullfile(parameters.path.localite, sprintf('sub-%03d', subject_id));
        end
        if ~exist(parameters.path.localite); mkdir(parameters.path.localite); end
    end

    % specify dedicated subfolder for debugging contents
    if isfield(parameters.io, 'output_dir') && ~isempty(parameters.io.output_dir)
        parameters.io.debug_dir = fullfile(parameters.io.output_dir, 'debug');
        if ~isfolder(parameters.io.debug_dir); mkdir(parameters.io.debug_dir); end
    end

    % create segmentation folder if it doesn't exist
    if isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg) && ~isfolder(parameters.path.seg)
        mkdir(parameters.path.seg);
    end
    
    if isfield(parameters.io, 'output_dir') && ~isempty(parameters.io.output_dir)
        % Save parameters
        filename_parameters = fullfile(parameters.io.output_dir, ...
            sprintf('sub-%03d_parameters_%s%s_%s.mat', ...
            subject_id, parameters.simulation.medium, parameters.io.output_affix, ...
            string(datetime('now'), 'yyMMdd_HHmm')));
        save(filename_parameters, 'parameters');
        clear filename_parameters;

        % Create a log
        filename_log = fullfile(parameters.io.output_dir, ...
            sprintf('sub-%03d_%s%s_%s.txt', ...
            subject_id, parameters.simulation.medium, parameters.io.output_affix, ...
            string(datetime('now'), 'yyMMdd_HHmm')));
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