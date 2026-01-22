function [parameters] = path_log_setup(parameters, currentLoc, subject_id)

    % currentLoc | Root of PRESTUS' 'functions' folder
    
    disp(['Location of PRESTUS: ', currentLoc])

    % Add paths to the 'functions' and toolbox folders
    functionsLoc = fullfile(currentLoc, '..', 'functions');
    toolboxesLoc = fullfile(currentLoc, '..', 'toolboxes');
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
    if isfield(parameters,'paths_to_add') && ~isempty(parameters.paths_to_add)
        for nPaths = 1:length(parameters.paths_to_add)
            addpath(parameters.paths_to_add{nPaths})
            disp(['Adding ', parameters.paths_to_add{nPaths}]);
        end
    end

    % If the path and subpaths need to be added, use this instead
    if isfield(parameters,'subpaths_to_add') && ~isempty(parameters.subpaths_to_add)
        for nPaths = 1:length(parameters.subpaths_to_add)
            addpath(genpath(parameters.subpaths_to_add{nPaths}))
            disp(['Adding ', parameters.subpaths_to_add{nPaths}, 'and subfolders']);
        end
    end
    
    % Verify that kwave is added and print k-wave information
    if ~exist('makeBowl','file')
        error('kwave not added');
    else
        kwave_version(getkWavePath);
    end

    % Make subfolder (if enabled) and check if directory exists
    if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
        parameters.output_dir = fullfile(parameters.sim_path, sprintf('sub-%03d', subject_id));
    else 
        parameters.output_dir = parameters.sim_path;
    end

    % specify dedicated subfolder for debugging contents
    parameters.debug_dir = fullfile(parameters.output_dir, 'debug');

    if ~isfolder(parameters.output_dir)
        mkdir(parameters.output_dir);
    end
    if ~isfolder(parameters.debug_dir)
        mkdir(parameters.debug_dir);
    end
    if isfield(parameters,'seg_path') && ~isfolder(parameters.seg_path)
        mkdir(parameters.seg_path);
    end
    
    % Save parameters
    filename_parameters = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_parameters%s_%s_%s.mat', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix, ...
        string(datetime('now'), 'yyMMdd_HHmm')));
    save(filename_parameters, 'parameters');
    clear filename_parameters;

    % Create a log
    filename_log = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d%s_%s_%s.txt', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix, ...
        string(datetime('now'), 'yyMMdd_HHmm')));
    diary(filename_log);

    % summarize parameters
    print_parameter_summary(parameters)

    % Define the filename of the summary table
    parameters.filename_output_table = ...
        fullfile(parameters.output_dir,sprintf('sub-%03d_%s_output_table%s.csv', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    % Add subject_id to parameters
    parameters.subject_id = subject_id;
    
    % suppress unneccessary warnings from export_fig when running without OpenGL
    warning('off','MATLAB:prnRenderer:opengl');
    
    % display GPU information (if requested)
    if strcmp(parameters.code_type, 'cuda') || strcmp(parameters.code_type, 'matlab_gpu')
        fprintf('========================================\n');
        fprintf('GPU INFO \n');
        fprintf('========================================\n\n');
        gpuDevice()
        fprintf('========================================\n\n');
    end

    % set initial time, RAM, GB state
    log_timer('start','single_subject_pipeline', parameters.output_dir);