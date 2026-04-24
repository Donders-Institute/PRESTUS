function job_id = prestus_pipeline_start(parameters, options)
% PRESTUS_PIPELINE_START  Universal PRESTUS pipeline launcher for MATLAB, SLURM, and qsub
%
% Auto-detects the execution platform (matlab / slurm / qsub) via
% hpc_detect_system and dispatches accordingly. On MATLAB, calls
% prestus_pipeline directly. On HPC, serialises parameters to a temp .mat
% file, generates a bootstrap .m script, and submits via hpc_submit_job.
% When parameters.simulation.uncertainty is true, delegates to
% uncertainty_pipeline instead of direct submission.
%
% Use as:
%   prestus_pipeline_start(parameters)
%   prestus_pipeline_start(parameters, options)
%   job_id = prestus_pipeline_start(parameters)
%   job_id = prestus_pipeline_start(parameters, options)
%
% Input:
%   parameters - PRESTUS config; must contain subject_id; platform field
%                controls dispatch ('matlab','slurm','qsub','auto')
%   options    - pipeline options (optional, default: struct());
%                .sequential_configs — cell array of config paths for sequential runs
%
% Output:
%   job_id     - HPC scheduler job ID; [] when platform is 'matlab'
%
% See also: PRESTUS_PIPELINE, UNCERTAINTY_PIPELINE, LOAD_PARAMETERS, PATH_LOG_SETUP

    arguments
        parameters struct
        options struct = struct()
    end

    if ~isfield(parameters, 'subject_id')
        error('parameters.subject_id must be set before calling prestus_pipeline_start.');
    end

    % Ensure helper functions are accessible
    helpers_path = fileparts(mfilename('fullpath'));
    if ~contains(path, helpers_path)
        addpath(helpers_path);
    end

    % ========== PLATFORM SELECTION ==========
    if ~isfield(parameters, 'platform') || strcmp(parameters.platform, 'auto')
        platform = hpc_detect_system();
        parameters.platform = platform;
        fprintf('➤ auto-detected: %s\n', upper(platform));
    else
        platform = parameters.platform;
        fprintf('➤ deploying: %s\n', upper(platform));
    end

    % ========== DISPATCH EXECUTION ==========
    switch parameters.platform
        case 'matlab'
            fprintf('🖥️  Running in MATLAB\n\n');
            if is_multi_isppa_mode(parameters)
                multi_isppa_pipeline(parameters, options);
            else
                prestus_pipeline(parameters, options);
            end
            job_id = [];

        case {'auto', 'slurm', 'qsub'}
            % ========== HPC EXECUTION ==========
            % Intercept uncertainty mode before single-job submission:
            % uncertainty_pipeline manages its own multi-stage job graph.
            if isfield(parameters.simulation, 'uncertainty') && parameters.simulation.uncertainty
                uncertainty_pipeline(parameters, options);
                job_id = [];
                return;
            end

            % Intercept multi-ISPPA mode before single-job submission:
            % multi_isppa_pipeline manages its own multi-stage job graph.
            if is_multi_isppa_mode(parameters)
                multi_isppa_pipeline(parameters, options);
                job_id = [];
                return;
            end
            hpc_validate_parameters(parameters, platform);
            [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = ...
                hpc_setup_temp_files(parameters);

            % Create job files
            save(temp_data_path, 'parameters');

            % Generate MATLAB call
            fid = fopen(temp_m_path, 'w+');
            fprintf(fid, 'load(''%s'');\n', temp_data_path);
            % Bootstrap safe_addpath, then use it to add all PRESTUS paths
            % without including hidden directories (e.g. .claude worktrees).
            fprintf(fid, 'addpath(''%s'');\n', fullfile(prestus_path, 'functions', 'helper'));
            fprintf(fid, 'safe_addpath(''%s'');\n', fullfile(prestus_path, 'functions'));
            fprintf(fid, 'safe_addpath(''%s'');\n', fullfile(prestus_path, 'toolboxes'));
            if isfield(options, 'sequential_configs')
                save(temp_data_path, 'options', '-append');
                fprintf(fid, 'prestus_pipeline(parameters, options);\n');
            else
                fprintf(fid, 'prestus_pipeline(parameters);\n');
            end
            fprintf(fid, 'delete(''%s'');\n', temp_data_path);
            fprintf(fid, 'delete(''%s'');\n', temp_m_path);
            fclose(fid);

            % Submit job
            [job_id, parameters] = hpc_submit_job(platform, temp_m_file, parameters, log_dir);

            % Optional wait
            if isfield(parameters.hpc, 'wait_for_job') && parameters.hpc.wait_for_job
                fprintf('⏳ Waiting for job completion...\n');
                fprintf('═══════════════════════════════\n');
                hpc_wait_for_completion(job_id, platform, parameters.hpc.max_wait_checks);
            else
                fprintf('➡️  Continuing in MATLAB ...\n\n');
            end
            
            % job_id is already set by hpc_submit_job above

        otherwise
            error('Unknown platform: %s. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', ...
                parameters.platform);
    end
end

function tf = is_multi_isppa_mode(parameters)
    tf = isfield(parameters, 'calibration') && ...
         isfield(parameters.calibration, 'target_isppa_wcm2') && ...
         numel(parameters.calibration.target_isppa_wcm2) > 1;
end