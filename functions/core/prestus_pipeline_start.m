function prestus_pipeline_start(parameters, options)
%% PRESTUS_PIPELINE_START  Universal PRESTUS pipeline launcher
%
%   prestus_pipeline_start(parameters)
%   prestus_pipeline_start(parameters, options)
%
%   Auto-detects platform and handles direct MATLAB, SLURM, or qsub execution.
%   Subject ID must be set as parameters.subject_id before calling.
%
%   Inputs:
%     parameters  - Struct; must contain parameters.subject_id
%     options     - Struct with sequential_configs (default: empty)

    arguments
        parameters struct
        options struct = struct()
    end

    if ~isfield(parameters, 'subject_id')
        error('parameters.subject_id must be set before calling prestus_pipeline_start.');
    end
    subject_id = parameters.subject_id;

    % Ensure helper functions are accessible
    helpers_path = fileparts(mfilename('fullpath'));
    if ~contains(path, helpers_path)
        addpath(helpers_path);
    end

    % ========== PLATFORM SELECTION ==========
    if ~isfield(parameters, 'hpc') || ~isfield(parameters, 'platform') || strcmp(parameters.platform, 'auto')
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
            prestus_pipeline(parameters, options);

        case {'slurm', 'qsub'}
            % ========== HPC EXECUTION ==========
            hpc_validate_parameters(parameters, platform);
            [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = ...
                hpc_setup_temp_files(parameters);

            % Create job files
            save(temp_data_path, 'parameters');

            % Generate MATLAB call
            fid = fopen(temp_m_path, 'w+');
            fprintf(fid, 'load(''%s'');\n', temp_data_path);
            fprintf(fid, 'addpath(genpath(''%s''));\n', prestus_path);
            if ismember(fieldnames(options), 'sequential_configs')
                sequential_configs = options.sequential_configs;
                save(temp_data_path, 'sequential_configs', '-append');
                fprintf(fid, 'prestus_pipeline(parameters, options);\n');
            else
                fprintf(fid, 'prestus_pipeline(parameters);\n');
            end
            fprintf(fid, 'delete(''%s'');\n', temp_data_path);
            fprintf(fid, 'delete(''%s'');\n', temp_m_path);
            fclose(fid);

            % Job name
            job_name = hpc_job_name(parameters);

            % Submit job
            [job_id, parameters] = hpc_submit_job(platform, temp_m_file, parameters, log_dir);

            % Display job info
            job_info = hpc_job_info(platform, job_id, job_name, ...
                parameters.hpc.memorylimit, parameters.hpc.timelimit, log_dir, true);

            % Optional wait
            if isfield(parameters.hpc, 'wait_for_job') && parameters.hpc.wait_for_job
                fprintf('⏳ Waiting for job completion...\n');
                fprintf('═══════════════════════════════\n');
                hpc_wait_for_completion(job_id, platform, parameters.hpc.max_wait_checks);
            else
                fprintf('➡️  Continuing in MATLAB ...\n\n');
            end
            
            % Save job ID for chaining
            parameters.hpc.job_id = job_id;

        otherwise
            error('Unknown platform: %s. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', ...
                parameters.platform);
    end
end