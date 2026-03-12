function prestus_pipeline_start(subject_id, parameters, options)
%% PRESTUS_PIPELINE_START  Universal PRESTUS pipeline launcher
%
%   prestus_pipeline_start(subject_id, parameters, options)
%
%   Auto-detects platform and handles direct MATLAB, SLURM, or qsub execution.
%
%   Inputs:
%     subject_id  - Subject number (double)
%     parameters  - Struct with sim_path, platform, hpc_* settings
%     options     - Struct with sequential_configs (default: empty)

    arguments
        subject_id double
        parameters struct
        options struct = struct()
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
            prestus_pipeline(subject_id, parameters, options);
            
        case {'slurm', 'qsub'}
            % ========== HPC EXECUTION ==========
            hpc_validate_parameters(parameters, platform);
            [log_dir, path_to_pipeline, temp_data_path, temp_m_path, temp_m_file] = ...
                hpc_setup_temp_files(parameters, subject_id);
            
            % Create job files
            save(temp_data_path, 'subject_id', 'parameters');
            hpc_matlab_pipeline(temp_m_path, temp_data_path, path_to_pipeline, options);
            
            % Job name
            job_name = hpc_job_name(parameters, subject_id);
            
            % Submit job
            job_id = hpc_submit_job(platform, temp_m_file, parameters, subject_id, log_dir);
            
            % Display job info
            job_info = hpc_job_info(platform, job_id, job_name, subject_id, ...
                parameters.hpc_memorylimit, parameters.hpc_timelimit, log_dir, true);
            
            % Optional wait
            if isfield(parameters, 'hpc_wait_for_job') && parameters.hpc_wait_for_job
                fprintf('⏳ Waiting for job completion...\n');
                fprintf('═══════════════════════════════\n');
                hpc_wait_for_completion(job_id, platform);
                fprintf('✅ Job %s completed\n\n', job_id_display);
            else
                fprintf('➡️  Continuing in MATLAB ...\n\n');
            end
            
            % Save job ID for chaining
            parameters.job_id = job_id;
            
        otherwise
            error('Unknown platform: %s. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', ...
                parameters.platform);
    end
end