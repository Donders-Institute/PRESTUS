function prestus_pipeline_start(subject_id, parameters, options)
%% PRESTUS_PIPELINE_START  Universal PRESTUS pipeline launcher
%
%   prestus_pipeline_start(subject_id, parameters, options)
%
%   Auto-detects platform and handles direct MATLAB, SLURM, or qsub execution.
%
%   Inputs:
%     subject_id  - Subject number (double)
%     parameters  - Struct with sim_path, submit_medium, hpc_* settings
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

    % ========== STEP 1: PLATFORM AUTO-DETECTION ==========
    if ~isfield(parameters, 'submit_medium') || strcmp(parameters.submit_medium, 'auto')
        submit_medium = hpc_detect_system();
        parameters.submit_medium = submit_medium;
        fprintf('➤ auto-detected: %s\n', upper(submit_medium));
    else
        submit_medium = parameters.submit_medium;
        fprintf('➤ deploying: %s\n', upper(submit_medium));
    end
    
    % ========== DISPATCH EXECUTION ==========
    switch parameters.submit_medium
        case 'matlab'
            fprintf('🖥️  Running in MATLAB\n\n');
            prestus_pipeline(subject_id, parameters, options);
            
        case {'slurm', 'qsub'}
            % ========== HPC EXECUTION ==========
            hpc_validate_parameters(parameters, submit_medium);
            [log_dir, path_to_pipeline, temp_data_path, temp_m_path, temp_m_file] = ...
                hpc_setup_temp_files(parameters, subject_id);
            
            % Create job files
            save(temp_data_path, 'subject_id', 'parameters');
            hpc_matlab_pipeline(temp_m_path, temp_data_path, path_to_pipeline, options);
            
            % Job name
            job_name = hpc_job_name(submit_medium, parameters, subject_id);
            
            % Submit job
            job_id = hpc_submit_job(submit_medium, temp_m_file, parameters, subject_id, log_dir);
            
            % Display job info
            job_info = hpc_job_info(submit_medium, job_id, job_name, subject_id, ...
                parameters.hpc_memorylimit, parameters.hpc_timelimit, log_dir, true);
            
            % Optional wait
            if isfield(parameters, 'hpc_wait_for_job') && parameters.hpc_wait_for_job
                fprintf('⏳ Waiting for job completion...\n');
                fprintf('═══════════════════════════════\n');
                hpc_wait_for_completion(job_id, submit_medium);
                fprintf('✅ Job %s completed\n\n', job_id_display);
            else
                fprintf('➡️  Continuing in MATLAB ...\n\n');
            end
            
            % Save job ID for chaining
            parameters.job_id = job_id;
            
        otherwise
            error('Unknown submit_medium: %s. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', ...
                parameters.submit_medium);
    end
end