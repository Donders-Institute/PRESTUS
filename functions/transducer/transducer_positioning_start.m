function transducer_positioning_start(parameters, pn, target_name, mni_targets)
    arguments
        parameters struct
        pn struct
        target_name string
        mni_targets struct
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
            transducer_positioning(parameters, pn, target_name, mni_targets);
            
        case {'slurm', 'qsub'}
            % ========== HPC EXECUTION ==========
            hpc_validate_parameters(parameters, platform)

            [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = ...
                hpc_setup_temp_files(parameters);

            % populate temp data
            save(temp_data_path, "parameters", "pn", "target_name", "mni_targets");

            fid = fopen(temp_m_path, 'w');
            fprintf(fid, 'load(''%s'');\n', temp_data_path);
            fprintf(fid, 'addpath(genpath(''%s''));\n', prestus_path);
            fprintf(fid, 'transducer_positioning(parameters, pn, target_name, mni_targets);\n');
            fprintf(fid, 'delete(''%s'');\n', temp_data_path);
            fprintf(fid, 'delete(''%s'');\n', temp_m_path);
            fclose(fid);

            parameters.hpc.job_prefix = 'TP';
            job_name = hpc_job_name(parameters);
            job_id = hpc_submit_job(platform, temp_m_file, parameters, log_dir);
            job_info = hpc_job_info(platform, job_id, job_name, ...
                parameters.hpc.memorylimit, parameters.hpc.timelimit, log_dir, 1);

            if isfield(parameters.hpc, 'wait_for_job') && parameters.hpc.wait_for_job
                fprintf('⏳ Waiting for job completion...\n');
                fprintf('═══════════════════════════════\n');
                hpc_wait_for_completion(job_id, platform, parameters.hpc.max_wait_checks);
                fprintf('✅ Job %s completed\n\n', job_id_display);
            end
            
            % Save job ID for chaining
            parameters.hpc.job_id = job_id;
            
        otherwise
            error('Unknown platform: %s. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', ...
                parameters.platform);
    end


end
