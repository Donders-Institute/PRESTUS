function transducer_positioning_start(subject_id, parameters, pn, target_name, mni_targets)
    arguments
        subject_id double
        parameters struct
        pn struct
        target_name string
        mni_targets struct
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
            transducer_positioning(parameters, pn, subject_id, target_name, mni_targets);
            
        case {'slurm', 'qsub'}
            % ========== HPC EXECUTION ==========
            hpc_validate_parameters(parameters, platform)

            [log_dir, path_to_pipeline, temp_data_path, temp_m_path, temp_m_file] = ...
                hpc_setup_temp_files(parameters, subject_id);
            
            % populate temp data
            save(temp_data_path, "subject_id", "parameters", "pn", "target_name", "mni_targets");
        
            fid = fopen(temp_m_path, 'w');
            fprintf(fid, 'load(''%s'');\n', temp_data_path);
            fprintf(fid, 'cd(''%s'');\n', path_to_pipeline);
            fprintf(fid, 'transducer_positioning(parameters, pn, subject_id, target_name, mni_targets);\n');
            fprintf(fid, 'delete(''%s'');\n', temp_data_path);
            fprintf(fid, 'delete(''%s'');\n', temp_m_path);
            fclose(fid);
        
            parameters.job_prefix = 'TP';
            job_name = hpc_job_name(parameters, subject_id);
            job_id = hpc_submit_job(platform, temp_m_file, parameters, subject_id, log_dir);
            job_info = hpc_job_info(platform, job_id, job_name, subject_id, ...
                parameters.hpc_memorylimit, parameters.hpc_timelimit, log_dir, 1);
        
            if isfield(parameters, 'hpc_wait_for_job') && parameters.hpc_wait_for_job
                fprintf('⏳ Waiting for job completion...\n');
                fprintf('═══════════════════════════════\n');
                hpc_wait_for_completion(job_id, platform);
                fprintf('✅ Job %s completed\n\n', job_id_display);
            end
            
            % Save job ID for chaining
            parameters.job_id = job_id;
            
        otherwise
            error('Unknown platform: %s. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', ...
                parameters.platform);
    end


end
