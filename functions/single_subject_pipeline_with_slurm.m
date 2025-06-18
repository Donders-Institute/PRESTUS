function single_subject_pipeline_with_slurm(subject_id, parameters, wait_for_job, timelimit, memorylimit, options)
    arguments
        subject_id double
        parameters struct
        wait_for_job logical = false % boolean for waiting for job to finish before continuing the code
        timelimit string = "04:00:00" % time limit for a job in seconds (4 hours by default)
        memorylimit (1,1) double = 40 % memory limit for a job in Gb (40 Gb by default)
        options.adopted_heatmap (:,:,:) = []
        options.sequential_configs struct = struct()
    end

    % Save that this parameter set is using slurm for further branching
    parameters.submit_medium = 'slurm';

    if parameters.interactive
        warning('Processing is set to interactive mode, this is not supported when running jobs with qsub, switching off interactive mode.')
        parameters.interactive = 0;
    end
    assert(matches(parameters.overwrite_files,["always","never"]), "When running jobs with qsub, it is not possible to create dialog windows to ask for a confirmation when a file already exists. Set parameters.overwrite_files to 'always' or 'never'");

    % Make subfolder (if enabled) and check if directory exists
    % This ensures that log files are saved in the subject subdirectory
    if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
        output_dir = fullfile(parameters.sim_path, sprintf('sub-%03d', subject_id));
    else
        output_dir = fullfile(parameters.sim_path);
    end
    
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end    

    log_dir = fullfile(output_dir, 'batch_job_logs');
    if ~exist(log_dir, 'dir' )
        mkdir(log_dir)
    end

    [path_to_pipeline, ~, ~] = fileparts(which('single_subject_pipeline'));
    
    subj_id_string = sprintf('sub-%03d', subject_id);

    % save inputs in the temp file
    temp_data_path = tempname(log_dir);
    [tempdir,tempfile] = fileparts(temp_data_path);
    tempfile = [tempfile '.mat'];
    temp_data_path = [temp_data_path '.mat'];
    save(temp_data_path, "subject_id", "parameters");

    temp_m_file = tempname(log_dir);
    fid = fopen([temp_m_file '.m'], 'w+');

    % Depending on the input, determine what to save and what to submit to the pipeline
    if ~isempty(fieldnames(options.sequential_configs)) && ~isempty(options.adopted_heatmap)
        sequential_configs = options.sequential_configs;
        adopted_heatmap = options.adopted_heatmap;
        save(temp_data_path, "subject_id", "parameters", "sequential_configs", "adopted_heatmap");
        fprintf(fid, "load '%s'; cd '%s'; single_subject_pipeline(subject_id, parameters, 'sequential_configs', sequential_configs, 'adopted_heatmap', adopted_heatmap); delete '%s'; delete '%s';", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    elseif ~isempty(fieldnames(options.sequential_configs)) && isempty(options.adopted_heatmap)
        sequential_configs = options.sequential_configs;
        save(temp_data_path, "subject_id", "parameters", "sequential_configs");
        fprintf(fid, "load '%s'; cd '%s'; single_subject_pipeline(subject_id, parameters, 'sequential_configs', sequential_configs); delete '%s'; delete '%s';", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    elseif isempty(fieldnames(options.sequential_configs)) && ~isempty(options.adopted_heatmap)
        adopted_heatmap = options.adopted_heatmap;
        save(temp_data_path, "subject_id", "parameters", "adopted_heatmap");
        fprintf(fid, "load '%s'; cd '%s'; single_subject_pipeline(subject_id, parameters, 'adopted_heatmap', adopted_heatmap); delete '%s'; delete '%s';", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    elseif isempty(fieldnames(options.sequential_configs)) && isempty(options.adopted_heatmap)
        save(temp_data_path, "subject_id", "parameters");
        fprintf(fid, "load '%s'; cd '%s'; single_subject_pipeline(subject_id, parameters); delete '%s'; delete '%s';", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    end
    fclose(fid);
    [~,temp_m_file_name,~] = fileparts(temp_m_file);

    if ~isfield(parameters, 'slurm_job_prefix')
        parameters.slurm_job_prefix = 'PRESTUS';
    end

    % Create a temporary SLURM batch script file
    temp_slurm_file = tempname(log_dir);
    job_name = [parameters.slurm_job_prefix '_' subj_id_string];
    fid = fopen([temp_slurm_file '.sh'], 'w+');
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
    if isfield(parameters, 'hpc_partition') && ~isempty(parameters.hpc_partition)
        fprintf(fid, '#SBATCH --partition=%s\n', parameters.hpc_partition);
    elseif strcmp(parameters.code_type, 'matlab_gpu') || strcmp(parameters.code_type, 'cuda')
        fprintf(fid, '#SBATCH --partition=gpu\n');
    end
    if isfield(parameters, 'hpc_gpu') && ~isempty(parameters.hpc_gpu)
        fprintf(fid, '#SBATCH --gres=%s\n', parameters.hpc_gpu);
    elseif strcmp(parameters.code_type, 'matlab_gpu') || strcmp(parameters.code_type, 'cuda')
        fprintf(fid, '#SBATCH --gres=gpu:1\n');
    end
    if isfield(parameters, 'hpc_reservation') && ~isempty(parameters.hpc_reservation)
        fprintf(fid, '#SBATCH --reservation=%s\n', parameters.hpc_reservation);
    end
    fprintf(fid, '#SBATCH --mem=%iG\n', memorylimit);
    fprintf(fid, '#SBATCH --time=%s\n', timelimit);
    fprintf(fid, '#SBATCH --output=%s\n', sprintf('%s_slurm_output_%%j.log', subj_id_string));
    fprintf(fid, '#SBATCH --error=%s\n', sprintf('%s_slurm_error_%%j.log', subj_id_string));
    fprintf(fid, '#SBATCH --chdir=%s\n', log_dir);
    fprintf(fid, 'nvidia-smi\n');
    fprintf(fid, 'module load matlab/R2023b\n');
    fprintf(fid, 'matlab -batch "%s"\n', temp_m_file_name);
    fclose(fid);

    % Create the full command to submit the batch script
    sbatch_call = sprintf('sbatch %s.sh', temp_slurm_file);

    % Execute the full command
    full_cmd = sprintf('cd %s; %s', log_dir, sbatch_call);
	fprintf('Submitted the job to the cluster with a command \n%s \nSee logs in %s in case there are errors. \n', full_cmd, log_dir)
    [status, out] = system(full_cmd);

    job_id = regexp(out, '\d+', 'match');
    job_id = str2double(job_id{1});
    fprintf('Job name: %s; job ID: %i\n', job_name, job_id)
    
    if status == 0
        disp('Job submitted successfully');
        
        % Polling for job status
        job_completed = true;
        if wait_for_job
            disp('User has chosen to wait until job is finished...');
            job_completed = false;
        end
    
        while ~job_completed
            % SLURM equivalent of qstat to check job state
            check_cmd = sprintf('sacct -j %i -o State --noheader | tail -n 1', job_id);
            [status, out] = system(check_cmd);
            
            n_sec = 20; % Pause for n seconds before checking again
            if status == 0
                % Extract the job state (e.g., "RUNNING", "PENDING", "COMPLETED")
                job_state = strtrim(out);
                
                if strcmp(job_state, 'RUNNING') == true
                    disp('Job is still running...');
                    pause(n_sec); % Pause for n seconds before checking again
                elseif strcmp(job_state, 'PENDING') == true
                    disp('Job is still queued...');
                    pause(n_sec); % Pause for n seconds before checking again
                elseif strcmp(job_state, 'COMPLETED') == true
                    disp('Job completed successfully.');
                    job_completed = true;
                else
                    disp(['Job status: ', job_state]);
                    pause(n_sec); % Pause for n seconds before checking again
                    % Additional states: "FAILED", "CANCELLED", etc.
                end
            else
                disp('Failed to check job status.');
                disp(out);
    
                % Handle case where job might have completed and dropped from squeue
                check_cmd = sprintf('scontrol show job %s', job_id);
                [status, out] = system(check_cmd);
                if status ~= 0
                    disp('Job is no longer listed in squeue. Assuming it completed.');
                    job_completed = true;
                else
                    break;
                end
            end
        end
    else
        disp('Command failed to submit the job.');
        disp(out); % Display the error message
    end
    
    % Continue with MATLAB script
    disp('Continuing with the MATLAB script...');
end