function single_subject_pipeline_with_slurm(subject_id, parameters, wait_for_job, timelimit, memorylimit, options)
    arguments
        subject_id double
        parameters struct
        wait_for_job logical = false
        timelimit string = "04:00:00"
        memorylimit (1,1) double = 40
        options.sequential_configs struct = struct()
    end

    % Setup and validation
    parameters.submit_medium = 'slurm';
    validate_slurm_parameters(parameters);
	define_hpc_type(parameters);
    [log_dir, path_to_pipeline] = setup_output_directory(parameters, subject_id);
    
    % Generate consistent temp filenames
    [temp_data_path, temp_m_path, temp_slurm_path, temp_m_file] = generate_temp_files(log_dir);
    
    % Create data file and MATLAB script
    save(temp_data_path, 'subject_id', 'parameters');
    write_matlab_script(temp_m_path, temp_data_path, path_to_pipeline, options.sequential_configs);
    
    % Create and submit SLURM job
    write_slurm_script(temp_slurm_path, parameters, subject_id, timelimit, memorylimit, temp_m_file, log_dir);
    job_id = submit_slurm_job(temp_slurm_path, log_dir, parameters, subject_id);
    
    % Wait for completion if requested
    if wait_for_job
        wait_for_job_completion(job_id);
    end
    
    disp('Continuing with the MATLAB script...');
end

%% Helper functions

function validate_slurm_parameters(parameters)
    if parameters.interactive
        warning('Interactive mode disabled for SLURM jobs.');
        parameters.interactive = false;
    end
    assert(matches(parameters.overwrite_files, ["always", "never"]), ...
        'overwrite_files must be "always" or "never" for SLURM jobs.');
end

function [log_dir, path_to_pipeline] = setup_output_directory(parameters, subject_id)
    if isfield(parameters, 'subject_subfolder') && parameters.subject_subfolder
        output_dir = fullfile(parameters.sim_path, sprintf('sub-%03d', subject_id));
    else
        output_dir = parameters.sim_path;
    end
    
    if ~isfolder(output_dir), mkdir(output_dir); end
    log_dir = fullfile(output_dir, 'batch_job_logs');
    if ~isfolder(log_dir), mkdir(log_dir); end
    
    [path_to_pipeline, ~, ~] = fileparts(which('single_subject_pipeline'));
end

function [temp_data_path, temp_m_path, temp_slurm_path, temp_m_file] = generate_temp_files(log_dir)
    % get timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    temp_base = tempname(log_dir);
    
    % get temporary id (to avoid competing parallel calls)
    [~, temp_base_name] = fileparts(temp_base);
    temp_base_name = temp_base_name(end-7:end);  % Last 8 chars only

    % .mat
    temp_data_path = fullfile(log_dir, sprintf('temp_data_%s_%s.mat', timestamp, temp_base_name));
    % .m
    temp_m_file = sprintf('temp_matlab_%s_%s', timestamp, temp_base_name);
    temp_m_path = fullfile(log_dir, [temp_m_file, '.m']);
    % .sh
    temp_slurm_file = sprintf('temp_slurm_%s_%s.sh', timestamp, temp_base_name);
    temp_slurm_path = fullfile(log_dir, temp_slurm_file);
end

function write_matlab_script(temp_m_path, temp_data_path, path_to_pipeline, sequential_configs)
    fid = fopen(temp_m_path, 'w+');
    
    if ~isempty(fieldnames(sequential_configs))
        save(temp_data_path, 'sequential_configs', '-append');
        fprintf(fid, "load '%s'; cd '%s'; single_subject_pipeline(subject_id, parameters, 'sequential_configs', sequential_configs); delete '%s'; delete '%s';", ...
            temp_data_path, path_to_pipeline, temp_data_path, temp_m_path);
    else
        fprintf(fid, "load '%s'; cd '%s'; single_subject_pipeline(subject_id, parameters); delete '%s'; delete '%s';", ...
            temp_data_path, path_to_pipeline, temp_data_path, temp_m_path);
    end
    fclose(fid);
end

function write_slurm_script(temp_slurm_path, parameters, subject_id, timelimit, memorylimit, temp_m_file, log_dir)
    subj_id_string = sprintf('sub-%03d', subject_id);
    job_name = get_job_name(parameters, subj_id_string);
	
    fid = fopen(temp_slurm_path, 'w+');
    fprintf_slurm_header(fid, job_name, parameters, subj_id_string, log_dir, timelimit, memorylimit);
    if get_gpu_request(parameters)
        fprintf(fid, 'nvidia-smi\n');
    end
	
	if strcmp(parameters.hpc_name, 'snellius')
        fprintf(fid, 'module load 2024\n');
        fprintf(fid, 'module load MATLAB/2024b\n');
    else
		fprintf(fid, 'module load matlab/R2023b\n');
	end
    fprintf(fid, 'matlab -batch "%s"\n', temp_m_file);
    fclose(fid);
end

function job_id = submit_slurm_job(temp_slurm_path, log_dir, parameters, subject_id)
    sbatch_call = sprintf('sbatch %s', temp_slurm_path);
    full_cmd = sprintf('cd %s; %s', log_dir, sbatch_call);
    
    subj_id_string = sprintf('sub-%03d', subject_id);
    job_name = get_job_name(parameters, subj_id_string);
    
    fprintf('Submitted the job to the cluster with a command \n%s \nSee logs in %s\n', full_cmd, log_dir);
    [status, out] = system(full_cmd);
    
    if status ~= 0
        error('SLURM submission failed: %s', out);
    end
    
    job_ids = regexp(out, '\d+', 'match');
    if isempty(job_ids)
        disp(out);
        error('No job ID returned from SLURM');
    end
    job_id = str2double(job_ids{1});
    fprintf('Job "%s" (ID: %i) submitted successfully\n', job_name, job_id);
end

function wait_for_job_completion(job_id)
    disp('User has chosen to wait until job is finished...');
    job_completed = false;
    
    while ~job_completed
        check_cmd = sprintf('sacct -j %i -o State --noheader | tail -n 1', job_id);
        [status_check, out] = system(check_cmd);
        
        n_sec = 20;
        if status_check == 0
            job_state = strtrim(out);
            disp(['Job status: ', job_state]);
            
            if strcmp(job_state, 'RUNNING')
                pause(n_sec);
            elseif strcmp(job_state, 'PENDING')
                pause(n_sec);
            elseif strcmp(job_state, 'COMPLETED')
                disp('Job completed successfully.');
                job_completed = true;
            else
                pause(n_sec);
            end
        else
            disp('Failed to check job status.');
            disp(out);
            check_cmd = sprintf('scontrol show job %s', job_id);
            [status_check, ~] = system(check_cmd);
            if status_check ~= 0
                disp('Job is no longer listed. Assuming it completed.');
                job_completed = true;
            else
                break;
            end
        end
    end
end

function job_name = get_job_name(parameters, subj_id_string)
    if ~isfield(parameters, 'slurm_job_prefix')
        parameters.slurm_job_prefix = 'PRESTUS';
    end
    job_name = [parameters.slurm_job_prefix '_' subj_id_string];
end

function request_gpu = get_gpu_request(parameters)
    if isfield(parameters, 'hpc_partition') && ~isempty(parameters.hpc_partition) && ~strcmp(parameters.hpc_partition, '')
        request_gpu = true;
    elseif strcmp(parameters.code_type, 'matlab_gpu') || strcmp(parameters.code_type, 'cpp_gpu')
        request_gpu = true;
    else
        request_gpu = false;
    end
end

function define_hpc_type(parameters) 

    if ~isfield(parameters, 'hpc_name')
        parameters.hpc_name = 'default';
    end
end

function [memorylimit, cores, n_gpu] = extract_snellius_parameters(parameters)
	
    switch parameters.hpc_partition
        case 'gpu_a100'
            memorylimit = parameters.snellius.a100.memorylimit;
            cores = parameters.snellius.a100.cores;
            max_timelimit = parameters.snellius.a100.timelimit;
            n_gpu = parameters.snellius.a100.n_gpu;
        case 'gpu_h100'
            memorylimit = parameters.snellius.h100.memorylimit;
            cores = parameters.snellius.h100.cores;
            max_timelimit = parameters.snellius.h100.timelimit;
            n_gpu = parameters.snellius.h100.n_gpu;
        otherwise
            error('GPU %s is unknown or not implemented for Snellius.', parameters.hpc_partition)
    end

    assert(timelimit <= max_timelimit, ...
        'Maximum wall time of %s is exceeded (%s).', ...
        max_timelimit, ...
        timelimit)
end

function fprintf_slurm_header(fid, job_name, parameters, subj_id_string, log_dir, timelimit, memorylimit)
	
	% Overwrite and extract additional parameters if Snellius HPC
	is_snellius = strcmp(parameters.hpc_name, 'snellius');
	if is_snellius
		[memorylimit, cores, n_gpu] = extract_snellius_parameters(parameters);
	end
	
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
    
    % Partition
    if isfield(parameters, 'hpc_partition') && ~isempty(parameters.hpc_partition) && ~strcmp(parameters.hpc_partition, '')
		fprintf(fid, '#SBATCH --partition=%s\n', parameters.hpc_partition);
    elseif get_gpu_request(parameters)
        fprintf(fid, '#SBATCH --partition=gpu\n');
    end
    
    % GPU resources
	if is_snellius
		fprintf(fid, '#SBATCH --gpus=%i\n', n_gpu);
    elseif isfield(parameters, 'hpc_gpu') && ~isempty(parameters.hpc_gpu) && ~strcmp(parameters.hpc_gpu, '')
        fprintf(fid, '#SBATCH --gres=%s\n', parameters.hpc_gpu);
    elseif get_gpu_request(parameters)
        fprintf(fid, '#SBATCH --gres=gpu:1\n');
    end
    
    % Reservation
    if isfield(parameters, 'hpc_reservation') && ~isempty(parameters.hpc_reservation) && ~strcmp(parameters.hpc_reservation, '')
        fprintf(fid, '#SBATCH --reservation=%s\n', parameters.hpc_reservation);
    end
    
    % Resources
	if is_snellius
        fprintf(fid, '#SBATCH --cpus-per-task=%i\n', cores);
    end
    fprintf(fid, '#SBATCH --mem=%iG\n', memorylimit);
    fprintf(fid, '#SBATCH --time=%s\n', timelimit);
    fprintf(fid, '#SBATCH --output=%s_slurm_output_%%j.log\n', subj_id_string);
    fprintf(fid, '#SBATCH --error=%s_slurm_error_%%j.log\n', subj_id_string);
    fprintf(fid, '#SBATCH --chdir=%s\n', log_dir);
end