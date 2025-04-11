function transducer_positioning_with_slurm(subject_id, parameters, pn, target_name, mni_targets, timelimit, memorylimit)

% TRANSDUCER_POSITIONING_WITH_SLURM Submits transducer positioning jobs to a SLURM cluster.
%
% This function prepares and submits a transducer positioning job to a SLURM-based 
% high-performance computing (HPC) cluster. It generates temporary MATLAB and SLURM 
% batch script files to execute the `transducer_positioning` function for a given 
% subject and target.
%
% Input:
%   subject_id  - Integer specifying the subject ID.
%   parameters  - Struct containing simulation parameters (e.g., paths, SLURM settings).
%   pn          - Struct containing subject-specific paths (e.g., segmentation folder).
%   target_name - String specifying the name of the target (e.g., 'motor_cortex').
%   mni_targets - Struct containing MNI coordinates for each target.
%   timelimit   - String specifying the time limit for the SLURM job (default: '01:00:00').
%   memorylimit - Scalar specifying the memory limit for the SLURM job in GB (default: 12).
%
% Output:
%   None. The function submits the job to the cluster and provides feedback on submission status.

    arguments
        subject_id double
        parameters struct
        pn struct
        target_name string
        mni_targets struct
        timelimit string = "01:00:00" % Time limit for a job in seconds (default: 1 hour)
        memorylimit (1,1) double = 12 % Memory limit for a job in GB (default: 12 GB)
    end

    %% Check interactive mode and overwrite settings
    if parameters.interactive
        warning('Interactive mode is not supported when submitting jobs with SLURM. Switching off interactive mode.');
        parameters.interactive = 0;
    end

    assert(matches(parameters.overwrite_files, ["always", "never"]), ...
           "When running jobs with SLURM, dialog windows cannot be used. Set parameters.overwrite_files to 'always' or 'never'.");

    %% Create log directory if it does not exist
    log_dir = fullfile(parameters.output_dir, 'batch_job_logs');
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end

    %% Determine pipeline location and prepare temporary files
    [path_to_pipeline, ~, ~] = fileparts(which('transducer_positioning.m'));

    subj_id_string = sprintf('sub-%03d', subject_id);

    % Save inputs in a temporary MAT file
    temp_data_path = tempname(log_dir);
    temp_data_path = [temp_data_path '.mat'];
    save(temp_data_path, "subject_id", "parameters", "pn", "target_name", "mni_targets");

    % Create temporary MATLAB script file
    temp_m_file = tempname(log_dir);
    fid = fopen([temp_m_file '.m'], 'w+');
    fprintf(fid, "load %s; cd %s; transducer_positioning(parameters, pn, subject_id, target_name, mni_targets); delete %s; delete %s;", ...
            temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    fclose(fid);

    [~, temp_m_file_name, ~] = fileparts(temp_m_file);

    %% Prepare SLURM batch script file
    if ~isfield(parameters, 'slurm_job_prefix')
        parameters.slurm_job_prefix = 'PRESTUS';
    end

    temp_slurm_file = tempname(log_dir);
    job_name = [parameters.slurm_job_prefix '_' subj_id_string];
    
    fid = fopen([temp_slurm_file '.sh'], 'w+');
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
    if isfield(parameters, 'hcp_partition') && ~isempty(parameters.hcp_partition)
        fprintf(fid, '#SBATCH --partition=%s\n', parameters.hcp_partition);
    else
        fprintf(fid, '#SBATCH --partition=gpu\n');
    end
    if isfield(parameters, 'hcp_gpu') && ~isempty(parameters.hcp_gpu)
        fprintf(fid, '#SBATCH --gres=%s\n', parameters.hcp_gpu);
    else
        fprintf(fid, '#SBATCH --gres=gpu:1\n');
    end
    fprintf(fid, '#SBATCH --mem=%iG\n', memorylimit);
    fprintf(fid, '#SBATCH --time=%s\n', timelimit);
    fprintf(fid, '#SBATCH --output=%s\n', sprintf('%s_slurm_output_%%j.log', subj_id_string));
    fprintf(fid, '#SBATCH --error=%s\n', sprintf('%s_slurm_error_%%j.log', subj_id_string));
    fprintf(fid, '#SBATCH --chdir=%s\n', log_dir);
    fprintf(fid, 'module load matlab\n');
    fprintf(fid, 'matlab -batch "%s"\n', temp_m_file_name);
    fclose(fid);

    % Create the full command to submit the batch script
    sbatch_call = sprintf('sbatch %s.sh', temp_slurm_file);

    % Execute the full command
    full_cmd = sprintf('cd %s; %s', log_dir, sbatch_call);
	fprintf('Submitted the job to the cluster with a command \n%s \nSee logs in %s in case there are errors. \n', full_cmd, log_dir)
	[res, out] = system(full_cmd);
    
    fprintf('Job name: %s; job ID: %s', job_name, out)
end