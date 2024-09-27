function transducer_positioning_with_slurm(subject_id, parameters, pn, target_name, mni_targets, timelimit, memorylimit)
    arguments
        subject_id double
        parameters struct
        pn struct
        target_name string
        mni_targets struct
        timelimit string = "01:00:00" % time limit for a job in seconds (1 hours by default)
        memorylimit (1,1) double = 12 % memory limit for a job in Gb (12 Gb by default)
    end
    if parameters.interactive
        warning('Processing is set to interactive mode, this is not supported when running jobs with qsub, switching off interactive mode.')
        parameters.interactive = 0;
    end
    assert(matches(parameters.overwrite_files,["always","never"]), "When running jobs with qsub, it is not possible to create dialog windows to ask for a confirmation when a file already exists. Set parameters.overwrite_files to 'always' or 'never'");

    log_dir = fullfile(parameters.output_dir, 'batch_job_logs');
    if ~exist(log_dir, 'dir' )
        mkdir(log_dir)
    end
    
    [path_to_pipeline, ~, ~] = fileparts(which('transducer_positioning.m'));
    
    subj_id_string = sprintf('sub-%03d', subject_id);
    % save inputs in the temp file

    temp_data_path = tempname(log_dir);
    [tempdir,tempfile] = fileparts(temp_data_path);
    tempfile = [tempfile '.mat'];
    temp_data_path = [temp_data_path '.mat'];
    save(temp_data_path, "subject_id", "parameters", "pn", "target_name", "mni_targets");

    temp_m_file = tempname(log_dir);
    fid = fopen([temp_m_file '.m'], 'w+');
    fprintf(fid, "load %s; cd %s; transducer_positioning(parameters, pn, subject_id, target_name, mni_targets); delete %s; delete %s;", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
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
    fprintf(fid, '#SBATCH --partition=gpu\n');
    %fprintf(fid, '#SBATCH --gpus=nvidia_a100-sxm4-40gb:2\n');
    %fprintf(fid, '#SBATCH --gpus=nvidia_rtx_a6000:1\n');
    fprintf(fid, '#SBATCH --gres=gpu:1\n');
    %fprintf(fid, '#SBATCH --constraint="cudacap>=8.0"\n');
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