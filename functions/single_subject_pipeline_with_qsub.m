function single_subject_pipeline_with_qsub(subject_id, parameters, timelimit, memorylimit)
    arguments
        subject_id double
        parameters struct
        timelimit (1,1) double = 60*60*4 % time limit for a job in seconds (4 hours by default)
        memorylimit (1,1) double = 40 % memory limit for a job in Gb (40 Gb by default)
    end
    if parameters.interactive
        warning('Processing is set to interactive mode, this is not supported when running jobs with qsub, switching off interactive mode.')
        parameters.interactive = 0;
    end
    assert(matches(parameters.overwrite_files,["always","never"]), "When running jobs with qsub, it is not possible to create dialog windows to ask for a confirmation when a file already exists. Set parameters.overwrite_files to 'always' or 'never'");

    % Make subfolder (if enabled) and check if directory exists
    % This ensures that log files are saved in the subject subdirectory
    if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
        output_dir = fullfile(parameters.temp_output_dir, sprintf('sub-%03d', subject_id));
    else
        output_dir = fullfile(parameters.temp_output_dir);
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
    fprintf(fid, "load %s; cd %s; single_subject_pipeline(subject_id, parameters); delete %s; delete %s;", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    fclose(fid);
    [~,temp_m_file_name,~] = fileparts(temp_m_file);
    
    matlab_cmd = sprintf('matlab -batch "%s"', temp_m_file_name);

    if ~isfield(parameters, 'qsub_job_prefix')
        parameters.qsub_job_prefix = 'PRESTUS';
    end
    job_name = [parameters.qsub_job_prefix '_' subj_id_string];
    qsub_call = sprintf('qsub -N %s -l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=8.0,mem=%iGb,walltime=%i" -o %s -e %s -d %s', ...
        job_name,...
        memorylimit, timelimit, ...
		sprintf('%s_qsub_pipeline_output_$timestamp.log', subj_id_string),...
		sprintf('%s_qsub_pipeline_error_$timestamp.log', subj_id_string),...
		log_dir);


	full_cmd = sprintf('cd %s; timestamp=$(date +%%Y%%m%%d_%%H%%M%%S); echo ''%s'' | %s',  log_dir, matlab_cmd, qsub_call);
	
	%fprintf('Submitted the job to the cluster with a command \n%s \nSee logs in %s in case there are errors. \n', full_cmd, log_dir)
    fprintf('Submitted the job to the cluster.\nSee logs in %s in case there are errors. \n', log_dir)
	[res, out] = system(full_cmd);
    
    fprintf('Job name: %s; job ID: %s', job_name, out)
end