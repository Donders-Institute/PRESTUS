function transducer_positioning_with_qsub(subject_id, parameters, pn, target_name, mni_targets, timelimit, memorylimit)

% TRANSDUCER_POSITIONING_WITH_QSUB Submits transducer positioning jobs to a cluster using Qsub.
%
% This function prepares and submits a transducer positioning job to a cluster using Qsub. 
% It generates temporary MATLAB and Qsub batch script files to execute the 
% `transducer_positioning` function for a given subject and target.
%
% Input:
%   subject_id  - Integer specifying the subject ID.
%   parameters  - Struct containing simulation parameters (e.g., paths, Qsub settings).
%   pn          - Struct containing subject-specific paths (e.g., segmentation folder).
%   target_name - String specifying the name of the target (e.g., 'motor_cortex').
%   mni_targets - Struct containing MNI coordinates for each target.
%   timelimit   - Scalar specifying the time limit for the Qsub job in seconds (default: 3600 seconds or 1 hour).
%   memorylimit - Scalar specifying the memory limit for the Qsub job in GB (default: 12 GB).
%
% Output:
%   None. The function submits the job to the cluster and provides feedback on submission status.

    arguments
        subject_id double
        parameters struct
        pn struct
        target_name string
        mni_targets struct
        timelimit (1,1) double = 60*60*1 % Time limit for a job in seconds (default: 1 hour)
        memorylimit (1,1) double = 12 % Memory limit for a job in GB (default: 12 GB)
    end

    %% Check interactive mode and overwrite settings
    if parameters.interactive
        warning('Interactive mode is not supported when submitting jobs with Qsub. Switching off interactive mode.');
        parameters.interactive = 0;
    end

    assert(matches(parameters.overwrite_files, ["always", "never"]), ...
           "When running jobs with Qsub, dialog windows cannot be used. Set parameters.overwrite_files to 'always' or 'never'.");

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

    matlab_cmd = sprintf('matlab -batch "%s"', temp_m_file_name);

    %% Prepare Qsub submission command
    if ~isfield(parameters, 'qsub_job_prefix')
        parameters.qsub_job_prefix = 'tusim_tp';
    end

    job_name = [parameters.qsub_job_prefix '_' subj_id_string];
    
    qsub_call = sprintf('qsub -N %s -l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=8.0,mem=%iGb,walltime=%i" -o %s -e %s -d %s', ...
        job_name,...
        memorylimit, timelimit, ...
		sprintf('%s_qsub_pipeline_output_$timestamp.log', subj_id_string),...
		sprintf('%s_qsub_pipeline_error_$timestamp.log', subj_id_string),...
		log_dir);

	full_cmd = sprintf('cd %s; timestamp=$(date +%%Y%%m%%d_%%H%%M%%S); echo ''%s'' | %s',  log_dir, matlab_cmd, qsub_call);

	fprintf('Submitted the job to the cluster with a command \n%s \nSee logs in %s in case there are errors. \n', full_cmd, log_dir);
	[res, out] = system(full_cmd);
    
	fprintf('Job name: %s; job ID: %s\n', job_name, out);
end