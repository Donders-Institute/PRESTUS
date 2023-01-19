function run_segmentation(data_path, subject_id, filename_t1, filename_t2, parameters)

    % 3) make logs directory if does not exist
    
    log_dir = fullfile(data_path, 'batch_job_logs');
    if ~exist(log_dir, 'file' )
        mkdir(log_dir)
    end
    
    log_dir_short = 'batch_job_logs/';

    subj_id_string = sprintf('sub-%03d', subject_id);

    if ~isfield(parameters, 'segmentation_software')
        parameters.segmentation_software = 'charm';
    
    end
    if strcmp(parameters.segmentation_software, 'charm')
        segment_call = sprintf('charm %s %s %s',subj_id_string, sprintf(parameters.t1_path_template, subject_id), sprintf(parameters.t1_path_template, subject_id));
    else
        segment_call = sprintf('headreco all %s %s %s -d no-conform',subj_id_string, filename_t1, filename_t2);
    end
    % 4) submit segmentation job
	if (parameters.using_donders_hpc)
		qsub_call = sprintf('qsub -l "nodes=1:ppn=1,mem=20Gb,walltime=12:00:00" -v MANPATH -o %s -e %s -d %s', ...
			fullfile(log_dir_short, sprintf('%s_qsub_segment_output_$timestamp.log', subj_id_string)),...
			fullfile(log_dir_short, sprintf('%s_qsub_segment_error_$timestamp.log', subj_id_string)),...
			data_path);
		% 1) activate anaconda module
		% 2) activate simnibs environment
	%     [res, out] = system(sprintf('module load anaconda3; source activate %s', simnibs_env_path))

		
		full_cmd = sprintf('cd %s; timestamp=$(date +%%Y%%m%%d_%%H%%M%%S); echo "%s/%s" | %s', data_path, parameters.simnibs_bin_path, segment_call, qsub_call);
		
		fprintf('Running segmentation with a command \n%s\n', full_cmd)
		[res, out] = system(full_cmd);
		display(res);
		display(out);
		
		fprintf('Now wait for the job with the id listed above to finish')
    else
	    fprintf('To get segmented head files, you need to run the segmentation software with a command: \n%s\n The script will stop for now, rerun it when the segmentation has finished.', segment_call)
	end    
end