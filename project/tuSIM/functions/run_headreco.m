function run_headreco(data_path, subject_id, filename_t1, filename_t2, simnibs_env_path, parameters)

    % 3) make logs directory if does not exist
    
    log_dir = fullfile(data_path, 'batch_job_logs');
    if ~exist(log_dir, 'file' )
        mkdir(log_dir)
    end
    
    subj_id_string = sprintf('sub-%03d', subject_id);
    headreco_call = sprintf('headreco all %s %s %s -d no-conform',subj_id_string, filename_t1, filename_t2);
    
    % 4) submit headreco job
	if (parameters.using_donders_hpc)
		qsub_call = sprintf('qsub -l "nodes=1:ppn=1,mem=20Gb,walltime=12:00:00" -o %s -e %s -w %s', ...
			fullfile(log_dir, sprintf('%s_qsub_output.log', subj_id_string)),...
			fullfile(log_dir, sprintf('%s_qsub_error.log', subj_id_string)),...
			data_path);
		% 1) activate anaconda module
		% 2) activate simnibs environment
	%     [res, out] = system(sprintf('module load anaconda3; source activate %s', simnibs_env_path))

		
		full_cmd = sprintf('cd %s; module load anaconda3; source activate %s; echo "cd %s; %s" | %s', data_path, simnibs_env_path, data_path, headreco_call, qsub_call);
		
		fprintf('Running headreco with a command \n%s', full_cmd)
		[res, out] = system(full_cmd);
		display(res);
		display(out);
		
		fprintf('Now wait for the job with the id listed above to finish')
    else
	    fprintf('To get segmented head files, you need to run headreco with a command: \n%s\n The script will stop for now, rerun it when headreco has finished.', headreco_call)
	end    
end