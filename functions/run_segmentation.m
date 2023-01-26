function run_segmentation(data_path, subject_id, filename_t1, filename_t2, parameters)

    % set segmentation path to data_path if no specific seg_path is defined
    if ~isfield(parameters, 'seg_path') || isempty(parameters.seg_path)
        parameters.seg_path = data_path;
    end

    % 1) make logs directory if does not exist
    
    log_dir = fullfile(parameters.seg_path, 'batch_job_logs');
    if ~isfolder(log_dir)
        mkdir(log_dir)
    end
    
    subj_id_string = sprintf('sub-%03d', subject_id);

    % 2) specify segmentation job
    if ~isfield(parameters, 'segmentation_software')
        parameters.segmentation_software = 'charm';
    end

    if strcmp(parameters.segmentation_software, 'charm')
        segment_call = sprintf('charm %s %s %s --forceqform --forcerun',...
            subj_id_string,filename_t1,filename_t2);
        if isfield(parameters, 'use_forceqform') && parameters.use_forceqform == 1
            segment_call = [segment_call '--forceqform'];
        end
    else
        segment_call = sprintf('headreco all %s %s %s -d no-conform',...
            subj_id_string, filename_t1, filename_t2);
    end
    

	if (parameters.using_donders_hpc)
  	qsub_call = sprintf('qsub -N %s -l "nodes=1:ppn=1,mem=20Gb,walltime=24:00:00" -v MANPATH -o %s -e %s -d %s', ...
      ['simnibs-', subj_id_string], ...
          fullfile(log_dir, sprintf('%s_qsub_segment_output_$timestamp.log', subj_id_string)),...
      fullfile(log_dir, sprintf('%s_qsub_segment_error_$timestamp.log', subj_id_string)), ...
          parameters.seg_path);
		

  	% execute simnibs call in segmentation directory
		full_cmd = sprintf('cd %s; timestamp=$(date +%%Y%%m%%d_%%H%%M%%S); echo "%s/%s" | %s', parameters.seg_path, parameters.simnibs_bin_path, segment_call, qsub_call);
		
	  % 3) submit segmentation job
		fprintf('Running segmentation with a command \n%s\n', full_cmd)
		[res, out] = system(full_cmd);
		display(res);
		display(out);
		
		fprintf('Now wait for the job with the id listed above to finish')
    else
	    fprintf('To get segmented head files, you need to run the segmentation software with a command: \n%s\n The script will stop for now, rerun it when the segmentation has finished.', segment_call)
	end    

end
