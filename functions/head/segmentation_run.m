function segmentation_run(data_path, subject_id, filename_t1, filename_t2, parameters)
% SEGMENTATION_RUN  Submit a SimNIBS charm segmentation job
%
% Constructs the SimNIBS charm command for the given subject and dispatches
% it either directly in MATLAB or via an HPC batch scheduler (Slurm or
% PBS/qsub) depending on parameters.platform.
%
% Use as:
%   segmentation_run(data_path, subject_id, filename_t1, filename_t2, parameters)
%
% Input:
%   data_path   - base output directory for segmentation results
%   subject_id  - numeric subject identifier
%   filename_t1 - path to T1-weighted NIfTI file
%   filename_t2 - path to T2-weighted NIfTI file (empty string if not available)
%   parameters  - (1,1) simulation configuration struct
%
% See also: PREPROC_SEGMENTATION, PREPROC_HEAD

    % set segmentation path to data_path if no specific seg_path is defined
    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'seg') || isempty(parameters.path.seg)
        parameters.path.seg = data_path;
    end

    % Create log directory in segmentaion folder (if it does not exist)
    
    log_dir = fullfile(parameters.path.seg, 'batch_job_logs');
    if ~isfolder(log_dir)
        mkdir(log_dir)
    end
    
    subj_id_string = sprintf('sub-%03d', subject_id);

    % Check if the last file produced in the charm pipeline exists. If
    % other files are missing charm will produce the '--forcerun has to
    % be set' error. Setting 'overwrite_simnibs' to 1 will resolve this.
    result_simnibs = sprintf('%sm2m_sub-%03d/final_tissues.nii.gz', parameters.path.seg, subject_id);
    if ~exist(result_simnibs, 'file')
        parameters.io.overwrite_simnibs = 1;
    end
    if ~isempty(filename_t2)
        segment_call = sprintf('charm %s %s %s', subj_id_string, filename_t1, filename_t2);
    else
        segment_call = sprintf('charm %s %s', subj_id_string, filename_t1);
    end
    if isfield(parameters, 'overwrite_simnibs') && parameters.io.overwrite_simnibs == 1
        segment_call = [segment_call ' --forcerun'];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'use_qform') && parameters.segmentation.use_qform == 1
        segment_call = [segment_call ' --forceqform'];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'debug') && parameters.segmentation.debug == 1
        segment_call = [segment_call ' --debug'];
    end
    
    % Platform selection
    if strcmp(parameters.platform, 'auto')
        platform = hpc_detect_system();
        parameters.platform = platform;
        fprintf('➤ auto-detected: %s\n', upper(platform));
    else
        platform = parameters.platform;
        fprintf('➤ deploying: %s\n', upper(platform));
    end

    % Deploy on selected platform
    if strcmp(parameters.platform, 'qsub')
        qsub_call = sprintf('qsub -N %s -l "nodes=1:ppn=1,mem=20Gb,walltime=24:00:00" -v MANPATH -o %s -e %s -d %s', ...
            ['simnibs-', subj_id_string], ...
            fullfile(log_dir, sprintf('%s_qsub_segment_output_$timestamp.log', subj_id_string)),...
            fullfile(log_dir, sprintf('%s_qsub_segment_error_$timestamp.log', subj_id_string)), ...
            parameters.path.seg);

        % execute simnibs call in segmentation directory
        charm_path = fullfile(parameters.startup.simnibs_bin_path, 'charm');
        % Anchor the replacement to the start of segment_call so that 'charm'
        % appearing inside T1/T2 paths (e.g. /home/user/charm_study/...) isn't mangled.
        segment_call_full = regexprep(segment_call, '^charm\b', charm_path);
        full_cmd = sprintf('cd "%s"; timestamp=$(date +%%Y%%m%%d_%%H%%M%%S); echo "%s" | %s', ...
            parameters.path.seg, segment_call_full, qsub_call);
        
        % 3) submit segmentation job
        fprintf('Running segmentation with a command \n%s\n', full_cmd)
        [res, out] = system(full_cmd);
        display(res);
        display(out);
        
        fprintf('Now wait for the job with the id listed above to finish')

    elseif strcmp(parameters.platform, 'slurm')

        % Create a temporary SLURM batch script file
        temp_slurm_file = tempname(log_dir);
        job_name = ['simnibs-', subj_id_string];
        
        fid = fopen([temp_slurm_file '.sh'], 'w+');
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
        fprintf(fid, '#SBATCH --nodes=1\n');
        fprintf(fid, '#SBATCH --ntasks=1\n');
        fprintf(fid, '#SBATCH --mem=20G\n');
        fprintf(fid, '#SBATCH --time=24:00:00\n');
        fprintf(fid, '#SBATCH --output=%s\n', fullfile(log_dir, sprintf('%s_slurm_segment_output_%%j.log', subj_id_string)));
        fprintf(fid, '#SBATCH --error=%s\n', fullfile(log_dir, sprintf('%s_slurm_segment_error_%%j.log', subj_id_string)));
        
        % Add environment setup
        fprintf(fid, 'source /etc/profile\n'); % load system-wide environment variable setups and shell initialization commands
        fprintf(fid, 'export PATH=%s:$PATH\n', parameters.startup.simnibs_bin_path);
        fprintf(fid, 'export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n', parameters.hpc.ld_library_path);
        
        % Add segmentation call
        fprintf(fid, 'cd %s\n', parameters.path.seg);
        fprintf(fid, 'charm --version\n');
        charm_path = fullfile(parameters.startup.simnibs_bin_path, 'charm');
        % Anchor replacement to the start so 'charm' in T1/T2 paths isn't mangled.
        segment_call_full = regexprep(segment_call, '^charm\b', sprintf('"%s"', charm_path));
        fprintf(fid, '%s\n', segment_call_full);
        fclose(fid);
    
        % Ensure script is executable
        system(sprintf('chmod +x %s', [temp_slurm_file,'.sh']));

        % Create the full command to submit the batch script
        sbatch_call = sprintf('sbatch --export=MANPATH %s.sh', temp_slurm_file);
    
        % 4) Submit segmentation job
        fprintf('Running segmentation with a command \n%s\n', sbatch_call);
        [res, out] = system(sbatch_call);
        display(res);
        display(out);
        
        fprintf('Now wait for the job with the id listed above to finish')

    elseif strcmp(parameters.platform, 'matlab')
        fprintf('Running segmentation locally:\n%s\n', segment_call);
        
        if ~isfield(parameters.startup, 'simnibs_bin_path') || isempty(parameters.startup.simnibs_bin_path)
            error('simnibs_bin_path required for local execution');
        end
        
        % Use FULL PATH to charm (don't rely on PATH).
        % Anchor replacement to the start so 'charm' in T1/T2 paths isn't mangled.
        charm_path = fullfile(parameters.startup.simnibs_bin_path, 'charm');
        full_segment_call = regexprep(segment_call, '^charm\b', sprintf('"%s"', charm_path));
        fprintf('Full command: %s\n', full_segment_call);
        
        orig_dir = pwd;
        try
            cd(parameters.path.seg);
            [res, out] = system(full_segment_call);
            cd(orig_dir);
            
            if res == 0
                fprintf('Segmentation completed successfully.\n');
            else
                error('Segmentation failed (exit code %d):\n%s', res, out);
            end
        catch ME
            cd(orig_dir);
            rethrow(ME);
        end

    else
        fprintf('To get segmented head files, you need to run the segmentation software with a command: \n%s\n The script will stop for now, rerun it when the segmentation has finished.\n', segment_call)
        error('Submission medium %s is not available for automatic segmentation.', parameters.platform);
    end 

end
