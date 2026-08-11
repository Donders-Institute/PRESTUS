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

    log_dir = fullfile(parameters.path.seg, 'log_hpc');
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
    % Build the charm arguments as an argv cell array rather than a
    % shell command string. For the 'matlab' branch below this is handed
    % straight to run_command(), which invokes charm directly with no
    % shell at all -- filename_t1/filename_t2 (which come from scanned
    % raw-data filenames) reach the process as literal argument values
    % regardless of their contents, with nothing to quote or escape.
    % The 'slurm'/'qsub' branches still need this as shell-safe text
    % (their bash scripts are executed on a compute node later), so
    % shell_quote_arg() is applied there when the line is composed.
    charm_args = {subj_id_string, filename_t1};
    if ~isempty(filename_t2)
        charm_args = [charm_args, {filename_t2}];
    end
    if isfield(parameters, 'overwrite_simnibs') && parameters.io.overwrite_simnibs == 1
        charm_args = [charm_args, {'--forcerun'}];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'use_qform') && parameters.segmentation.use_qform == 1
        charm_args = [charm_args, {'--forceqform'}];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'debug') && parameters.segmentation.debug == 1
        charm_args = [charm_args, {'--debug'}];
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

    if ~isfield(parameters, 'startup') || ~isfield(parameters.startup, 'simnibs_bin_path') ...
            || isempty(parameters.startup.simnibs_bin_path)
        error('simnibs_bin_path required for segmentation');
    end
    charm_path = fullfile(parameters.startup.simnibs_bin_path, 'charm');
    % Shell-safe text form of the charm call, for the two branches below
    % that write it into a bash script executed later (by a compute-node
    % shell, for slurm/qsub) rather than run directly.
    charm_line = strjoin([{shell_quote_arg(charm_path)}, ...
        cellfun(@shell_quote_arg, charm_args, 'UniformOutput', false)], ' ');

    % Deploy on selected platform
    if strcmp(parameters.platform, 'qsub')
        % Job output/error log filenames used to embed a shell-side
        % `$timestamp` set via `timestamp=$(date ...)` in the same command
        % line as `echo ... | qsub ...`. Computing it directly in MATLAB
        % removes the need for that shell altogether.
        timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        qsub_argv = {'qsub', ...
            '-N', ['simnibs-', subj_id_string], ...
            '-l', 'nodes=1:ppn=1,mem=20Gb,walltime=24:00:00', ...
            '-v', 'MANPATH', ...
            '-o', fullfile(log_dir, sprintf('%s_qsub_segment_output_%s.log', subj_id_string, timestamp)), ...
            '-e', fullfile(log_dir, sprintf('%s_qsub_segment_error_%s.log', subj_id_string, timestamp)), ...
            '-d', parameters.path.seg};

        % qsub reads the job script from stdin when given no file operand
        % (the original `echo <script> | qsub ...` pattern). run_command's
        % Stdin option replicates this directly via ProcessBuilder, with
        % no shell pipe and no outer shell parsing the qsub invocation
        % itself; charm_line still needs its own shell_quote_arg()
        % quoting because it IS executed by a shell -- just the one on
        % the compute node when the job actually runs, not a local one.
        fprintf('Running segmentation with a command \n%s\n', charm_line)
        [res, out] = run_command(qsub_argv, 'Stdin', charm_line);
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
        fprintf(fid, 'cd %s\n', shell_quote_arg(parameters.path.seg));
        fprintf(fid, 'charm --version\n');
        fprintf(fid, '%s\n', charm_line);
        fclose(fid);

        % Ensure script is executable. java.io.File.setExecutable (not
        % system('chmod +x ...')) needs no shell for a plain, internally
        % generated path.
        java.io.File([temp_slurm_file '.sh']).setExecutable(true);

        % Submit the batch script via argv -- no shell involved locally.
        fprintf('Running segmentation with a command \nsbatch --export=MANPATH %s.sh\n', temp_slurm_file);
        [res, out] = run_command({'sbatch', '--export=MANPATH', [temp_slurm_file '.sh']});
        display(res);
        display(out);

        fprintf('Now wait for the job with the id listed above to finish')

    elseif strcmp(parameters.platform, 'matlab')
        % run_command invokes charm directly via argv, in parameters.path.seg,
        % with no shell at all -- filename_t1/filename_t2 reach the process
        % as literal argument values regardless of their contents.
        fprintf('Running segmentation locally:\n%s\n', charm_line);
        [res, out] = run_command([{charm_path}, charm_args], 'WorkingDir', parameters.path.seg);

        if res == 0
            fprintf('Segmentation completed successfully.\n');
        else
            error('Segmentation failed (exit code %d):\n%s', res, out);
        end

    else
        fprintf('To get segmented head files, you need to run the segmentation software with a command: \n%s\n The script will stop for now, rerun it when the segmentation has finished.\n', charm_line)
        error('Submission medium %s is not available for automatic segmentation.', parameters.platform);
    end

end
