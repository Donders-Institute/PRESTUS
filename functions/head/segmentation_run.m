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

    % Create log directory in segmentation folder (if it does not exist)
    log_dir = fullfile(parameters.path.seg, 'log_hpc');
    if ~isfolder(log_dir)
        mkdir(log_dir)
    end

    subj_id_string = sprintf('sub-%03d', subject_id);

    % Check if the last file produced in the charm pipeline exists. If
    % other files are missing charm will produce the '--forcerun has to
    % be set' error. Setting 'overwrite_simnibs' to 1 will resolve this.
    result_simnibs = fullfile(parameters.path.seg, ['m2m_' subj_id_string], 'final_tissues.nii.gz');
    if ~exist(result_simnibs, 'file')
        parameters.io.overwrite_simnibs = 1;
    end
    % charm arguments as an argv cell array. The 'matlab' branch passes it
    % to run_command, which runs charm without a shell, so filenames need no
    % quoting. The 'slurm'/'qsub' branches write it into a bash script and
    % quote each argument with shell_quote_arg.
    charm_args = {subj_id_string, filename_t1};
    if ~isempty(filename_t2)
        charm_args = [charm_args, {filename_t2}];
    end
    if isfield(parameters.io, 'overwrite_simnibs') && parameters.io.overwrite_simnibs == 1
        charm_args = [charm_args, {'--forcerun'}];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'use_qform') && parameters.segmentation.use_qform == 1
        charm_args = [charm_args, {'--forceqform'}];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'debug') && parameters.segmentation.debug == 1
        charm_args = [charm_args, {'--debug'}];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'noneck') && parameters.segmentation.noneck == 1
        charm_args = [charm_args, {'--noneck'}];
    end
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'settings_file') ...
            && ~isempty(parameters.segmentation.settings_file)
        settings_file = char(parameters.segmentation.settings_file);
        if ~isfile(settings_file)
            error('segmentation_run:settingsMissing', 'segmentation.settings_file not found: %s', settings_file);
        end
        % charm --forcerun deletes the m2m folder, and refuses a settings file inside it
        % (absolute: charm runs in path.seg, where a relative path would point elsewhere)
        d = dir(settings_file);
        settings_file = fullfile(d(1).folder, d(1).name);
        m2m_folder = fullfile(parameters.path.seg, ['m2m_' subj_id_string]);
        dm = dir(m2m_folder);                   % entries carry the folder's absolute path
        if ~isempty(dm) && isfolder(m2m_folder) && strcmp(d(1).folder, dm(1).folder)
            error('segmentation_run:settingsInM2m', ...
                'Keep segmentation.settings_file outside %s: charm deletes that folder when it re-runs.', m2m_folder);
        end
        charm_args = [charm_args, {'--usesettings', settings_file}];
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
    % Segmentation labels (charm_seg_labels.m) are those of charm, SimNIBS 4 and later
    info = simnibs_install_info(char(parameters.startup.simnibs_bin_path));
    if ~isnan(info.major) && info.major < 4
        error('segmentation_run:simnibsVersion', ...
            'SimNIBS %s in %s: PRESTUS segments with charm from SimNIBS 4 or later.', ...
            info.version, parameters.startup.simnibs_bin_path);
    elseif isempty(info.version)
        warn('SimNIBS version of %s not found; PRESTUS needs SimNIBS 4 or later.', ...
            parameters.startup.simnibs_bin_path);
    else
        fprintf('➤ SimNIBS %s (%s)\n', info.version, parameters.startup.simnibs_bin_path);
    end
    charm_path = fullfile(parameters.startup.simnibs_bin_path, 'charm');
    % Shell-safe text form of the charm call, for the two branches below
    % that write it into a bash script executed later (by a compute-node
    % shell, for slurm/qsub) rather than run directly.
    charm_line = strjoin([{shell_quote_arg(charm_path)}, ...
        cellfun(@shell_quote_arg, charm_args, 'UniformOutput', false)], ' ');

    % Deploy on selected platform
    if strcmp(parameters.platform, 'qsub')
        % Job output/error log filenames carry a timestamp set in MATLAB.
        timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        qsub_argv = {'qsub', ...
            '-N', ['simnibs-', subj_id_string], ...
            '-l', 'nodes=1:ppn=1,mem=20Gb,walltime=24:00:00', ...
            '-v', 'MANPATH', ...
            '-o', fullfile(log_dir, sprintf('%s_qsub_segment_output_%s.log', subj_id_string, timestamp)), ...
            '-e', fullfile(log_dir, sprintf('%s_qsub_segment_error_%s.log', subj_id_string, timestamp)), ...
            '-d', parameters.path.seg};

        % qsub reads the job script from stdin when given no file operand;
        % run_command passes charm_line as Stdin, with no local shell.
        % charm_line is quoted because the compute-node shell runs it.
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
        % hpc.ld_library_path, else the selected install's libraries (SIMNIBS_LD_LIBRARY_PATH)
        ld = simnibs_ld_library_path(parameters);
        if ~isempty(ld)
            fprintf(fid, 'export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n', shell_quote_arg(ld));
        end

        % Add segmentation call
        fprintf(fid, 'cd %s\n', shell_quote_arg(parameters.path.seg));
        fprintf(fid, 'charm --version\n');
        fprintf(fid, '%s\n', charm_line);
        fclose(fid);

        % Make the script executable (no shell needed)
        java.io.File([temp_slurm_file '.sh']).setExecutable(true);

        % Submit the batch script via argv (no local shell)
        fprintf('Running segmentation with a command \nsbatch --export=MANPATH %s.sh\n', temp_slurm_file);
        [res, out] = run_command({'sbatch', '--export=MANPATH', [temp_slurm_file '.sh']});
        display(res);
        display(out);

        fprintf('Now wait for the job with the id listed above to finish')

    elseif strcmp(parameters.platform, 'matlab')
        % run_command runs charm via argv in parameters.path.seg, without a shell
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
