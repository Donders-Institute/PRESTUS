function job_id = calibration_pipeline_start(...
    profile_empirical, ...
    equipment_name, ...
    desired_intensity, ...
    desired_focal_distance_ep, ...
    parameters, ...
    sim_id)
% CALIBRATION_PIPELINE_START  Universal launcher for calibration_transducer (MATLAB / SLURM / qsub)
%
% Mirrors prestus_pipeline_start for the calibration module. On MATLAB,
% calls calibration_transducer directly. On HPC, serialises all calibration
% arguments to a temp .mat file, generates a bootstrap .m script, and submits
% via hpc_submit_job.
%
% When running on HPC, the inner water simulations executed inside
% calibration_transducer run locally on the allocated node (parameters.platform
% is forced to 'matlab' in the serialised payload), so the job must request the
% GPU resources required by parameters.simulation.code_type.
%
% Use as:
%   calibration_pipeline_start(profile_empirical, equipment_name, ...
%       desired_intensity, desired_focal_distance_ep, parameters, sim_id)
%   job_id = calibration_pipeline_start(...)
%
% Input:
%   profile_empirical         - struct with axial_intensity and axial_distance_bowl [mm]
%   equipment_name            - name/identifier for the transducer or equipment
%   desired_intensity         - target focal intensity/intensities [W/cm²] (scalar or row vector)
%   desired_focal_distance_ep - desired focal distance from transducer exit plane [mm]
%   parameters                - PRESTUS config with calibration settings;
%                               parameters.platform controls dispatch
%                               ('matlab' | 'slurm' | 'qsub' | 'auto')
%   sim_id                    - numeric simulation ID (used for job naming and output folders)
%
% Output:
%   job_id - HPC scheduler job ID; [] when platform is 'matlab'
%
% See also: CALIBRATION_TRANSDUCER, PRESTUS_PIPELINE_START, HPC_SUBMIT_JOB

    %% Platform selection
    if ~isfield(parameters, 'platform') || strcmp(parameters.platform, 'auto')
        platform = hpc_detect_system();
        parameters.platform = platform;
        fprintf('calibration_pipeline_start: auto-detected platform: %s\n', upper(platform));
    else
        platform = parameters.platform;
        fprintf('calibration_pipeline_start: deploying on: %s\n', upper(platform));
    end

    %% Dispatch

    switch platform

        % ------------------------------------------------------------------
        case 'matlab'
            calibration_transducer(...
                profile_empirical, ...
                equipment_name, ...
                desired_intensity, ...
                desired_focal_distance_ep, ...
                parameters, ...
                sim_id);
            job_id = [];

        % ------------------------------------------------------------------
        case {'slurm', 'qsub'}

            % subject_id is needed by hpc_submit_job / hpc_job_name
            parameters.subject_id = sim_id;

            % Non-interactive is required for unattended HPC execution
            parameters.simulation.interactive = 0;

            % Use calibration-specific job prefix so jobs are distinguishable
            if ~isfield(parameters, 'hpc') || ~isfield(parameters.hpc, 'job_prefix') ...
                    || isempty(parameters.hpc.job_prefix)
                parameters.hpc.job_prefix = 'CAL';
            end

            hpc_validate_parameters(parameters, platform);

            %% Set up log directory under calibration output
            cal_output = parameters.calibration.path_output;
            log_dir    = fullfile(cal_output, sprintf('sub-%03d', sim_id), 'log_hpc');
            if ~isfolder(log_dir), mkdir(log_dir); end

            prestus_path = get_prestus_path();

            %% Generate temp file names
            timestamp      = datestr(now, 'yyyymmdd_HHMMSS');
            temp_base      = tempname(log_dir);
            [~, base_name] = fileparts(temp_base);
            base_name      = base_name(end-7:end);

            temp_data_path = fullfile(log_dir, ...
                sprintf('temp_data_%s_%s.mat', timestamp, base_name));
            temp_m_file    = sprintf('temp_matlab_%s_%s', timestamp, base_name);
            temp_m_path    = fullfile(log_dir, [temp_m_file, '.m']);

            %% Serialise calibration arguments
            % Force inner water simulations to run locally on the HPC node —
            % we do not want nested HPC submission from within the job.
            parameters.platform = 'matlab';

            cal_args.profile_empirical         = profile_empirical;
            cal_args.equipment_name            = equipment_name;
            cal_args.desired_intensity         = desired_intensity;
            cal_args.desired_focal_distance_ep = desired_focal_distance_ep;
            cal_args.sim_id                    = sim_id;

            save(temp_data_path, 'parameters', 'cal_args');

            %% Generate bootstrap MATLAB script
            fid = fopen(temp_m_path, 'w+');
            fprintf(fid, 'load(''%s'');\n', temp_data_path);
            % Bootstrap safe_addpath, then add all PRESTUS paths (excluding
            % hidden directories such as .claude worktrees).
            fprintf(fid, 'addpath(''%s'');\n', fullfile(prestus_path, 'functions', 'helper'));
            fprintf(fid, 'safe_addpath(''%s'');\n', fullfile(prestus_path, 'functions'));
            fprintf(fid, 'safe_addpath(''%s'');\n', fullfile(prestus_path, 'external'));
            fprintf(fid, 'calibration_transducer(cal_args.profile_empirical, cal_args.equipment_name, cal_args.desired_intensity, cal_args.desired_focal_distance_ep, parameters, cal_args.sim_id);\n');
            fprintf(fid, 'delete(''%s'');\n', temp_data_path);
            fprintf(fid, 'delete(''%s'');\n', temp_m_path);
            fclose(fid);

            %% Submit job
            [job_id, parameters] = hpc_submit_job(platform, temp_m_file, parameters, log_dir);

            %% Optional blocking wait
            if isfield(parameters.hpc, 'wait_for_job') && parameters.hpc.wait_for_job
                fprintf('Waiting for calibration job %d to complete...\n', job_id);
                if isfield(parameters.hpc, 'max_wait_checks') && ~isempty(parameters.hpc.max_wait_checks)
                    hpc_wait_for_completion(job_id, platform, parameters.hpc.max_wait_checks);
                else
                    hpc_wait_for_completion(job_id, platform);
                end
            else
                fprintf('Calibration job submitted (ID: %d). Continuing.\n', job_id);
            end

        % ------------------------------------------------------------------
        otherwise
            error('calibration_pipeline_start: unknown platform ''%s''. Use ''matlab'', ''slurm'', ''qsub'', or ''auto''.', platform);
    end
end
