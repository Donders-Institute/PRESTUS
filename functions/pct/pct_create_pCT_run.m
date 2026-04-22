function pct_create_pCT_run(parameters)
% PCT_CREATE_PCT_RUN  Submit a pseudoCT generation job
%
% Sources pct_create_pseudoCT.sh and calls the pct_create_pseudoCT
% function for the given subject. Dispatches the job either directly in
% MATLAB or via an HPC batch scheduler (Slurm or PBS/qsub) depending on
% parameters.platform.
%
% SimNIBS charm must have completed before calling this function (the UTE
% image must have been passed as the T2 input so that UTE_reg.nii.gz or
% T2_reg.nii.gz is present in the m2m folder).
%
% The pseudoCT pipeline requires FSL (fslmaths), ANTs (N4BiasFieldCorrection,
% SmoothImage), and MATLAB on PATH. Set startup.fsl_bin_path and
% startup.ants_bin_path in your config so the pipeline can prepend these
% directories to PATH at runtime. On HPC, also set startup.matlab_bin_path.
%
% Use as:
%   pct_create_pCT_run(parameters)
%
% Input:
%   parameters - (1,1) simulation configuration struct with fields:
%     .subject_id               - numeric subject identifier
%     .path.seg                 - base SimNIBS output directory
%     .path.pct                 - (optional) base directory for pCT outputs; empty →
%                                 path.seg/m2m_sub-NNN/, non-empty → path.pct/sub-NNN/
%     .io.pct_dir               - (optional) pre-resolved subject pCT dir set by path_log_setup
%     .pct.skull_mapping        - UTE→HU algorithm (default: 'kosciessa')
%     .pct.debug                - keep intermediate images (default: 1)
%     .platform                 - 'matlab'/'slurm'/'qsub'/'auto'
%     .startup.matlab_bin_path  - MATLAB binary (required on HPC; local falls back to matlabroot)
%     .startup.fsl_bin_path     - FSL bin directory (prepended to PATH)
%     .startup.ants_bin_path    - ANTs bin directory (prepended to PATH)
%     .startup.simnibs_bin_path - SimNIBS bin dir (prepended to PATH on HPC)
%     .hpc.ld_library_path      - (optional) LD_LIBRARY_PATH override for HPC
%
% See also: PREPROC_SEGMENTATION, PCT_CREATE_PSEUDOCT

arguments
    parameters (1,1) struct
end

    subj_id_string = sprintf('%03d', parameters.subject_id);
    path_seg       = parameters.path.seg;

    % Resolve the subject-specific pCT output directory.
    % Use io.pct_dir when path_log_setup has already resolved it; otherwise apply the
    % same resolution logic: empty path.pct → m2m_sub-NNN, non-empty → sub-NNN.
    if isfield(parameters, 'io') && isfield(parameters.io, 'pct_dir') && ~isempty(parameters.io.pct_dir)
        path_pct_out = parameters.io.pct_dir;
    elseif isfield(parameters.path, 'pct') && ~isempty(parameters.path.pct)
        path_pct_out = fullfile(parameters.path.pct, sprintf('sub-%s', subj_id_string));
    else
        path_pct_out = fullfile(path_seg, sprintf('m2m_sub-%s', subj_id_string));
    end
    if ~isfolder(path_pct_out)
        mkdir(path_pct_out);
    end

    % Resolve optional pCT parameters with safe defaults
    skull_mapping = 'kosciessa';
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'skull_mapping') ...
            && ~isempty(parameters.pct.skull_mapping)
        skull_mapping = parameters.pct.skull_mapping;
    end

    debug_flag = '1';
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'debug')
        debug_flag = num2str(parameters.pct.debug);
    end

    % Path to the PRESTUS functions directory (needed inside the bash script)
    path_fun = fileparts(fileparts(mfilename('fullpath'))); % functions/

    % Path to the bash script
    script_path = fullfile(path_fun, 'pct', 'pct_create_pseudoCT.sh');

    % Log directory
    log_dir = fullfile(path_seg, 'batch_job_logs');
    if ~isfolder(log_dir)
        mkdir(log_dir);
    end

    % Detect platform
    if strcmp(parameters.platform, 'auto')
        platform = hpc_detect_system();
        parameters.platform = platform;
        fprintf('➤ auto-detected: %s\n', upper(platform));
    else
        platform = parameters.platform;
        fprintf('➤ deploying: %s\n', upper(platform));
    end

    % ── Resolve MATLAB binary ─────────────────────────────────────────────
    if isfield(parameters, 'startup') && isfield(parameters.startup, 'matlab_bin_path') ...
            && ~isempty(parameters.startup.matlab_bin_path)
        matlab_bin = parameters.startup.matlab_bin_path;
    elseif strcmp(platform, 'matlab')
        matlab_bin = fullfile(matlabroot, 'bin', 'matlab');
    else
        error(['parameters.startup.matlab_bin_path must be set for HPC pseudoCT generation. ' ...
               'Set it in your config (e.g. startup.matlab_bin_path: ''/opt/matlab/R2024a/bin/matlab'').']);
    end

    % ── Build environment preamble ────────────────────────────────────────
    % Two mechanisms are supported and can be combined:
    %   1. startup.modules_to_load  — cell array of module names, loaded via
    %      "module load" (preferred on HPC with environment-module systems).
    %   2. startup.*_bin_path       — directories prepended to PATH directly
    %      (fallback for systems without modules or fixed install paths).
    %   3. startup.ants_simg_path   — Singularity image; sets ANTS_SIMG so
    %      the bash script calls ANTs tools via "singularity exec".

    % module load lines (written verbatim into job scripts; inlined for local)
    modules_to_load = {};
    if isfield(parameters, 'startup') && isfield(parameters.startup, 'modules_to_load') ...
            && ~isempty(parameters.startup.modules_to_load)
        modules_to_load = parameters.startup.modules_to_load;
        if ischar(modules_to_load)
            modules_to_load = {modules_to_load};
        end
    end
    % Inline form used for local bash -c execution
    if ~isempty(modules_to_load)
        module_prefix = ['source /etc/profile && module load ' ...
            strjoin(modules_to_load, ' ') ' && '];
    else
        module_prefix = '';
    end

    % PATH-based fallback (only used when modules_to_load is empty)
    extra_paths = {};
    if isempty(modules_to_load) && isfield(parameters, 'startup')
        for field = {'fsl_bin_path', 'ants_bin_path', 'simnibs_bin_path'}
            f = field{1};
            if isfield(parameters.startup, f) && ~isempty(parameters.startup.(f))
                extra_paths{end+1} = parameters.startup.(f); %#ok<AGROW>
            end
        end
    end
    if ~isempty(extra_paths)
        path_prefix = sprintf('export PATH=%s:$PATH && ', strjoin(extra_paths, ':'));
    else
        path_prefix = '';
    end

    % Singularity image for ANTs (sets ANTS_SIMG for pct_create_pseudoCT.sh)
    ants_simg_export = '';
    if isfield(parameters, 'startup') && isfield(parameters.startup, 'ants_simg_path') ...
            && ~isempty(parameters.startup.ants_simg_path)
        ants_simg_export = sprintf('export ANTS_SIMG="%s" && ', parameters.startup.ants_simg_path);
    end

    % ── Build the core bash call (no PATH prefix — each platform sets PATH in its own env block)
    % Arguments: sub_id  path_simnibs  path_matlab  skullmapping  debug  path_fun  path_pct_out
    pct_call = sprintf('source "%s" && pct_create_pseudoCT "%s" "%s" "%s" "%s" "%s" "%s" "%s"', ...
        script_path, subj_id_string, path_seg, matlab_bin, ...
        skull_mapping, debug_flag, path_fun, path_pct_out);

    % ── Deploy ────────────────────────────────────────────────────────────
    if strcmp(platform, 'matlab')
        % For local execution prepend dependency paths inline so they are
        % visible to fslmaths / N4BiasFieldCorrection / SmoothImage called
        % from within pct_create_pseudoCT.sh.
        full_call = [module_prefix path_prefix ants_simg_export pct_call];
        fprintf('Running pseudoCT generation locally:\n%s\n', full_call);
        [res, out] = system(sprintf('bash -c %s', shellescape(full_call)));
        if res ~= 0
            error('pseudoCT generation failed (exit code %d):\n%s', res, out);
        end
        fprintf('%s\n', out);
        fprintf('pseudoCT generation completed successfully.\n');

    elseif strcmp(platform, 'slurm')
        temp_slurm_file = tempname(log_dir);
        job_name = sprintf('pct-%s', subj_id_string);

        fid = fopen([temp_slurm_file '.sh'], 'w+');
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
        fprintf(fid, '#SBATCH --nodes=1\n');
        fprintf(fid, '#SBATCH --ntasks=1\n');
        fprintf(fid, '#SBATCH --mem=8G\n');
        fprintf(fid, '#SBATCH --time=02:00:00\n');
        fprintf(fid, '#SBATCH --output=%s\n', ...
            fullfile(log_dir, sprintf('sub-%s_slurm_pct_output_%%j.log', subj_id_string)));
        fprintf(fid, '#SBATCH --error=%s\n', ...
            fullfile(log_dir, sprintf('sub-%s_slurm_pct_error_%%j.log', subj_id_string)));
        fprintf(fid, 'source /etc/profile\n');
        for mi = 1:numel(modules_to_load)
            fprintf(fid, 'module load %s\n', modules_to_load{mi});
        end
        if ~isempty(extra_paths)
            fprintf(fid, 'export PATH=%s:$PATH\n', strjoin(extra_paths, ':'));
        end
        if ~isempty(ants_simg_export)
            fprintf(fid, 'export ANTS_SIMG="%s"\n', parameters.startup.ants_simg_path);
        end
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') ...
                && ~isempty(parameters.hpc.ld_library_path)
            fprintf(fid, 'export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n', parameters.hpc.ld_library_path);
        end
        fprintf(fid, '%s\n', pct_call);
        fclose(fid);

        system(sprintf('chmod +x %s.sh', temp_slurm_file));
        sbatch_call = sprintf('sbatch --export=MANPATH %s.sh', temp_slurm_file);
        fprintf('Submitting pseudoCT job:\n%s\n', sbatch_call);
        [res, out] = system(sbatch_call);
        if res ~= 0
            error('sbatch submission failed (exit code %d):\n%s', res, out);
        end
        % Extract first run of digits — matches any sbatch output format
        job_id_tok = regexp(strtrim(out), '\d+', 'match', 'once');
        max_checks = 540;
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'max_wait_checks')
            max_checks = parameters.hpc.max_wait_checks;
        end
        pseudoCT_file = fullfile(path_pct_out, 'pseudoCT.nii.gz');
        if ~isempty(job_id_tok)
            pct_job_id = str2double(job_id_tok);
            fprintf('pseudoCT SLURM job %d submitted — waiting for completion.\n', pct_job_id);
            hpc_wait_for_completion(pct_job_id, 'slurm', max_checks);
        else
            % Fallback: no job ID parseable — poll for the output file directly
            warning('Could not parse SLURM job ID from sbatch output:\n%s\nPolling for pseudoCT.nii.gz instead.', out);
            pct_poll_for_file(pseudoCT_file, max_checks);
        end
        % Final guard: error if the file still does not exist after waiting
        if ~exist(pseudoCT_file, 'file')
            error('pseudoCT generation did not produce %s — check the pct-009 job log.', pseudoCT_file);
        end

    elseif strcmp(platform, 'qsub')
        % Use a temp script file (same pattern as SLURM) to avoid quoting
        % problems when piping multi-word commands through echo | qsub.
        temp_qsub_file = tempname(log_dir);
        job_name = sprintf('pct-%s', subj_id_string);

        fid = fopen([temp_qsub_file '.sh'], 'w+');
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, '#PBS -N %s\n', job_name);
        fprintf(fid, '#PBS -l nodes=1:ppn=1,mem=8gb,walltime=02:00:00\n');
        fprintf(fid, '#PBS -o %s\n', fullfile(log_dir, sprintf('sub-%s_qsub_pct_output.log', subj_id_string)));
        fprintf(fid, '#PBS -e %s\n', fullfile(log_dir, sprintf('sub-%s_qsub_pct_error.log', subj_id_string)));
        fprintf(fid, 'source /etc/profile\n');
        for mi = 1:numel(modules_to_load)
            fprintf(fid, 'module load %s\n', modules_to_load{mi});
        end
        if ~isempty(extra_paths)
            fprintf(fid, 'export PATH=%s:$PATH\n', strjoin(extra_paths, ':'));
        end
        if ~isempty(ants_simg_export)
            fprintf(fid, 'export ANTS_SIMG="%s"\n', parameters.startup.ants_simg_path);
        end
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') ...
                && ~isempty(parameters.hpc.ld_library_path)
            fprintf(fid, 'export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n', parameters.hpc.ld_library_path);
        end
        fprintf(fid, 'cd %s\n', path_seg);
        fprintf(fid, '%s\n', pct_call);
        fclose(fid);

        system(sprintf('chmod +x %s.sh', temp_qsub_file));
        qsub_submit = sprintf('qsub -v MANPATH %s.sh', temp_qsub_file);
        fprintf('Submitting pseudoCT job:\n%s\n', qsub_submit);
        [res, out] = system(qsub_submit);
        if res ~= 0
            error('qsub submission failed (exit code %d):\n%s', res, out);
        end
        pct_job_id = strtrim(out);
        max_checks = 540;
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'max_wait_checks')
            max_checks = parameters.hpc.max_wait_checks;
        end
        pseudoCT_file = fullfile(path_pct_out, 'pseudoCT.nii.gz');
        if ~isempty(pct_job_id)
            fprintf('pseudoCT qsub job %s submitted — waiting for completion.\n', pct_job_id);
            hpc_wait_for_completion(pct_job_id, 'qsub', max_checks);
        else
            warning('Empty qsub job ID — polling for pseudoCT.nii.gz instead.');
            pct_poll_for_file(pseudoCT_file, max_checks);
        end
        if ~exist(pseudoCT_file, 'file')
            error('pseudoCT generation did not produce %s — check the pct job log.', pseudoCT_file);
        end

    else
        fprintf(['To generate a pseudoCT, run the following command in bash:\n' ...
                 'bash -c %s\n' ...
                 'Rerun the pipeline once pseudoCT.nii.gz exists in the m2m folder.\n'], ...
                 shellescape(pct_call));
        error('Platform ''%s'' is not supported for automatic pseudoCT generation.', platform);
    end

end

% ── Helper: quote a string for safe bash -c argument passing ─────────────
function s = shellescape(str)
    s = ['"' strrep(str, '"', '\"') '"'];
end

% ── Helper: poll for a file to appear (fallback when job ID is unavailable) ──
function pct_poll_for_file(filepath, max_checks)
    fprintf('Polling for %s (up to %d checks, 20 s apart)...\n', filepath, max_checks);
    for k = 1:max_checks
        if exist(filepath, 'file')
            fprintf('✓ %s found after check %d.\n', filepath, k);
            return;
        end
        fprintf('  check %d/%d — not yet present\n', k, max_checks);
        pause(20);
    end
    warning('File %s not found after %d checks.', filepath, max_checks);
end
