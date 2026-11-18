function pct_create_pCT_run(parameters)
% PCT_CREATE_PCT_RUN  Submit a pseudoCT generation job
%
% Sources pct_create_pseudoCT.sh and calls the pct_create_pseudoCT
% function for the given subject. Dispatches the job either directly in
% MATLAB or via an HPC batch scheduler (Slurm or PBS/qsub) depending on
% parameters.platform.
%
% SimNIBS charm must have completed before calling this function (the UTE
% image must have been passed as the T2 input so that T2_reg.nii.gz, or the
% UTE_reg.nii.gz older versions renamed it to, is present in the m2m folder;
% the m2m folder is only read).
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
%     .path.seg                 - base SimNIBS output directory (only read)
%     .path.pct                 - (optional) base directory for pCT outputs (GET_PCT_DIR):
%                                 path.pct/sub-NNN/, or <parent of path.seg>/pct/sub-NNN/
%                                 when empty
%     .pct.skull_mapping        - UTE→HU algorithm (default: 'kosciessa'); 'custom' uses
%                                 .pct.custom_slope / .pct.custom_intercept (HU = slope·UTE + intercept)
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

    % Subject-specific pCT output directory (GET_PCT_DIR); the SimNIBS m2m
    % folder is only read.
    path_pct_out = get_pct_dir(parameters);
    if ~isfolder(path_pct_out)
        mkdir(path_pct_out);
    end

    % Resolve the skull mapping (also recorded in the sidecar, PCT_SIGNATURE)
    [pct_sig, skull_mapping] = pct_signature(parameters);

    debug_flag = '1';
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'debug')
        debug_flag = num2str(parameters.pct.debug);
    end

    % Path to the PRESTUS functions directory (needed inside the bash script)
    path_fun = fileparts(fileparts(mfilename('fullpath'))); % functions/

    % Path to the bash script
    script_path = fullfile(path_fun, 'pct', 'pct_create_pseudoCT.sh');

    % Log directory (with the pCT outputs, not in the segmentation folder)
    log_dir = fullfile(path_pct_out, 'log_hpc');
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
               'Set it in your config (e.g. startup.matlab_bin_path: ''/opt/matlab/R2023b/bin/matlab'').']);
    end

    % ── Build environment preamble ────────────────────────────────────────
    % Three mechanisms:
    %   1. startup.modules_to_load  — cell array of module names, loaded via
    %      "module load" (preferred on HPC with environment-module systems).
    %   2. startup.*_bin_path       — directories prepended to PATH directly
    %      (used only when modules_to_load is empty).
    %   3. startup.ants_simg_path   — Singularity image; sets ANTS_SIMG so
    %      the bash script calls ANTs tools via "singularity exec".

    % module load lines (written verbatim into job scripts; inlined for local)
    modules_to_load = {};
    if isfield(parameters, 'startup') && isfield(parameters.startup, 'modules_to_load') ...
            && ~isempty(parameters.startup.modules_to_load)
        modules_to_load = cellstr(parameters.startup.modules_to_load);
    end
    % Inline form used for local bash -c execution
    if ~isempty(modules_to_load)
        module_prefix = ['source /etc/profile && module load ' ...
            strjoin(cellfun(@shell_quote_arg, modules_to_load, 'UniformOutput', false), ' ') ' && '];
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
        path_prefix = sprintf('export PATH=%s:"$PATH" && ', shell_quote_arg(strjoin(extra_paths, ':')));
    else
        path_prefix = '';
    end

    % Singularity image for ANTs (sets ANTS_SIMG for pct_create_pseudoCT.sh)
    ants_simg_export = '';
    if isfield(parameters, 'startup') && isfield(parameters.startup, 'ants_simg_path') ...
            && ~isempty(parameters.startup.ants_simg_path)
        ants_simg_export = sprintf('export ANTS_SIMG=%s && ', shell_quote_arg(parameters.startup.ants_simg_path));
    end

    % ── Build the core bash call (no PATH prefix — each platform sets PATH in its own env block)
    % Arguments: sub_id  path_simnibs  path_matlab  skullmapping  debug  path_fun  path_pct_out
    % Every value is shell_quote_arg()-quoted rather than "%s"-wrapped: a
    % path containing e.g. $(curl ...|sh) would still be expanded by the
    % shell inside double quotes, but not inside shell_quote_arg's single
    % quotes.
    pct_call = sprintf('source %s && pct_create_pseudoCT %s %s %s %s %s %s %s', ...
        shell_quote_arg(script_path), shell_quote_arg(subj_id_string), shell_quote_arg(path_seg), ...
        shell_quote_arg(matlab_bin), shell_quote_arg(skull_mapping), shell_quote_arg(debug_flag), ...
        shell_quote_arg(path_fun), shell_quote_arg(path_pct_out));

    % ── Deploy ────────────────────────────────────────────────────────────
    if strcmp(platform, 'matlab')
        % For local execution prepend dependency paths inline so they are
        % visible to fslmaths / N4BiasFieldCorrection / SmoothImage called
        % from within pct_create_pseudoCT.sh.
        full_call = [module_prefix path_prefix ants_simg_export pct_call];
        fprintf('Running pseudoCT generation locally:\n%s\n', full_call);
        % run_shell_script (not system('bash -c ...')) launches bash
        % directly via ProcessBuilder as a two-element argv, so full_call
        % reaches it as one atomic argument -- never re-tokenized by an
        % outer system()-level shell first.
        [res, out] = run_shell_script(full_call);
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
            fprintf(fid, 'module load %s\n', shell_quote_arg(modules_to_load{mi}));
        end
        if ~isempty(extra_paths)
            fprintf(fid, 'export PATH=%s:"$PATH"\n', shell_quote_arg(strjoin(extra_paths, ':')));
        end
        if ~isempty(ants_simg_export)
            fprintf(fid, 'export ANTS_SIMG=%s\n', shell_quote_arg(parameters.startup.ants_simg_path));
        end
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') ...
                && ~isempty(parameters.hpc.ld_library_path)
            fprintf(fid, 'export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n', parameters.hpc.ld_library_path);
        end
        fprintf(fid, '%s\n', pct_call);
        fclose(fid);

        java.io.File([temp_slurm_file '.sh']).setExecutable(true);
        fprintf('Submitting pseudoCT job:\nsbatch --export=MANPATH %s.sh\n', temp_slurm_file);
        [res, out] = run_command({'sbatch', '--parsable', '--export=MANPATH', [temp_slurm_file '.sh']});
        if res ~= 0
            error('sbatch submission failed (exit code %d):\n%s', res, out);
        end
        job_id_tok = parse_sbatch_job_id(out);
        if isnan(job_id_tok); job_id_tok = []; end
        max_checks = 540;
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'max_wait_checks')
            max_checks = parameters.hpc.max_wait_checks;
        end
        pseudoCT_file = fullfile(path_pct_out, 'pseudoCT.nii.gz');
        if ~isempty(job_id_tok)
            pct_job_id = job_id_tok;
            fprintf('pseudoCT SLURM job %d submitted — waiting for completion.\n', pct_job_id);
            hpc_wait_for_completion(pct_job_id, 'slurm', max_checks);
        else
            % Fallback: no job ID parseable — poll for the output file directly
            warn('Could not parse SLURM job ID from sbatch output:\n%s\nPolling for pseudoCT.nii.gz instead.', out);
            pct_poll_for_file(pseudoCT_file, max_checks);
        end
        % Final guard: error if the file still does not exist after waiting
        if ~exist(pseudoCT_file, 'file')
            error('pseudoCT generation did not produce %s — check the %s job log.', pseudoCT_file, job_name);
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
            fprintf(fid, 'module load %s\n', shell_quote_arg(modules_to_load{mi}));
        end
        if ~isempty(extra_paths)
            fprintf(fid, 'export PATH=%s:"$PATH"\n', shell_quote_arg(strjoin(extra_paths, ':')));
        end
        if ~isempty(ants_simg_export)
            fprintf(fid, 'export ANTS_SIMG=%s\n', shell_quote_arg(parameters.startup.ants_simg_path));
        end
        if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') ...
                && ~isempty(parameters.hpc.ld_library_path)
            fprintf(fid, 'export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n', parameters.hpc.ld_library_path);
        end
        fprintf(fid, 'cd %s\n', shell_quote_arg(path_seg));
        fprintf(fid, '%s\n', pct_call);
        fclose(fid);

        java.io.File([temp_qsub_file '.sh']).setExecutable(true);
        fprintf('Submitting pseudoCT job:\nqsub -v MANPATH %s.sh\n', temp_qsub_file);
        [res, out] = run_command({'qsub', '-v', 'MANPATH', [temp_qsub_file '.sh']});
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
            warn('Empty qsub job ID — polling for pseudoCT.nii.gz instead.');
            pct_poll_for_file(pseudoCT_file, max_checks);
        end
        if ~exist(pseudoCT_file, 'file')
            error('pseudoCT generation did not produce %s — check the pct job log.', pseudoCT_file);
        end

    else
        fprintf(['To generate a pseudoCT, run the following command in bash:\n' ...
                 'bash -c %s\n' ...
                 'Rerun the pipeline once pseudoCT.nii.gz exists in the pCT folder.\n'], ...
                 shell_quote_arg(pct_call));
        error('Platform ''%s'' is not supported for automatic pseudoCT generation.', platform);
    end

    % Sidecar next to pseudoCT.nii.gz recording how it was made; a later run
    % that reuses this pseudo-CT checks its signature (PREPROC_SEGMENTATION).
    write_pct_sidecar(path_pct_out, parameters, pct_sig);

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
    warn('File %s not found after %d checks.', filepath, max_checks);
end

% ── Helper: write pseudoCT.json sidecar with generation provenance ────────
function write_pct_sidecar(pct_dir, parameters, pct_sig)
    m2m = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
    content = struct( ...
        'files', struct('pseudoCT', 'pseudoCT.nii.gz', 'tissues_mask', 'tissues_mask.nii.gz'), ...
        'units', 'HU', ...
        'space', 'T1w (SimNIBS m2m T1.nii.gz grid)', ...
        'inputs', struct('segmentation', m2m), ...
        'skull_mapping', pct_sig.skull_mapping, ...
        'signature', pct_sig);
    sidecar_path = write_sidecar(fullfile(pct_dir, 'pseudoCT.json'), 'pseudoCT', content, parameters);
    if ~isempty(sidecar_path)
        fprintf('pseudoCT.json written to %s\n', sidecar_path);
    end
end
