function [parameters, cleanup] = hpc_scratch_setup(parameters)
% HPC_SCRATCH_SETUP  Put the k-Wave .h5 files of a SLURM job on node-local scratch
%
% With hpc.tmp_gb > 0 (the job requested node-local disk, HPC_SUBMIT_JOB) and
% MATLAB running inside a SLURM job, creates a job-specific folder on scratch
% and sets io.dir_scratch. The C++ k-Wave codes (cpp_cpu, cpp_gpu) write their
% input and output .h5 files there (ACOUSTIC_SIMULATION). The folder is
% removed when the pipeline returns or errors. cache/ stays in the output
% folder: it holds results that later stages and runs read.
%
% Scratch root (hpc.scratch_dir): 'auto' takes the first existing of
% $PRESTUS_SCRATCH, $SLURM_TMPDIR, <base>/$USER/$SLURM_JOB_NAME and
% <base>/$USER/$SLURM_JOB_ID; if none exists but <base>/$USER does, it creates
% <base>/$USER/$SLURM_JOB_NAME (removed again at the end if empty); last
% $TMPDIR. <base> is $PRESTUS_SCRATCH_BASE, default /scratch. Any other value
% is a path; $VAR and ${VAR} are expanded and the folder is created if missing.
%
% Outside a SLURM job, with hpc.tmp_gb = 0, or when no root is found,
% parameters are returned unchanged (the .h5 files then go to $TMPDIR or
% cache/, ACOUSTIC_SIMULATION).
%
% Use as:
%   [parameters, cleanup] = hpc_scratch_setup(parameters)
%
% Input:
%   parameters - PRESTUS config after PATH_LOG_SETUP; reads hpc.tmp_gb,
%                hpc.scratch_dir, subject_id, io.output_affix
%
% Output:
%   parameters - io.dir_scratch set when scratch is used
%   cleanup    - onCleanup object that calls HPC_SCRATCH_FINALIZE; [] when
%                scratch is not used or a nested pipeline call in the same
%                job reuses a folder it does not own
%
% See also: HPC_SCRATCH_FINALIZE, HPC_SUBMIT_JOB, ACOUSTIC_SIMULATION

cleanup = [];
tmp_gb = get_nonempty_field(parameters, {'hpc', 'tmp_gb'}, 0);
if ~isnumeric(tmp_gb) || ~isscalar(tmp_gb) || tmp_gb <= 0
    return
end
job_id = getenv('SLURM_JOB_ID');
if isempty(job_id)
    return
end
job_id = regexprep(job_id, '[^\w.-]', '_');

[root, created_root] = resolve_scratch_root(parameters);
if isempty(root)
    warning('prestus:scratch', ['hpc.tmp_gb = %g but no scratch folder was found (hpc.scratch_dir, ' ...
        '$PRESTUS_SCRATCH, $SLURM_TMPDIR, /scratch/$USER, $TMPDIR); the .h5 files stay in the default place.'], tmp_gb);
    return
end
if ~isempty(created_root)
    fprintf('Created per-job scratch folder: %s\n', created_root);
end

affix = char(get_nonempty_field(parameters, {'io', 'output_affix'}, ''));
scratch_dir = fullfile(root, sprintf('prestus_%s', job_id), sprintf('sub-%03d%s', parameters.subject_id, affix));

% A nested pipeline call in the same job reuses the folder but does not own
% it, so it never deletes it from under the outer run.
owns_dir = ~isfolder(scratch_dir);
if owns_dir
    [ok, msg] = mkdir(scratch_dir);
    if ~ok
        warning('prestus:scratch', 'Could not create scratch folder %s (%s); the .h5 files stay in the default place.', ...
            scratch_dir, msg);
        if ~isempty(created_root), [~, ~] = rmdir(created_root); end
        return
    end
end
parameters.io.dir_scratch = scratch_dir;

free_gb = scratch_free_gb(scratch_dir);
if isnan(free_gb)
    fprintf('Node-local scratch: %s\n', scratch_dir);
else
    fprintf('Node-local scratch: %s (%.1f GB free, %g GB requested)\n', scratch_dir, free_gb, tmp_gb);
    if free_gb < tmp_gb
        warning('prestus:scratch', 'Only %.1f GB free on scratch, less than hpc.tmp_gb = %g GB.', free_gb, tmp_gb);
    end
end

if owns_dir
    cleanup = onCleanup(@() hpc_scratch_finalize(scratch_dir, created_root));
end
end

% ========== LOCAL FUNCTIONS ==========
function [root, created] = resolve_scratch_root(parameters)
    root = '';
    created = '';
    spec = strtrim(char(get_nonempty_field(parameters, {'hpc', 'scratch_dir'}, 'auto')));
    if ~strcmpi(spec, 'auto')
        cand = expand_env(spec);
        if ~isempty(cand) && (isfolder(cand) || mkdir(cand))
            root = cand;
        end
        return
    end

    base = getenv('PRESTUS_SCRATCH_BASE');
    if isempty(base), base = '/scratch'; end
    user_dir = [base '/$USER'];
    job_dir  = [user_dir '/$SLURM_JOB_NAME'];
    for c = {'$PRESTUS_SCRATCH', '$SLURM_TMPDIR', job_dir, [user_dir '/$SLURM_JOB_ID']}
        cand = expand_env(c{1});
        if ~isempty(cand) && isfolder(cand)
            root = cand;
            return
        end
    end
    % The per-job folder does not exist yet: create it under the user's scratch
    cand = expand_env(job_dir);
    user = expand_env(user_dir);
    if ~isempty(cand) && ~isempty(user) && isfolder(user) && mkdir(cand)
        root = cand;
        created = cand;
        return
    end
    cand = expand_env('$TMPDIR');
    if ~isempty(cand) && isfolder(cand)
        root = cand;
    end
end

function s = expand_env(s)
% Expand $VAR and ${VAR}; '' if a referenced variable is unset, so a path like
% /scratch/$USER/$SLURM_JOB_NAME never collapses to /scratch/<user>/.
    vars = regexp(s, '\$\{?(\w+)\}?', 'tokens');
    for k = 1:numel(vars)
        if isempty(getenv(vars{k}{1}))
            s = '';
            return
        end
    end
    s = regexprep(s, '\$\{?(\w+)\}?', '${getenv($1)}');
end

function free_gb = scratch_free_gb(d)
% Free space [GB] from df (POSIX); NaN when unavailable.
    free_gb = NaN;
    if ~isunix
        return
    end
    try
        [st, out] = run_command({'df', '-Pk', d});
        if st == 0
            lines = strsplit(strtrim(out), newline);
            tok = regexp(strtrim(lines{end}), '\s+', 'split');
            free_gb = str2double(tok{4}) / 1024^2;
        end
    catch
    end
end
