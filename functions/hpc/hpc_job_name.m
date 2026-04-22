function job_name = hpc_job_name(parameters)
% HPC_JOB_NAME  Generate a standardised HPC job name for a PRESTUS pipeline run
%
% Constructs a job name in the format PREFIX_sub-NNN[_affix], where PREFIX
% defaults to 'PRESTUS' (overridable via parameters.hpc.job_prefix) and affix
% is appended from io.output_affix when set (e.g. '_liberal' → 'liberal').
% If parameters.hpc.job_name is already set, it is returned unchanged.
%
% Use as:
%   job_name = hpc_job_name(parameters)
%
% Input:
%   parameters - PRESTUS config struct; relevant fields:
%                  .subject_id        — used to form sub-NNN suffix
%                  .hpc.job_name      — if set, returned as-is (bypass)
%                  .hpc.job_prefix    — prefix string (default: 'PRESTUS')
%                  .io.output_affix   — appended as disambiguator (optional)
%
% Output:
%   job_name - char; scheduler job name (keep under 20 chars for compatibility)
%
% See also: HPC_SUBMIT_JOB, PRESTUS_PIPELINE_START

    % Allow callers to supply a fully-formed job name (e.g. uncertainty pipeline).
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'job_name') && ...
            ~isempty(parameters.hpc.job_name)
        job_name = parameters.hpc.job_name;
        return;
    end

    subj_id_string = sprintf('sub-%03d', parameters.subject_id);

    if ~isfield(parameters.hpc, 'job_prefix')
        parameters.hpc.job_prefix = 'PS';
    end
    % YAML loading may return strings as cell arrays — coerce to char
    prefix = char(parameters.hpc.job_prefix);

    job_name = [prefix '_' subj_id_string];

    % Append output affix when present so that uncertainty pipeline variants
    % (default / _liberal / _conservative) produce distinguishable job names.
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
            ~isempty(parameters.io.output_affix)
        % Strip leading underscore for readability: '_liberal' -> 'liberal'
        affix = char(parameters.io.output_affix);
        affix = regexprep(affix, '^_', '');
        job_name = [job_name '_' affix];
    end

    % SLURM silently truncates beyond 20 chars; enforce here so the printed
    % name matches what the scheduler stores.
    if numel(job_name) > 20
        job_name = job_name(1:20);
    end
end
