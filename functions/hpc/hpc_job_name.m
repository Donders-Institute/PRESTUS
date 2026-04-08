function job_name = hpc_job_name(parameters)
%% HPC_JOB_NAME  Generate standardized HPC job name
%
%   job_name = hpc_job_name(parameters)
%
%   Creates job name in format: PREFIX_sub-XXX where PREFIX is
%   job_prefix from parameters (default: PRESTUS).
%   Subject ID is read from parameters.subject_id.
%
%   Inputs:
%     parameters    - Struct; must contain parameters.subject_id
%
%   Output:
%     job_name    - String for scheduler job name (max 20 chars recommended)
%
%   See also PRESTUS_PIPELINE_START, HPC_SUBMIT_JOB.

    % Allow callers to supply a fully-formed job name (e.g. uncertainty pipeline).
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'job_name') && ...
            ~isempty(parameters.hpc.job_name)
        job_name = parameters.hpc.job_name;
        return;
    end

    subj_id_string = sprintf('sub-%03d', parameters.subject_id);

    if ~isfield(parameters.hpc, 'job_prefix')
        parameters.hpc.job_prefix = 'PRESTUS';
    end

    job_name = [parameters.hpc.job_prefix '_' subj_id_string];

    % Append output affix when present so that uncertainty pipeline variants
    % (default / _liberal / _conservative) produce distinguishable job names.
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
            ~isempty(parameters.io.output_affix)
        % Strip leading underscore for readability: '_liberal' -> 'liberal'
        affix = regexprep(parameters.io.output_affix, '^_', '');
        job_name = [job_name '_' affix];
    end
end
