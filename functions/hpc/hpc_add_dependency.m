function parameters = hpc_add_dependency(parameters, job_ids)
% HPC_ADD_DEPENDENCY  Add scheduler job IDs to the jobs a submission waits for
%
% Merges job_ids with any hpc.depend_job_id / hpc.depend_job_ids already set
% and stores the union in hpc.depend_job_ids (which HPC_SUBMIT_JOB prefers).
% Lets a chained pipeline nested inside another one (e.g. multi-ISPPA inside
% an uncertainty variant) keep its own stage dependencies while also waiting
% for the outer pipeline's upstream job.
%
% Use as:
%   parameters = hpc_add_dependency(parameters, job_ids)
%
% Input:
%   parameters - PRESTUS config struct
%   job_ids    - numeric vector of job IDs; [] and non-positive IDs are ignored
%
% Output:
%   parameters - with hpc.depend_job_ids set (hpc.depend_job_id removed)
%
% See also: HPC_SUBMIT_JOB, UNCERTAINTY_PIPELINE, MULTI_ISPPA_PIPELINE

    ids = job_ids(:)';
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'depend_job_ids')
        ids = [parameters.hpc.depend_job_ids(:)', ids];
    end
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'depend_job_id')
        ids = [parameters.hpc.depend_job_id(:)', ids];
        parameters.hpc = rmfield(parameters.hpc, 'depend_job_id');
    end
    ids = unique(ids(ids > 0), 'stable');
    if ~isempty(ids)
        parameters.hpc.depend_job_ids = ids;
    end
end
