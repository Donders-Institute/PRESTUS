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

    subj_id_string = sprintf('sub-%03d', parameters.subject_id);
    
    if ~isfield(parameters.hpc, 'job_prefix')
        parameters.hpc.job_prefix = 'PRESTUS';
    end
    
    job_name = [parameters.hpc.job_prefix '_' subj_id_string];
end
