function job_name = hpc_job_name(parameters, subject_id)
%% HPC_JOB_NAME  Generate standardized HPC job name
%
%   job_name = hpc_job_name(parameters, subject_id)
%
%   Creates job name in format: PREFIX_sub-XXX where PREFIX is
%   {slurm|qsub}_job_prefix from parameters (default: PRESTUS).
%
%   Inputs:
%     parameters    - Struct (may contain slurm_job_prefix/qsub_job_prefix)
%     subject_id    - Subject number (double)
%
%   Output:
%     job_name    - String for scheduler job name (max 20 chars recommended)
%
%   See also PRESTUS_PIPELINE_START, HPC_SUBMIT_JOB.

    subj_id_string = sprintf('sub-%03d', subject_id);
    prefix_field = sprintf('job_prefix');
    
    if ~isfield(parameters, 'job_prefix')
        parameters.job_prefix = 'PRESTUS';
    end
    
    job_name = [parameters.job_prefix '_' subj_id_string];
end
