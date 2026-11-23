function tf = hpc_job_active(platform, job_id)
% HPC_JOB_ACTIVE  True while a scheduler job is queued or running
%
% Asks the scheduler once (squeue for SLURM, qstat for qsub). A job that has
% finished, failed, or is unknown to the scheduler, or a query that fails,
% gives false, so the caller submits the work again rather than waiting on a
% job that will never finish.
%
% Use as:
%   tf = hpc_job_active('slurm', '12345')
%
% Input:
%   platform - 'slurm' or 'qsub'
%   job_id   - job ID (char or numeric)
%
% Output:
%   tf - logical
%
% See also: HPC_WAIT_FOR_COMPLETION, ASYNC_TRANSDUCER_PIPELINE

arguments
    platform {mustBeTextScalar}
    job_id
end

    tf = false;
    if isnumeric(job_id)
        job_id = num2str(job_id);
    end
    job_id = strtrim(char(job_id));
    if isempty(job_id) || ~all(isstrprop(strrep(job_id, '.', ''), 'alphanum'))
        return
    end
    try
        switch platform
            case 'slurm'
                [status, out] = run_command({'squeue', '--noheader', '-o', '%T', '-j', job_id});
                state = strtrim(out);
                tf = status == 0 && any(startsWith(state, ...
                    {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'REQUEUED', 'RESIZING', 'SUSPENDED'}));
            case 'qsub'
                [status, out] = run_command({'qstat', '-f', job_id});
                tf = status == 0 && ~isempty(regexp(out, 'job_state\s*=\s*[QRHTW]', 'once'));
        end
    catch
        tf = false;
    end
end
