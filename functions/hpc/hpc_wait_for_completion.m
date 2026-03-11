function hpc_wait_for_completion(job_id, hpc_type)
%% HPC_WAIT_FOR_COMPLETION  Monitor HPC job until completion
%
%   hpc_wait_for_completion(job_id, hpc_type)
%
%   Polls job status every 20s until COMPLETED/FAILED (SLURM) or C/no-longer-listed (qsub).
%
%   Inputs:
%     job_id    - Job ID (numeric string for SLURM, string for qsub)
%     hpc_type  - 'slurm' or 'qsub'
%
%   See also HPC_SUBMIT_JOB.

disp('User has chosen to wait until job is finished...');
job_completed = false;

while ~job_completed
    switch hpc_type
        case 'slurm'
            [status, out] = system(sprintf('sacct -j %s -o State --noheader | tail -n 1', job_id));
            if status == 0
                job_state = strtrim(out);
                fprintf('SLURM Job status: %s\n', job_state);
                if any(strcmp(job_state, {'COMPLETED', 'FAILED'}))
                    job_completed = true;
                else
                    pause(20);
                end
            else
                pause(20);
            end
            
        case 'qsub'
            [status, out] = system(sprintf('qstat -f %s | grep job_state', job_id));
            if status == 0
                parts = strsplit(out, '=');
                if numel(parts) == 2
                    job_state = strtrim(parts{2});
                    fprintf('qsub Job status: %s\n', job_state);
                    if strcmp(job_state, 'C')
                        disp('Job completed successfully.');
                        job_completed = true;
                    else
                        pause(20);
                    end
                end
            else
                [status_check, ~] = system(sprintf('qstat %s', job_id));
                if status_check ~= 0
                    disp('Job no longer listed. Assuming completed.');
                    job_completed = true;
                end
            end
    end
end
end
