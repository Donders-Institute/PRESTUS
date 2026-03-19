function display_info = hpc_job_info(platform, job_id, job_name, ...
    memory_gb, timelimit, log_dir, visualize)
%% HPC_JOB_INFO  Generate formatted job display information
%
%   Creates structured display info for job submission feedback.
%
%   Inputs:
%     platform      - 'slurm' or 'qsub'
%     job_id        - Raw job ID from submission
%     job_name      - Job name string (encodes subject ID)
%     memory_gb     - Memory allocation (GB)
%     timelimit     - Time limit string
%     log_dir       - Log directory path
%     visualize     - Display job info?
%
%   Output:
%     display_info  - Struct with all formatted display fields

    % Format job ID
    if strcmp(platform, 'slurm') && isnumeric(job_id)
        display_info.id = sprintf('%d', job_id);
        display_info.check_cmd = sprintf('squeue -u $USER | grep %s', job_name);
        display_info.detail_cmd = sprintf('sacct -j %s', display_info.id);
    else
        display_info.id = string(job_id);
        display_info.check_cmd = sprintf('qstat %s', display_info.id);
        display_info.detail_cmd = display_info.check_cmd;
    end
    
    % Common fields
    display_info.name = job_name;
    display_info.memory_gb = memory_gb;
    display_info.timelimit = timelimit;
    display_info.log_dir = log_dir;

    if visualize == true
        fprintf('\n⚙️ JOB INFO\n');
        fprintf('═════════════════════════════\n');
        fprintf('Job ID:       %s\n', display_info.id);
        fprintf('Job name:     %s\n', display_info.name);
        fprintf('Subject:      %s\n', job_name);
        fprintf('Memory:       %.0f GB\n', display_info.memory_gb);
        fprintf('Time limit:   %s\n', display_info.timelimit);
        fprintf('Log dir:      %s\n', display_info.log_dir);
        fprintf('Check:        %s\n', display_info.check_cmd);
        fprintf('Details:      %s\n', display_info.detail_cmd);
        fprintf('\n');
    end
end
