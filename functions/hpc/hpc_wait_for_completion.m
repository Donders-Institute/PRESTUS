function hpc_wait_for_completion(job_id, hpc_type, max_checks)
% HPC_WAIT_FOR_COMPLETION  Poll an HPC job until it completes or times out
%
% Repeatedly checks job status via squeue/sacct (SLURM) or qstat (qsub)
% every 20 seconds until the job finishes or max_checks is reached.
% Calibrated for the Donders/DCCN SLURM cluster column format.
% Blocks the calling MATLAB session for the duration of polling.
%
% Use as:
%   hpc_wait_for_completion(job_id, hpc_type)
%   hpc_wait_for_completion(job_id, hpc_type, max_checks)
%
% Input:
%   job_id     - scheduler job ID (numeric for SLURM, string for qsub)
%   hpc_type   - 'slurm' or 'qsub'
%   max_checks - (optional) maximum number of status polls before giving up;
%                default 540 (~3 hours at one check per 20 s)
%
% Note:
%   This function blocks the MATLAB session. Only use when the caller
%   explicitly wants to wait (e.g. single-stage pipeline on SLURM).
%   For fire-and-forget submission use hpc_submit_job without calling this.
%
% See also: HPC_SUBMIT_JOB, HPC_DETECT_SYSTEM

if nargin < 3 || isempty(max_checks)
    max_checks = 540;  % default: ~3 hours at 1 check/20s
end
disp('User has chosen to wait until job is finished...');
job_completed = false;
checks = 0;
tic_start = tic;

while ~job_completed && checks < max_checks
    checks = checks + 1;
    
    switch lower(hpc_type)
        case 'slurm'
            job_id_str = sprintf('%.0f', job_id);
            
            % 1. squeue - simple column split for YOUR format
            [status_q, out_q] = system(sprintf('squeue --noheader -j %s 2>/dev/null', job_id_str));
            if status_q == 0 && ~isempty(strtrim(out_q))
                parts = strsplit(strtrim(out_q));
                if length(parts) >= 5
                    job_state = upper(parts{5});  % Column 5 = ST (PD, R, CG...)
                    fprintf('SLURM Job %s: %s (check %d/%d)\n', job_id_str, job_state, checks, max_checks);
                    
                    % Keep waiting for active states
                    if any(strcmp(job_state, {'PD','R','CG','CD','F','TO','S','CA'}))
                        pause(20);
                        continue;
                    end
                end
            end
            
            % 2. Job gone from squeue → check sacct
            [status, out] = system(sprintf('sacct -j %s -o State --noheader | tail -n 1 2>/dev/null', job_id_str));
            if status == 0 && ~isempty(strtrim(out))
                % Safe first word extraction
                first_word = strtok(strtrim(out));  % Safer than regex
                terminal_states = {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY'};
                if any(strcmpi(first_word, terminal_states))
                    fprintf('✓ SLURM Job %s FINISHED: %s\n', job_id_str, first_word);
                    job_completed = true;
                else
                    fprintf('SLURM Job %s sacct: %s (check %d)\n', job_id_str, first_word, checks);
                    pause(20);
                end
            else
                fprintf('✓ SLURM Job %s gone from queue → COMPLETE (check %d)\n', job_id_str, checks);
                job_completed = true;
            end
            
        case {'qsub', 'pbs', 'torque'}
            job_id_str = char(job_id);
            
            % qstat -f for detailed state
            [status, out] = system(sprintf('qstat -f %s 2>/dev/null | grep job_state', job_id_str));
            if status == 0 && ~isempty(out)
                parts = strsplit(strtrim(out), '=');
                if length(parts) >= 2
                    job_state = strtrim(parts{end});
                    fprintf('qsub Job %s: %s (check %d/%d)\n', job_id_str, job_state, checks, max_checks);
                    
                    if any(strcmpi(job_state, {'C','E','F'}))
                        fprintf('✓ qsub Job %s FINISHED: %s\n', job_id_str, job_state);
                        job_completed = true;
                    else
                        pause(20);
                    end
                    continue;
                end
            end
            
            % Fallback: simple qstat
            [status_check, ~] = system(sprintf('qstat %s >/dev/null 2>&1', job_id_str));
            persistent seen_once;
            if isempty(seen_once), seen_once = false; end
            
            if status_check ~= 0 && seen_once
                fprintf('✓ qsub Job %s gone from queue → COMPLETE\n', job_id_str);
                job_completed = true;
            elseif ~seen_once
                fprintf('qsub Waiting for job %s to appear... (check %d/%d)\n', job_id_str, checks, max_checks);
                seen_once = true;
            end
            pause(20);
    end
end

elapsed = toc(tic_start);
if ~job_completed
    warning('TIMEOUT job %s after %.1f min (%d checks)', job_id_str, elapsed/60, checks);
else
    fprintf('✓ Job %s (%s) complete after %.1f min (%d checks)\n', job_id_str, hpc_type, elapsed/60, checks);
end

end

