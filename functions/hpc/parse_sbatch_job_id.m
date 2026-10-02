function job_id = parse_sbatch_job_id(out)
% PARSE_SBATCH_JOB_ID  Job ID from the output of sbatch
%
% Submit with 'sbatch --parsable', which prints '<jobid>' or
% '<jobid>;<cluster>' on its last line. Site wrappers and SLURM may print
% warnings before it (e.g. 'sbatch: warning: ... 2 GPUs'), so the first number
% anywhere in the output is not reliable. Reads, in order: the last line in
% --parsable form, then 'Submitted batch job <id>' (sbatch without --parsable).
%
% Use as:
%   job_id = parse_sbatch_job_id(out)
%
% Input:
%   out - char/string output of sbatch
%
% Output:
%   job_id - double job ID; NaN when none is found
%
% See also: HPC_SUBMIT_JOB, HPC_SUBMIT_REPORT_JOB

    job_id = NaN;
    text = strtrim(char(out));
    if isempty(text)
        return
    end
    lines = strtrim(splitlines(text));
    lines = lines(~cellfun(@isempty, lines));
    tok = regexp(lines{end}, '^(\d+)(;\S+)?$', 'tokens', 'once');
    if ~isempty(tok)
        job_id = str2double(tok{1});
        return
    end
    tok = regexp(text, 'Submitted batch job (\d+)', 'tokens', 'once');
    if ~isempty(tok)
        job_id = str2double(tok{1});
    end
end
