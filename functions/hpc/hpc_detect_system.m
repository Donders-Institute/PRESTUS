function hpc_type = hpc_detect_system()
%% HPC_DETECT_SYSTEM  Detect SLURM or qsub HPC system
%
%   Detects available HPC scheduler by checking for sbatch (SLURM) or qsub.
%   Returns 'slurm', 'qsub', or throws error if neither found.
%
%   See also HPC_SUBMIT_JOB, HPC_WAIT_FOR_COMPLETION.

[sbatch_out, sbatch_path] = system('which sbatch');
[qsub_out, qsub_path] = system('which qsub');

if sbatch_out == 0
    hpc_type = 'slurm';
elseif qsub_out == 0
    hpc_type = 'qsub';
else
    error('No HPC system detected (sbatch/qsub)');
end
end
