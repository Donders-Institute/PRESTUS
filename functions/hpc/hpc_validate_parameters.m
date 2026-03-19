function hpc_validate_parameters(parameters, hpc_type)
%% HPC_VALIDATE_PARAMETERS  Validate HPC job parameters
%
%   hpc_validate_parameters(parameters, hpc_type)
%
%   Disables interactive mode and validates overwrite_files option for batch
%   HPC jobs.
%
%   Inputs:
%     parameters  - Struct with HPC parameters
%     hpc_type    - 'slurm' or 'qsub'
%
%   See also HPC_DETECT_SYSTEM, HPC_SUBMIT_JOB.

if parameters.simulation.interactive
    warning('Interactive mode disabled for %s jobs.', upper(hpc_type));
    parameters.simulation.interactive = false;
end
assert(matches(parameters.io.overwrite_files, ["always", "never"]), ...
    'overwrite_files must be "always" or "never" for %s jobs.', upper(hpc_type));
end
