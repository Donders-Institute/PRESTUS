function hpc_validate_parameters(parameters, hpc_type)
% HPC_VALIDATE_PARAMETERS  Validate and enforce settings required for HPC batch submission
%
% Checks that parameters are compatible with non-interactive HPC execution.
% Warns and disables interactive mode if enabled. Asserts that overwrite_files
% is set to 'always' or 'never' (not 'ask', which would block an unattended job).
% This function has no return value — it errors or warns rather than returning a flag.
%
% Use as:
%   hpc_validate_parameters(parameters, hpc_type)
%
% Input:
%   parameters - PRESTUS config struct; checked fields:
%                  .simulation.interactive — must be false for batch jobs
%                  .io.overwrite_files     — must be 'always' or 'never'
%   hpc_type   - 'slurm' or 'qsub' (used in warning/error messages)
%
% See also: HPC_DETECT_SYSTEM, HPC_SUBMIT_JOB

if parameters.simulation.interactive
    warning('Interactive mode disabled for %s jobs.', upper(hpc_type));
    parameters.simulation.interactive = false;
end
assert(matches(parameters.io.overwrite_files, ["always", "never"]), ...
    'overwrite_files must be "always" or "never" for %s jobs.', upper(hpc_type));
end
