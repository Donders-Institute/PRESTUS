function parameters = hpc_validate_parameters(parameters, hpc_type)
% HPC_VALIDATE_PARAMETERS  Validate and enforce settings required for HPC batch submission
%
% Checks that parameters are compatible with non-interactive HPC execution.
% Warns and disables interactive mode if enabled. Asserts that overwrite_files
% is set to 'always', 'force' or 'never' (not 'ask', which would block an unattended job).
%
% Use as:
%   parameters = hpc_validate_parameters(parameters, hpc_type)
%
% Input:
%   parameters - PRESTUS config struct; checked fields:
%                  .simulation.interactive — must be false for batch jobs
%                  .io.overwrite_files     — must be 'always', 'force' or 'never'
%                  .hpc.tmp_gb             — scalar >= 0 (node-local disk [GB])
%   hpc_type   - 'slurm' or 'qsub' (used in warning/error messages)
%
% Output:
%   parameters - input struct with .simulation.interactive set to false
%
% See also: HPC_DETECT_SYSTEM, HPC_SUBMIT_JOB

if parameters.simulation.interactive
    warn('Interactive mode disabled for %s jobs.', upper(hpc_type));
    parameters.simulation.interactive = false;
end
assert(matches(parameters.io.overwrite_files, ["always", "force", "never"]), ...
    'overwrite_files must be "always", "force" or "never" for %s jobs.', upper(hpc_type));
if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'tmp_gb') && ~isempty(parameters.hpc.tmp_gb)
    assert(isnumeric(parameters.hpc.tmp_gb) && isscalar(parameters.hpc.tmp_gb) && parameters.hpc.tmp_gb >= 0, ...
        'hpc.tmp_gb must be a scalar >= 0 (GB of node-local disk per job).');
end
end
