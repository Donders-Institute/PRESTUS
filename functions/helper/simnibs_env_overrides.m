function env = simnibs_env_overrides(parameters)
% SIMNIBS_ENV_OVERRIDES  Environment for a SimNIBS tool run with RUN_COMMAND
%
% Prepends SIMNIBS_LD_LIBRARY_PATH (hpc.ld_library_path, else the libraries of
% the SimNIBS install in startup.simnibs_bin_path) to LD_LIBRARY_PATH, keeping
% the current value behind it, as the SLURM job scripts do. Empty struct when
% there is nothing to prepend (and on Windows).
%
% Use as:
%   env = simnibs_env_overrides(parameters)
%   run_command(argv, 'Env', env)
%
% Input:
%   parameters - PRESTUS config
%
% Output:
%   env - struct for run_command's 'Env' option
%
% See also: SIMNIBS_LD_LIBRARY_PATH, RUN_COMMAND

    env = struct();
    ld = simnibs_ld_library_path(parameters);
    if isempty(ld) || ispc
        return
    end
    current = getenv('LD_LIBRARY_PATH');
    if isempty(current)
        env.LD_LIBRARY_PATH = ld;
    else
        env.LD_LIBRARY_PATH = [ld ':' current];
    end
end
