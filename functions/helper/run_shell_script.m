function [status, cmdout] = run_shell_script(script, options)
% RUN_SHELL_SCRIPT  Execute a shell script string via ProcessBuilder, with only ONE shell in the loop.
%
% For commands that need shell semantics (module load, chained
% `export VAR=... && cmd`, globbing). The script reaches bash -c (cmd.exe /c
% on Windows) as one argument through RUN_COMMAND, with no outer shell;
% quote every value interpolated into it with SHELL_QUOTE_ARG.
%
% Use as:
%   [status, cmdout] = run_shell_script(sprintf( ...
%       'module load %s && mytool %s', module_name, shell_quote_arg(input_path)));
%
% Input:
%   script - (1,:) char, the shell script/command line to execute
%
% Name-value options: same as run_command.m (WorkingDir, Env, ExtraPath,
%   Stdin, Echo, Background) — forwarded as-is.
%
% Output: [status, cmdout] — same semantics as run_command.m / system().
%
% See also: RUN_COMMAND, SHELL_QUOTE_ARG, SYSTEM

    arguments
        script (1,:) char
        options.WorkingDir (1,:) char = ''
        options.Env (1,1) struct = struct()
        options.ExtraPath (1,:) cell = {}
        options.Stdin (1,:) char = ''
        options.Echo (1,1) logical = false
        options.Background (1,1) logical = false
    end

    if ispc
        shell_argv = {'cmd.exe', '/c', script};
    else
        shell_argv = {'/bin/bash', '-c', script};
    end

    nv = namedargs2cell(options);
    [status, cmdout] = run_command(shell_argv, nv{:});
end
