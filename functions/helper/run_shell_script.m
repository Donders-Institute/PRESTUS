function [status, cmdout] = run_shell_script(script, options)
% RUN_SHELL_SCRIPT  Execute a shell script string via ProcessBuilder, with only ONE shell in the loop.
%
% Some workflows genuinely need real shell semantics — `module load`,
% chained `export VAR=... && cmd`, command substitution, filename globbing
% — that run_command.m's plain-argv model can't express. For those cases,
% use this function rather than MATLAB's system(): system(script) hands
% the *entire* string to a shell for parsing, which is an extra, redundant
% shell layer on top of whatever shell you actually intend to invoke — and
% is exactly the mechanism the BabelBrain report (PRs #137/#139) exploited
% (see shell_quote_arg.m for the full background).
%
% run_shell_script instead launches the shell itself directly through
% ProcessBuilder, as a two/three-element argv ({'/bin/bash','-c',script}),
% so `script` reaches bash as ONE atomic argument — never re-tokenized,
% re-quoted, or re-expanded by an outer shell first. It is parsed by
% exactly the one shell you asked for, once.
%
% This does NOT make the script text itself immune to injection: any
% untrusted value interpolated into `script` still needs shell_quote_arg()
% around it, the same as if you were building it for system(). What this
% function removes is the redundant *outer* shell layer, not the need for
% care inside the script you hand it.
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
