function [status, cmdout] = run_command(argv, options)
% RUN_COMMAND  Execute an external program via java.lang.ProcessBuilder.
%
% Drop-in emulation of MATLAB's [status, cmdout] = system(cmd) / system(cmd,
% '-echo'), but takes the program and its arguments as a cell array (argv)
% instead of a single command string, and NEVER hands that string to a
% shell for parsing.
%
% Each argv element reaches the program as one literal argument
% (execve/CreateProcess, no shell), so no value, e.g. a file name with
% $(...) in it, can inject shell syntax.
%
% Use run_command for any call that is really just "run this program with
% these arguments". For the rarer case that genuinely needs shell
% semantics (module load, chained `export ... &&`, globbing), use
% run_shell_script.m instead, which still avoids the redundant *outer*
% system()-level shell.
%
% Use as:
%   [status, cmdout] = run_command({'charm', 'sub-001', t1_path})
%   [status, cmdout] = run_command({python_bin, script, '--flag', val}, 'Echo', true)
%   [status, cmdout] = run_command({'charm', args{:}}, 'WorkingDir', seg_dir)
%   [status, cmdout] = run_command({'mni2subject_coords', args{:}}, ...
%       'Env', struct('LD_LIBRARY_PATH', ld_path))
%
% Input:
%   argv - (1,:) cell array of char row vectors: {program, arg1, arg2, ...}.
%          Each element is passed to the OS exactly as given — including
%          any spaces, quotes, or shell metacharacters it contains, which
%          are always treated as literal argument content, never syntax.
%
% Name-value options:
%   WorkingDir - (1,:) char, directory the child process runs in
%                (default: MATLAB's current pwd, matching system())
%   Env        - (1,1) struct of environment variables to set/override,
%                merged on top of the inherited environment. This is the
%                native replacement for shell patterns like
%                `export LD_LIBRARY_PATH="..." && cmd`.
%   ExtraPath  - (1,:) cell array of directories prepended to PATH
%                (convenience wrapper around Env.PATH)
%   Stdin      - (1,:) char, text piped to the child process's stdin
%                (replacement for shell patterns like `echo "..." | cmd`)
%   Echo       - (1,1) logical, stream output to the console live as it is
%                produced, matching system(cmd, '-echo') (default: false)
%   Background - (1,1) logical, start the process and return immediately
%                without waiting for it to exit or reading its output
%                (replacement for shell patterns like `cmd &` used to run
%                several jobs concurrently and poll for their output files
%                afterward). Output is redirected to a throwaway temp file
%                (there is nothing to read asynchronously without blocking
%                the caller); the file is intentionally left behind rather
%                than deleted, since there is no safe synchronous point at
%                which the still-running process is guaranteed done with
%                it. Default: false.
%
% Output:
%   status - process exit code (0 = success), or 127 if the program could
%            not be found/started — matching the POSIX shell's own
%            "command not found" convention, which is what system() itself
%            returns in that situation. With Background=true, this is NaN
%            (the exit code is not yet known) rather than a real status.
%   cmdout - combined stdout+stderr text (stderr is merged into stdout,
%            matching system()'s default combined-stream behavior). Line
%            terminators are normalized to the platform's newline; a
%            trailing newline is added whenever the process produced any
%            output, whether or not its last line ended in one — an
%            immaterial difference for the regexp/contains-style parsing
%            PRESTUS does on cmdout, but not a byte-exact reproduction.
%            With Background=true, this is always '' (output is discarded).
%
% See also: SYSTEM, RUN_SHELL_SCRIPT, SHELL_QUOTE_ARG

    arguments
        argv (1,:) cell
        options.WorkingDir (1,:) char = ''
        options.Env (1,1) struct = struct()
        options.ExtraPath (1,:) cell = {}
        options.Stdin (1,:) char = ''
        options.Echo (1,1) logical = false
        options.Background (1,1) logical = false
    end

    if isempty(argv)
        error('run_command:emptyArgv', 'argv must contain at least the program name.');
    end
    argv = cellfun(@(a) char(a), argv, 'UniformOutput', false);

    % Effective PATH for THIS call: Env.PATH if given, else MATLAB's live
    % getenv('PATH'); ExtraPath is prepended on top either way.
    if isfield(options.Env, 'PATH')
        effective_path = char(string(options.Env.PATH));
    else
        effective_path = getenv('PATH');
    end
    if ~isempty(options.ExtraPath)
        extra = cellfun(@char, options.ExtraPath, 'UniformOutput', false);
        effective_path = strjoin([extra, {effective_path}], char(java.io.File.pathSeparator));
    end

    % Resolve a bare program name (argv{1} with no path separator) to an
    % absolute path ourselves, using the effective PATH above, BEFORE
    % constructing ProcessBuilder. This is not optional politeness: Java's
    % ProcessBuilder resolves a bare executable name using the JVM's OWN
    % process environment (a snapshot from JVM startup), never the
    % environment map set via pb.environment() below -- that map only
    % controls what the CHILD process sees once launched, not how the
    % executable itself is located. Confirmed empirically: mutating
    % pb.environment().PATH has no effect on which binary gets found.
    % Without this, run_command could never find a tool whose directory
    % was added to PATH via setenv() after MATLAB/the JVM started -- which
    % is exactly how PRESTUS's own tests mock external tools, and a
    % realistic pattern in real module-load-style HPC setups too.
    argv{1} = resolve_executable(argv{1}, effective_path);

    pb = java.lang.ProcessBuilder(argv);
    pb.redirectErrorStream(true);   % merge stderr into stdout, matching system()'s combined output

    if ~isempty(options.WorkingDir)
        pb.directory(java.io.File(options.WorkingDir));
    end

    env_map = pb.environment();
    env_fields = fieldnames(options.Env);
    for i = 1:numel(env_fields)
        env_map.put(env_fields{i}, char(string(options.Env.(env_fields{i}))));
    end

    % Also sync the CHILD's own view of PATH (distinct from the resolution
    % step above): ProcessBuilder's inherited environment is a snapshot
    % Java took when the JVM started, NOT the live OS environment, so
    % without this the launched program's own getenv("PATH") would see a
    % stale value even though run_command correctly found and launched it.
    env_map.put('PATH', effective_path);

    stdin_file = '';
    if ~isempty(options.Stdin)
        stdin_file = [tempname() '.stdin'];
        fid = fopen(stdin_file, 'w');
        fwrite(fid, options.Stdin);
        fclose(fid);
        pb.redirectInput(java.io.File(stdin_file));
    end
    cleanup_stdin = onCleanup(@() delete_if_exists(stdin_file)); %#ok<NASGU>

    bg_log_file = '';
    if options.Background
        % Nothing will read the child's output pipe, so redirect it to a
        % throwaway file at the OS level rather than leaving it unconsumed,
        % which would otherwise risk the child blocking once the pipe
        % buffer fills.
        bg_log_file = [tempname() '.bglog'];
        pb.redirectOutput(java.io.File(bg_log_file));
    end

    try
        proc = pb.start();
    catch ME
        % Mirrors system()'s own failure mode: a shell that cannot find/exec
        % the target program reports exit status 127 ("command not found").
        % Extract just the IOException's own message (e.g. `Cannot run
        % program "foo": error=2, No such file or directory`) rather than
        % surfacing MATLAB's full wrapped Java stack trace.
        status = 127;
        msg_lines = strsplit(ME.message, newline);
        reason = msg_lines{1};
        for i = 1:numel(msg_lines)
            if contains(msg_lines{i}, 'IOException:')
                reason = extractAfter(msg_lines{i}, 'IOException:');
                break;
            end
        end
        cmdout = sprintf('%s: %s', argv{1}, strtrim(reason));
        return;
    end

    if options.Background
        status = NaN;
        cmdout = '';
        return;
    end

    reader = java.io.BufferedReader(java.io.InputStreamReader(proc.getInputStream(), 'UTF-8'));
    lines = {};
    while true
        line = reader.readLine();
        % readLine() returns a java.lang.String object for every real line
        % (ischar() on it is false — MATLAB does NOT auto-convert Java
        % strings to char). Only Java null (EOF) converts to a numeric [],
        % so that is the actual end-of-stream test.
        if isnumeric(line) && isempty(line)
            break;
        end
        lines{end+1} = char(line); %#ok<AGROW>
        if options.Echo
            fprintf('%s\n', lines{end});
        end
    end
    reader.close();

    status = double(proc.waitFor());
    if isempty(lines)
        cmdout = '';
    else
        cmdout = [strjoin(lines, newline) newline];
    end
end

function delete_if_exists(filepath)
    if ~isempty(filepath) && isfile(filepath)
        delete(filepath);
    end
end

function resolved = resolve_executable(name, path_str)
% Resolve a bare executable name against path_str (a PATH-style,
% pathsep-separated list of directories), mirroring what execvp()/a shell
% would do. If `name` already contains a path separator (i.e. it is an
% absolute or explicit relative path, e.g. from fullfile(bin_dir,'charm')),
% it is returned unchanged -- resolution only applies to bare names like
% 'sacct' or 'git'. Falls back to returning `name` unchanged if no match
% is found in path_str, so the caller's own error handling (status 127)
% still applies exactly as if this resolution step were not there.
    resolved = name;
    if isempty(name) || any(name == '/') || (ispc && any(name == '\'))
        return;
    end
    if ispc
        candidate_names = strcat(name, {'', '.exe', '.bat', '.cmd'});
    else
        candidate_names = {name};
    end
    dirs = strsplit(path_str, char(java.io.File.pathSeparator));
    for i = 1:numel(dirs)
        if isempty(dirs{i})
            continue;
        end
        for j = 1:numel(candidate_names)
            candidate = fullfile(dirs{i}, candidate_names{j});
            if isfile(candidate)
                resolved = candidate;
                return;
            end
        end
    end
end
