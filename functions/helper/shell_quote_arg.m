function quoted = shell_quote_arg(arg)
% SHELL_QUOTE_ARG  Quote a string as a single literal argument for a shell command line
%
% Wraps arg in single quotes (embedded quotes as '\''), which a POSIX shell
% never expands; double quotes are not enough, since $(...), backticks, and
% $VAR expand inside them. On Windows it strips cmd.exe metacharacters
% (" ^ & | < > %, newlines) and double-quotes the result, which is not
% injection-proof: prefer RUN_COMMAND, which uses no shell.
%
% Use as:
%   [status, out] = run_shell_script(sprintf('cp -r %s %s', ...
%       shell_quote_arg(src), shell_quote_arg(dst)));
%
% Input:
%   arg - (1,:) char or string scalar to be used as a single shell argument
%
% Output:
%   quoted - char row vector: arg safely quoted for splicing into a shell
%            command string
%
% See also: RUN_COMMAND, RUN_SHELL_SCRIPT

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
