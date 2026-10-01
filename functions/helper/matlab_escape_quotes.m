function s = matlab_escape_quotes(s)
% MATLAB_ESCAPE_QUOTES  Escape text for a single-quoted MATLAB char literal
%
% Doubles every single quote, so a path written into a generated .m script
% as '...' stays one literal (a path with an apostrophe would otherwise end
% the string and let the rest run as code).
%
% Use as:
%   fprintf(fid, 'load(''%s'');\n', matlab_escape_quotes(path))
%
% Input:
%   s - char or string scalar
%
% Output:
%   s - char with ' replaced by ''
%
% See also: SHELL_QUOTE_ARG, R_QUOTE_STRING

    s = strrep(char(s), '''', '''''');
end
