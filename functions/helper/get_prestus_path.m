function prestus_path = get_PRESTUSpath()
%get_PRESTUSpath Return pathname to the PRESTUS Toolbox.
%
% DESCRIPTION:
%     get_PRESTUSpath returns the full directory pathname to the root 
%     directory in the PRESTUS Toolbox using the slash direction native to
%     the users operating system. 
%
% USAGE:
%     path = get_PRESTUSpath()

% get the full pathname of the toolbox directory
full_path = mfilename('fullpath');

% Replace all separators with /
normalized_path = strrep(full_path, '\', '/');

% Extract up to PRESTUS (handles trailing / or not)
tokens = strsplit(normalized_path, '/');
prestus_idx = find(contains(tokens, 'PRESTUS'), 1);
if isempty(prestus_idx)
    error('No PRESTUS found in path');
end

% Rebuild with original separators
prestus_path = strjoin(tokens(1:prestus_idx), filesep);
prestus_path = strrep(prestus_path, '/', filesep);
