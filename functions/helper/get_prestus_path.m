function prestus_path = get_PRESTUSpath()
% GET_PRESTUSPATH  Return the absolute path to the PRESTUS toolbox root
%
% Resolves the PRESTUS root directory at runtime from mfilename('fullpath'),
% so it works regardless of the current working directory or how the toolbox
% was added to the MATLAB path. Errors if 'PRESTUS' does not appear in the
% resolved path (e.g. if the toolbox folder has been renamed).
%
% Use as:
%   prestus_path = get_PRESTUSpath()
%
% Output:
%   prestus_path - char; absolute path to the PRESTUS root directory
%                  using the OS-native file separator
%
% See also: MFILENAME, LOAD_PARAMETERS, SAFE_ADDPATH

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
