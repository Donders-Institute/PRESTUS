function prestus_path = get_prestus_path()
% GET_PRESTUSPATH  Return the absolute path to the PRESTUS toolbox root
%
% Resolves the PRESTUS root directory at runtime from this file's location
% (functions/helper/get_prestus_path.m → three levels up), so it works
% regardless of the current working directory, how the toolbox was added to
% the MATLAB path, or what the checkout folder is called (e.g. 'PRESTUS-main'
% from a GitHub zip, or 'prestus'). The root is confirmed by the presence of
% config/config_default.yaml.
%
% If this file has been copied outside a PRESTUS tree, the last path component
% named exactly 'PRESTUS' is used instead (not a substring: a parent folder
% like 'PRESTUS_itrusst' must not match).
%
% Use as:
%   prestus_path = get_prestus_path()
%
% Output:
%   prestus_path - char; absolute path to the PRESTUS root directory
%                  using the OS-native file separator
%
% See also: MFILENAME, LOAD_PARAMETERS, SAFE_ADDPATH

full_path = mfilename('fullpath');

% functions/helper/get_prestus_path → PRESTUS root
prestus_path = fileparts(fileparts(fileparts(full_path)));
if isfile(fullfile(prestus_path, 'config', 'config_default.yaml'))
    return
end

% Fallback: a folder named exactly 'PRESTUS' in the path
normalized_path = strrep(full_path, '\', '/');
tokens = strsplit(normalized_path, '/');
prestus_idx = find(strcmp(tokens, 'PRESTUS'), 1, 'last');
if isempty(prestus_idx)
    error('prestus:rootNotFound', ...
        ['Cannot locate the PRESTUS root from %s: expected config/config_default.yaml ' ...
         'two folders above functions/helper/.'], full_path);
end
prestus_path = strjoin(tokens(1:prestus_idx), '/');
prestus_path = strrep(prestus_path, '/', filesep);
