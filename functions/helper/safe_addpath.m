function safe_addpath(root)
% SAFE_ADDPATH  Add a directory tree to the MATLAB path, skipping hidden directories
%
% Equivalent to addpath(genpath(root)) but excludes any directory whose
% path contains a component starting with '.' (e.g., .claude, .git).
% Prevents stale copies in hidden directories from shadowing current files.
%
% Use as:
%   safe_addpath(root)
%
% Input:
%   root - path to the root directory to add recursively
%
% See also: ADDPATH, GENPATH, PATH_LOG_SETUP

arguments
    root (1,:) char
end

    dirs = strsplit(genpath(root), pathsep);
    dirs = dirs(~cellfun(@has_hidden_component, dirs));
    dirs = dirs(~cellfun(@isempty, dirs));
    if ~isempty(dirs)
        addpath(dirs{:});
    end
end

function tf = has_hidden_component(d)
    parts = strsplit(d, {filesep, '/'});
    tf = any(cellfun(@(c) ~isempty(c) && c(1) == '.', parts));
end
