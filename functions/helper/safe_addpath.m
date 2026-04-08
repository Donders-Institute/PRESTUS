function safe_addpath(root)
%SAFE_ADDPATH  Add a directory tree to the MATLAB path, skipping hidden dirs.
%
%   safe_addpath(root)
%
%   Like addpath(genpath(root)) but excludes any directory whose path
%   contains a component starting with '.' (e.g. .claude, .git).
%   Older MATLAB versions do not filter these automatically, which can
%   cause stale copies in hidden directories to shadow current files.

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
