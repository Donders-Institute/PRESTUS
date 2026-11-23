function folders = transducer_library_folders(parameters, t_i)
% TRANSDUCER_LIBRARY_FOLDERS  Your own transducer libraries, most specific first
%
% PRESTUS checks these folders before its shipped library
% (config/equipment, config/calibration), in this order:
%   1. transducer.library_path         this study (per transducer)
%   2. path.transducer_library         this project (e.g. in the project copy of
%                                      config_default.yaml); a folder or a list
%   3. PRESTUS_TRANSDUCER_LIBRARY      this machine or user: an environment
%                                      variable, folders separated by pathsep
% Each folder may hold equipment entries (*_transducer.yaml, *_driver.yaml),
% calibration entries and hydrophone profiles (profiles/).
%
% Use as:
%   folders = transducer_library_folders()                 % environment only
%   folders = transducer_library_folders(parameters)       % all transducers
%   folders = transducer_library_folders(parameters, t_i)  % transducer t_i
%
% Output:
%   folders - cellstr of existing folders, duplicates removed, most specific first
%
% See also: LOAD_EQUIPMENT_CONFIG, LIBRARY_ENTRY_FILE, RESOLVE_TRANSDUCER_FROM_SERIAL

arguments
    parameters (1,1) struct = struct()
    t_i {mustBeNumeric} = []
end

    if nargin < 1; parameters = struct(); end
    if nargin < 2; t_i = []; end
    folders = {};

    % 1. study
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        idx = 1:numel(parameters.transducer);
        if ~isempty(t_i); idx = t_i; end
        for k = idx
            if iscell(parameters.transducer); t = parameters.transducer{k}; else; t = parameters.transducer(k); end
            if isfield(t, 'library_path'); folders = [folders, as_list(t.library_path)]; end %#ok<AGROW>
        end
    end
    % 2. project
    if isfield(parameters, 'path') && isfield(parameters.path, 'transducer_library')
        folders = [folders, as_list(parameters.path.transducer_library)];
    end
    % 3. machine / user
    env = strtrim(getenv('PRESTUS_TRANSDUCER_LIBRARY'));
    if ~isempty(env)
        folders = [folders, strsplit(env, pathsep)];
    end

    folders = strtrim(folders);
    folders = folders(~cellfun(@isempty, folders));
    folders = folders(cellfun(@isfolder, folders));
    [~, keep] = unique(folders, 'stable');
    folders = folders(sort(keep));
end

function c = as_list(v)
    if isempty(v) || isa(v, 'yaml.Null')
        c = {};
    elseif ischar(v) || (isstring(v) && isscalar(v))
        c = {char(v)};
    elseif iscell(v) || isstring(v)
        c = cellfun(@char, cellstr(v), 'UniformOutput', false);
        c = c(:)';
    else
        c = {};
    end
end
