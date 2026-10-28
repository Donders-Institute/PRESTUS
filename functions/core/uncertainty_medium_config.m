function file = uncertainty_medium_config(parameters, variant)
% UNCERTAINTY_MEDIUM_CONFIG  Medium-property YAML of one uncertainty variant.
%
%   file = uncertainty_medium_config(parameters, 'liberal')
%   file = uncertainty_medium_config(parameters, 'conservative')
%
% Reads simulation.uncertainty_<variant>. An absolute path is used as is; a
% relative path or file name is looked up in config/uncertainty/. Empty or
% absent selects the shipped config/uncertainty/config_medium_<variant>.yaml.
% The file is not checked here; uncertainty_pipeline errors when it is missing.
%
% See also: UNCERTAINTY_PIPELINE

    if ~any(strcmp(variant, {'liberal', 'conservative'}))
        error('uncertainty_medium_config:badVariant', ...
            'Variant must be ''liberal'' or ''conservative'', not ''%s''.', variant);
    end
    folder = fullfile(get_prestus_path(), 'config', 'uncertainty');
    key = ['uncertainty_' variant];
    file = '';
    if isfield(parameters, 'simulation') && isfield(parameters.simulation, key)
        file = char(parameters.simulation.(key));
    end
    if isempty(file)
        file = fullfile(folder, ['config_medium_' variant '.yaml']);
    elseif ~is_absolute(file)
        file = fullfile(folder, file);
    end
end

function tf = is_absolute(p)
    tf = startsWith(p, '/') || startsWith(p, '\') || ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end
