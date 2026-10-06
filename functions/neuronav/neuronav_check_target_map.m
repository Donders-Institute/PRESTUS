function msg = neuronav_check_target_map(target_map)
% NEURONAV_CHECK_TARGET_MAP  Check a target_map before NEURONAV_INGEST_MARKERS writes it
%
% The checks neuronav_gui and the Localite page of the PRESTUS Nexus run on
% the series the user ticked and named: at least one series, every series
% named, names unique also after conversion to JSON keys
% (matlab.lang.makeValidName), colours empty or '#rrggbb'.
%
% Use as:
%   msg = neuronav_check_target_map(target_map)
%
% Input:
%   target_map - struct array with .series_index, .transducer_id,
%                .target_name and optional .color ([] / empty = none ticked)
%
% Output:
%   msg - '' when the map can be written, otherwise the first problem
%
% See also: NEURONAV_INGEST_MARKERS, NEURONAV_BUILD_TARGET_MAP, NEURONAV_NAME_SERIES

    msg = '';
    if isempty(target_map)
        msg = 'Tick at least one series to write.'; return;
    end
    names = arrayfun(@(e) strtrim(char(e.target_name)), target_map, 'UniformOutput', false);
    if any(cellfun(@isempty, names))
        msg = 'Every ticked series needs a target name.'; return;
    end
    keys = matlab.lang.makeValidName(names);
    if numel(unique(keys)) < numel(keys)
        msg = 'Target names must be unique (also after conversion to JSON keys).'; return;
    end
    if ~isfield(target_map, 'color'); return; end
    for k = 1:numel(target_map)
        col = strtrim(char(target_map(k).color));
        if ~isempty(col) && isempty(regexp(col, '^#[0-9a-fA-F]{6}$', 'once'))
            msg = sprintf('Colour "%s" is not a #rrggbb hex string.', col); return;
        end
    end
end
