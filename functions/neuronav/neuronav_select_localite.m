function localite = neuronav_select_localite(parameters, sub_id, ses_id, target_name)
% NEURONAV_SELECT_LOCALITE  Read the normalised InstrumentMarker XML for a session.
%
% Locates the InstrumentMarker XML written by NEURONAV_INGEST_MARKERS under
%   <parameters.path.localite_post>/<sub_id>/<ses_id>/
% and returns it as a parsed struct.  An optional target_name filter selects
% a single InstrumentMarker entry by its description attribute.
%
% This function is Stage 2 of the neuronav pipeline.  Raw marker selection,
% de-duplication, and series statistics are handled upstream by
% NEURONAV_INGEST_MARKERS; this function is a thin reader only.
%
% Use as:
%   localite = neuronav_select_localite(parameters, sub_id, ses_id)
%   localite = neuronav_select_localite(parameters, sub_id, ses_id, target_name)
%
% Input:
%   parameters  - PRESTUS config struct with parameters.path.localite_post field
%   sub_id      - subject identifier string (e.g. 'sub-020')
%   ses_id      - session identifier; numeric or string (e.g. 2 or 'ses-02')
%   target_name - (optional) filter to a single marker by description string
%
% Output:
%   localite    - struct parsed from the InstrumentMarker XML, or [] if not found
%
% See also: NEURONAV_INGEST_MARKERS, NEURONAV_COMPUTE_SERIES_STATISTICS

    localite = [];

    % Harmonise session string
    if isnumeric(ses_id)
        session = sprintf('ses-%02d', double(ses_id));
    else
        session = char(ses_id);
    end

    % Locate the normalised XML
    stem     = sprintf('InstrumentMarker_%s_%s', sub_id, session);
    xml_path = fullfile(parameters.path.localite_post, sub_id, session, [stem '.xml']);

    if ~isfile(xml_path)
        warning('neuronav_select_localite: XML not found: %s\n  Run neuronav_ingest_markers first.', xml_path);
        return;
    end

    try
        localite = readstruct(xml_path);
    catch ME
        warning('neuronav_select_localite: failed to read %s\n  %s', xml_path, ME.message);
        return;
    end

    % Optional: filter to a single named target
    if nargin >= 4 && ~isempty(target_name)
        localite = filter_by_target(localite, target_name, xml_path);
    end
end


function localite = filter_by_target(localite, target_name, xml_path)
% Return a localite struct containing only the InstrumentMarker whose
% Marker description matches target_name.

    if ~isfield(localite, 'InstrumentMarker')
        warning('neuronav_select_localite: no InstrumentMarker entries in %s', xml_path);
        localite = [];
        return;
    end

    markers = localite.InstrumentMarker;
    if ~iscell(markers), markers = num2cell(markers); end

    match = [];
    for i = 1:numel(markers)
        if isfield(markers{i}, 'Marker') && isfield(markers{i}.Marker, 'descriptionAttribute')
            if strcmp(strtrim(char(markers{i}.Marker.descriptionAttribute)), strtrim(target_name))
                match = markers{i};
                break;
            end
        end
    end

    if isempty(match)
        warning('neuronav_select_localite: target "%s" not found in %s', target_name, xml_path);
        localite = [];
        return;
    end

    localite.InstrumentMarker = match;
end
