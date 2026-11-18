function [coord_matrix, info] = neuronav_localite_pose(parameters)
% NEURONAV_LOCALITE_POSE  Resolve the Localite 4×4 instrument pose for placement.
%
% Single entry point for pipeline Localite placement. Returns the mean 4×4
% instrument matrix of one recorded series (column 4 = position, column 1 =
% beam direction; Localite planning-image mm), to be passed to
% LOCALITE_MATRIX_TO_POSITIONS.
%
% Two sources, selected by placement.localite.file:
%
%   File mode (placement.localite.file set)
%       The XML (TriggerMarkers, GUMMarkers or InstrumentMarker; detected
%       from its content) is parsed and one series is picked.
%
%   Automatic mode (placement.localite.file empty)
%       The ingested file <path.localite_post>/sub-XXX/ses-XX/
%       InstrumentMarker_sub-XXX_ses-XX.xml is read. If it does not exist and
%       path.localite_raw is set, NEURONAV_INGEST_MARKERS is run first (with
%       neuronav.coil_map if configured, else a single target for the
%       requested series). See NEURONAV_RESOLVE_PATHS for the path keys.
%
% Series selection (both modes):
%   placement.localite.target   - pick the entry whose description equals
%                                 this name (ingested target name, or the
%                                 InstrumentMarker/GUMMarkers description)
%   placement.localite.position - otherwise: the index of the detected
%                                 series (default 1; after the
%                                 expected_segment_length filter). In automatic mode
%                                 it is matched against series_index in the
%                                 JSON sidecar.
%
% Other keys: placement.localite.session (default 1; n → ses-0n),
% placement.localite.expected_segment_length (file mode, TriggerMarkers).
% The subject folder is sub-%03d of parameters.subject_id (or the string
% itself if subject_id is a string).
%
% Use as:
%   [coord_matrix, info] = neuronav_localite_pose(parameters)
%
% Output:
%   coord_matrix - [4×4] mean instrument matrix (not transposed)
%   info         - struct: .source (file read), .mode ('file'/'automatic'),
%                  .target (description, may be ''), .series_index
%
% See also: LOCALITE_MATRIX_TO_POSITIONS, NEURONAV_LOCALITE_FILE_SERIES, NEURONAV_INGEST_MARKERS,
%           NEURONAV_SELECT_LOCALITE, NEURONAV_RESOLVE_PATHS

    lc = struct();
    if isfield(parameters, 'placement') && isfield(parameters.placement, 'localite')
        lc = parameters.placement.localite;
    end
    target   = char(get_nonempty_field(lc, 'target', ''));
    position = get_nonempty_field(lc, 'position', 1);

    file = char(get_nonempty_field(lc, 'file', ''));
    if ~isempty(file)
        if ~isfile(file)
            error('PRESTUS:localite:fileNotFound', ...
                'placement.localite.file does not exist: %s', file);
        end
        [coord_matrix, info] = pose_from_file(file, target, position, ...
            get_nonempty_field(lc, 'expected_segment_length', []));
    else
        [coord_matrix, info] = pose_from_ingested(parameters, lc, target, position);
    end

    if all(abs(coord_matrix(1:3, 4)) < 1e-9)
        error('PRESTUS:localite:untracked', ...
            ['Selected Localite series (%s) is untracked: its position is at the origin. ' ...
             'Choose another placement.localite.position / target.'], describe(info));
    end
end


% ======================================================================
function [M, info] = pose_from_file(file, target, position, seg_len)
    series = neuronav_localite_file_series(file, seg_len);
    descr = {series.description};
    if ~isempty(target)
        k = find(strcmp(descr, target), 1);
        if isempty(k)
            error('PRESTUS:localite:targetNotFound', ...
                'placement.localite.target "%s" not found in %s (descriptions: %s)', ...
                target, file, strjoin(unique(descr(~cellfun(@isempty, descr))), ', '));
        end
    else
        k = position;
        if k < 1 || k > numel(series)
            error('PRESTUS:localite:noSeries', ...
                'placement.localite.position = %d, but %s has %d series.', ...
                position, file, numel(series));
        end
        if numel(series) > 1
            fprintf('[localite] %d series in %s; using series %d.\n', numel(series), file, k);
        end
    end

    if ~series(k).is_pose
        error('PRESTUS:localite:notAPose', ...
            ['Series %d of %s is a planned EntryTargetPair (a target point), not a ' ...
             'recorded instrument pose, and cannot place the transducer.'], k, file);
    end

    M = series(k).matrix;
    info = struct('source', file, 'mode', 'file', 'target', descr{k}, 'series_index', k);
end


% ======================================================================
function [M, info] = pose_from_ingested(parameters, lc, target, position)
    paths = neuronav_resolve_paths(parameters);
    if isempty(paths.post)
        error('PRESTUS:localite:noSource', ...
            ['Localite placement needs placement.localite.file (single XML) or ' ...
             'path.localite_post (ingested outputs, created from path.localite_raw ' ...
             'if missing). See doc_placement_neuronav.md, "Localite I/O".']);
    end

    if isnumeric(parameters.subject_id)
        sub_id = sprintf('sub-%03d', parameters.subject_id);
    else
        sub_id = char(parameters.subject_id);
    end
    ses = get_nonempty_field(lc, 'session', 1);
    if isnumeric(ses), ses = sprintf('ses-%02d', double(ses)); end
    ses = char(ses);

    stem      = sprintf('InstrumentMarker_%s_%s', sub_id, ses);
    xml_path  = fullfile(paths.post, sub_id, ses, [stem '.xml']);
    json_path = fullfile(paths.post, sub_id, ses, [stem '.json']);

    if ~isfile(xml_path)
        if isempty(paths.raw)
            error('PRESTUS:localite:notIngested', ...
                ['No ingested Localite markers at %s, and path.localite_raw is not set ' ...
                 'to create them. Run neuronav_ingest_markers / neuronav_gui first.'], xml_path);
        end
        p_ing = parameters;
        p_ing.path.localite_raw  = paths.raw;
        p_ing.path.localite_post = paths.post;
        coil_map = [];
        if isfield(parameters, 'neuronav') && isfield(parameters.neuronav, 'coil_map')
            coil_map = parameters.neuronav.coil_map;
        end
        if ~isempty(coil_map)
            target_map = neuronav_build_target_map(coil_map);
        else
            name = target;
            if isempty(name), name = sprintf('series_%d', position); end
            target_map = struct('series_index', position, 'transducer_id', 0, ...
                'target_name', name);
        end
        fprintf('[localite] No ingested markers for %s %s; ingesting from %s\n', ...
            sub_id, ses, paths.raw);
        xml_out = neuronav_ingest_markers(p_ing, sub_id, ses, target_map);
        if isempty(xml_out) || ~isfile(xml_path)
            error('PRESTUS:localite:noMarkers', ...
                'Ingestion found no usable Localite markers for %s %s under %s.', ...
                sub_id, ses, paths.raw);
        end
    end

    localite = readstruct(xml_path);
    if ~isfield(localite, 'InstrumentMarker')
        error('PRESTUS:localite:noMarkers', 'No InstrumentMarker entries in %s', xml_path);
    end
    marks = localite.InstrumentMarker;
    if ~iscell(marks), marks = num2cell(marks); end
    % string(): readstruct returns numeric-looking descriptions as numbers
    names = cellfun(@(m) strtrim(char(string(m.Marker.descriptionAttribute))), marks, 'UniformOutput', false);

    series_index = NaN;
    meta = struct();
    if isfile(json_path), meta = jsondecode(fileread(json_path)); end

    if ~isempty(target)
        k = find(strcmp(names, target), 1);
        if isempty(k)
            error('PRESTUS:localite:targetNotFound', ...
                'placement.localite.target "%s" not in %s (targets: %s)', ...
                target, xml_path, strjoin(names, ', '));
        end
    elseif ~isempty(fieldnames(meta))
        idx = cellfun(@(n) series_of(meta, n), names);
        k = find(idx == position, 1);
        if isempty(k)
            error('PRESTUS:localite:noSeries', ...
                ['Raw series %d is not among the ingested targets in %s (series: %s). ' ...
                 'Re-ingest it (neuronav_gui) or set placement.localite.target.'], ...
                position, xml_path, mat2str(idx(~isnan(idx))));
        end
    else
        warning('PRESTUS:localite:noSidecar', ...
            'No JSON sidecar next to %s; using entry %d in file order.', xml_path, position);
        k = position;
        if k < 1 || k > numel(marks)
            error('PRESTUS:localite:noSeries', '%s has %d entries.', xml_path, numel(marks));
        end
    end
    if ~isempty(fieldnames(meta)), series_index = series_of(meta, names{k}); end

    M = matrix4d_to_full(marks{k}.Marker.Matrix4D);
    info = struct('source', xml_path, 'mode', 'automatic', 'target', names{k}, ...
        'series_index', series_index);
    fprintf('[localite] Using target "%s" from %s\n', names{k}, xml_path);
end


% ======================================================================
function s = series_of(meta, name)
    key = matlab.lang.makeValidName(name);
    s = NaN;
    if isfield(meta, key) && isfield(meta.(key), 'series_index')
        s = meta.(key).series_index;
    end
end

function M = matrix4d_to_full(m)
    M = zeros(4);
    for row = 0:3
        for col = 0:3
            f = sprintf('data%d%dAttribute', row, col);
            if isfield(m, f), M(row+1, col+1) = double(m.(f)); end
        end
    end
end

function t = describe(info)
    t = sprintf('%s, series %g', info.source, info.series_index);
    if ~isempty(info.target), t = sprintf('%s, target "%s"', t, info.target); end
end
