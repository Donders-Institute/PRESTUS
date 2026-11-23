function D = neuronav_detect_series(parameters, sub_id, ses_id)
% NEURONAV_DETECT_SERIES  Detect the recorded series of one session without writing.
%
% Runs the file selection and series statistics of NEURONAV_INGEST_MARKERS
% into a temporary folder (nothing is written to path.localite_post) and
% collects the session's named planned targets. This is the "detect" step
% shared by neuronav_gui and NEURONAV_INGEST_BATCH; name the series with
% NEURONAV_NAME_SERIES, then write with NEURONAV_INGEST_MARKERS.
%
% Use as:
%   D = neuronav_detect_series(parameters, sub_id, ses_id)
%
% Input:
%   parameters - config struct; uses path.localite_raw, placement.localite.*
%                (loc_folder, session_folder_pattern, expected_segment_length,
%                voxel_size) and neuronav.coil_map
%   sub_id     - e.g. 'sub-001'
%   ses_id     - 'ses-02' or 2
%
% Output:
%   D - struct with fields
%         .sub_id, .session     identifiers ('sub-001', 'ses-02')
%         .results              1×N cell of series stats (empty = coil dropout)
%         .markertype           'TriggerMarkers' / 'GUMMarkers' ('' if none)
%         .sources              1×N cell, raw file of each series
%         .coil_ids             1×N coil index per series ([] single-coil)
%         .coil_map             normalized coil_map (struct array: .coil,
%                               .sequences(:).target), [] single-coil
%         .typed_names          1×N names given in coil_map ('' = none)
%         .candidates           named planned targets (NEURONAV_READ_LABEL_CANDIDATES)
%         .log                  text printed by the ingest step
%
% See also: NEURONAV_NAME_SERIES, NEURONAV_INGEST_BATCH, NEURONAV_INGEST_MARKERS

arguments
    parameters (1,1) struct
    sub_id     {mustBeTextScalar}
    ses_id
end

    if isnumeric(ses_id), ses_id = sprintf('ses-%02d', double(ses_id)); end
    ses_id = char(ses_id);
    sub_id = char(sub_id);

    paths = neuronav_resolve_paths(parameters);
    if isempty(paths.raw)
        error('neuronav_detect_series:paths', 'Set path.localite_raw (raw Localite exports).');
    end

    coil_map = normalize_coil_map(get_nested(parameters, {'neuronav', 'coil_map'}, []));
    loc = get_nested(parameters, {'placement', 'localite', 'loc_folder'}, 'loc');

    % Probe target_map: series 1 (single-coil) or one entry per expected
    % series (multi-coil), so the JSON sidecar records type and source files.
    series_coils = [];
    typed = {};
    if ~isempty(coil_map)
        tm = neuronav_build_target_map(coil_map);
        series_coils = [tm.transducer_id];
        typed = {tm.target_name};
        typed(~cellfun(@isempty, regexp(typed, '^coil\d+(_\d+)?$', 'once'))) = {''};  % defaults
        probe = struct('series_index', num2cell(1:numel(tm)), ...
            'transducer_id', num2cell(series_coils), ...
            'target_name', arrayfun(@(k) sprintf('probe_%d', k), 1:numel(tm), 'UniformOutput', false));
    else
        probe = struct('series_index', 1, 'transducer_id', 0, 'target_name', 'probe_1');
    end

    p = parameters;
    p.path.localite_raw  = paths.raw;
    p.path.localite_post = tempname;
    p.io.overwrite_files = 'always';
    if ~isempty(coil_map), p.neuronav.coil_map = coil_map; end
    cleanup = onCleanup(@() remove_dir(p.path.localite_post));

    out = cell(1, 3);
    bt = warning('off', 'backtrace');
    restore = onCleanup(@() warning(bt));
    log = evalc('[out{:}] = neuronav_ingest_markers(p, sub_id, ses_id, probe);');
    log = regexprep(strrep(log, p.path.localite_post, '<temp>'), '</?a[^>]*>', '');
    [~, json_path, results] = out{:};

    D = struct('sub_id', sub_id, 'session', ses_id, 'results', {results}, ...
        'markertype', '', 'sources', {repmat({''}, 1, numel(results))}, ...
        'coil_ids', series_coils, 'coil_map', coil_map, ...
        'typed_names', {typed}, 'candidates', [], 'log', log);
    if isempty(coil_map), D.typed_names = repmat({''}, 1, numel(results)); end

    if ~isempty(json_path) && isfile(json_path)
        meta = jsondecode(fileread(json_path));
        fn = fieldnames(meta);
        D.markertype = meta.(fn{1}).markertype;
        if isempty(coil_map)
            D.sources(:) = {meta.(fn{1}).source_file};
        else
            for f = 1:numel(fn)
                D.sources{meta.(fn{f}).series_index} = meta.(fn{f}).source_file;
            end
        end
    end

    D.candidates = neuronav_read_label_candidates(paths.raw, sub_id, ses_id, loc);
end


% ======================================================================
function cm = normalize_coil_map(cm)
% Struct array with .coil and .sequences(:).target, from a struct array or
% a YAML-loaded cell (entries with either .target or .sequences).
    if isempty(cm), cm = []; return; end
    if isstruct(cm), cm = num2cell(cm); end
    out = struct('coil', {}, 'sequences', {});
    for i = 1:numel(cm)
        e = cm{i};
        if isfield(e, 'sequences') && ~isempty(e.sequences)
            seqs = e.sequences;
            if iscell(seqs), seqs = [seqs{:}]; end
            names = arrayfun(@(q) char(string(q.target)), seqs, 'UniformOutput', false);
        elseif isfield(e, 'target') && ~isempty(e.target)
            names = cellstr(string(e.target));
        else
            names = {sprintf('coil%d', e.coil)};
        end
        out(i).coil = double(e.coil);
        out(i).sequences = struct('target', names);
    end
    cm = out;
end

function v = get_nested(s, fields, default)
    v = s;
    for i = 1:numel(fields)
        if isstruct(v) && isscalar(v) && isfield(v, fields{i})
            v = v.(fields{i});
        else
            v = default; return;
        end
    end
end

function remove_dir(d)
    if isfolder(d), try, rmdir(d, 's'); catch, end, end
end
