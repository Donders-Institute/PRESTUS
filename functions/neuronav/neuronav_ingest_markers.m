function [xml_path, json_path, results] = neuronav_ingest_markers(parameters, sub_id, ses_id, target_map)
% NEURONAV_INGEST_MARKERS  Ingest raw Localite markers and write normalised outputs.
%
% Reads TriggerMarkers or GUMMarkers from a raw session folder, computes
% per-series statistics, maps series to named targets via a caller-supplied
% target_map, then writes:
%   1. An InstrumentMarker XML in the Localite RAS format (one entry per
%      target, description = target_map.target_name).
%   2. A JSON sidecar with placement QC statistics (N_pulses, position and
%      rotation deviations, source file, marker type).
%
% Both files are written to
%   <parameters.path.localite_post>/<sub_id>/<ses_id>/
% and are named
%   InstrumentMarker_<sub_id>_<ses_id>.xml / .json
%
% This function is Stage 1 of the neuronav pipeline.  Downstream functions
% (neuronav_select_localite, neuronav_convert_trigger_to_voxels, …) read the
% normalised XML; they do NOT need to know the original marker type.
%
% Use as:
%   [xml_path, json_path] = neuronav_ingest_markers(parameters, sub_id, ses_id, target_map)
%
% Required inputs:
%   parameters - PRESTUS config struct.  Relevant fields:
%                  parameters.path.localite_raw          raw input root  (<sub>/<ses>/loc/…)
%                  parameters.path.localite_post         derived output root (created if absent)
%                  parameters.placement.localite.loc_folder           (default: 'loc')
%                  parameters.placement.localite.session_folder_pattern (default: '*')
%                  parameters.io.overwrite_files                      ('always'/'never'/other)
%   sub_id     - subject string, e.g. 'sub-020'
%   ses_id     - session string or number, e.g. 'ses-02' or 2
%   target_map - struct array with one entry per target.  Required fields:
%                  .series_index   integer — which results{} cell this target maps to
%                  .transducer_id  integer — transducer index (informational; written to JSON)
%                  .target_name    char    — written to XML description and JSON key
%                Optional fields:
%                  .color          char    — hex color string (default '#ffff00')
%
% Outputs:
%   xml_path   - full path to the written InstrumentMarker XML file
%   json_path  - full path to the written JSON sidecar file
%
% See also: NEURONAV_SELECT_LOCALITE, NEURONAV_COMPUTE_SERIES_STATISTICS

    % ------------------------------------------------------------------
    % Resolve settings from parameters
    % ------------------------------------------------------------------
    loc_folder             = get_nested(parameters, {'placement','localite','loc_folder'},             'loc');
    session_folder_pattern = get_nested(parameters, {'placement','localite','session_folder_pattern'}, '*');
    voxel_size             = get_nested(parameters, {'placement','localite','voxel_size'},             0.9);
    expected_segment_length= get_nested(parameters, {'placement','localite','expected_segment_length'},80);
    overwrite_flag         = get_nested(parameters, {'io','overwrite_files'},                          'always');
    overwrite              = ~strcmpi(overwrite_flag, 'never');
    coil_map               = get_nested(parameters, {'neuronav','coil_map'},                           []);

    % ------------------------------------------------------------------
    % Harmonise session string
    % ------------------------------------------------------------------
    if isnumeric(ses_id)
        session = sprintf('ses-%02d', double(ses_id));
    else
        session = char(ses_id);
    end

    % Previous session for de-duplication
    ses_num = str2double(regexp(session, '\d+', 'match', 'once'));
    if ses_num > 1
        session_prev = sprintf('ses-%02d', ses_num - 1);
    else
        session_prev = '';
    end

    % ------------------------------------------------------------------
    % Output paths
    % ------------------------------------------------------------------
    out_dir  = fullfile(parameters.path.localite_post, sub_id, session);
    stem     = sprintf('InstrumentMarker_%s_%s', sub_id, session);
    xml_path = fullfile(out_dir, [stem '.xml']);
    json_path= fullfile(out_dir, [stem '.json']);

    results = {};

    if ~overwrite && isfile(xml_path) && isfile(json_path)
        fprintf('  [ingest] outputs exist, skipping %s %s\n', sub_id, session);
        return;
    end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    % ------------------------------------------------------------------
    % Select raw marker file(s) and compute per-series statistics
    % ------------------------------------------------------------------
    raw_root = fullfile(parameters.path.localite_raw, sub_id, session, loc_folder);

    if ~isempty(coil_map)
        % --- Multi-coil path: one file per coil in coil_map --------------
        results    = {};
        markertype = 'TriggerMarkers';
        source_file = {};
        for c = 1:numel(coil_map)
            coil_idx = coil_map(c).coil;
            [localite_c, src_c] = select_coil_file(raw_root, coil_idx, sub_id, session);
            if isempty(localite_c)
                warning('neuronav_ingest_markers: no usable file for coil %d (%s %s)', ...
                    coil_idx, sub_id, session);
                results{end+1} = []; %#ok<AGROW>
                source_file{end+1} = ''; %#ok<AGROW>
                continue;
            end
            stats_c = neuronav_compute_series_statistics(localite_c, voxel_size, ...
                expected_segment_length, markertype);
            if isempty(stats_c)
                warning('neuronav_ingest_markers: no valid series for coil %d (%s %s)', ...
                    coil_idx, sub_id, session);
                results{end+1} = []; %#ok<AGROW>
            else
                results{end+1} = stats_c{1}; %#ok<AGROW>
            end
            source_file{end+1} = src_c; %#ok<AGROW>
        end
        if all(cellfun(@isempty, results))
            warning('neuronav_ingest_markers: all coils empty for %s %s', sub_id, session);
            xml_path = ''; json_path = ''; results = {}; return;
        end
    else
        % --- Legacy single-coil path (Coil0) -----------------------------
        [localite, markertype, source_file] = select_raw_file(raw_root, sub_id, session, ...
            session_prev, parameters.path.localite_raw, loc_folder, session_folder_pattern);

        if isempty(localite)
            warning('neuronav_ingest_markers: no usable marker file for %s %s', sub_id, session);
            xml_path = ''; json_path = ''; results = {}; return;
        end

        results = neuronav_compute_series_statistics(localite, voxel_size, ...
            expected_segment_length, markertype);

        if isempty(results)
            warning('neuronav_ingest_markers: no valid series for %s %s', sub_id, session);
            xml_path = ''; json_path = ''; results = {}; return;
        end
    end

    % ------------------------------------------------------------------
    % Validate target_map against available series
    % ------------------------------------------------------------------
    n_targets = numel(target_map);
    n_series  = numel(results);
    for t = 1:n_targets
        si = target_map(t).series_index;
        if si > n_series
            error('neuronav_ingest_markers: target_map(%d).series_index=%d but only %d series found for %s %s', ...
                t, si, n_series, sub_id, session);
        end
        if isempty(results{si})
            warning('neuronav_ingest_markers: series %d is empty (coil dropout?) for %s %s — skipping target "%s"', ...
                si, sub_id, session, target_map(t).target_name);
        end
    end

    % Remove targets whose series is empty
    valid_targets = arrayfun(@(tm) ~isempty(results{tm.series_index}), target_map);
    target_map = target_map(valid_targets);
    if isempty(target_map)
        warning('neuronav_ingest_markers: no valid targets remain for %s %s', sub_id, session);
        xml_path = ''; json_path = ''; return;
    end

    % Resolve source_file to a single string for JSON sidecar
    if iscell(source_file)
        src_str = strjoin(source_file(~cellfun(@isempty, source_file)), '; ');
    else
        src_str = source_file;
    end

    % ------------------------------------------------------------------
    % Write InstrumentMarker XML
    % ------------------------------------------------------------------
    write_instrument_marker_xml(xml_path, target_map, results);

    % ------------------------------------------------------------------
    % Write JSON sidecar
    % ------------------------------------------------------------------
    write_json_sidecar(json_path, target_map, results, markertype, src_str);

    fprintf('  [ingest] %s %s → %s\n', sub_id, session, xml_path);
    fprintf('  [ingest] %s %s → %s\n', sub_id, session, json_path);
end


% ======================================================================
%  LOCAL HELPER: select raw marker file
% ======================================================================
function [localite, markertype, source_name] = select_raw_file( ...
        raw_root, sub_id, session, session_prev, data_postlocalite, loc_folder, ses_pat)

    localite    = [];
    markertype  = '';
    source_name = '';

    for attempt = 1:2
        if attempt == 1
            mt      = 'TriggerMarkers';
            pattern = 'TriggerMarkers_Coil0*.xml';
            % Use ** to find TMSTrigger/ at any depth (handles direct
            % Session_*/ layout and nested <patient>/Sessions/Session_/).
            files = dir(fullfile(raw_root, '**', 'TMSTrigger', pattern));
            if ~isempty(files)
                files = files([files.bytes] > 10000);
            end
            for i = 1:numel(files)
                tok = regexp(files(i).name, '_(\d{17})', 'tokens', 'once');
                if ~isempty(tok)
                    try, files(i).dt = datetime(tok{1}, 'InputFormat', 'yyyyMMddHHmmssSSS');
                    catch, files(i).dt = NaT; end
                else
                    files(i).dt = NaT;
                end
            end
        else
            mt      = 'GUMMarkers';
            pattern = 'GUMMarkers*.xml';
            files = dir(fullfile(raw_root, '**', 'GUMMarkers', pattern));
            for i = 1:numel(files)
                files(i).dt = datetime(files(i).datenum, 'ConvertFrom', 'datenum');
            end
        end

        % Remove invalid timestamps
        if isempty(files), continue; end
        valid = ~arrayfun(@(f) isnat(f.dt), files);
        files = files(valid);
        if isempty(files), continue; end

        % Sort newest first
        [~, sidx] = sort([files.dt], 'descend');
        files = files(sidx);

        % Load previous-session files for de-duplication
        prev_files = [];
        if ~isempty(session_prev)
            prev_root = fullfile(data_postlocalite, sub_id, session_prev, loc_folder);
            if attempt == 1
                prev_files = dir(fullfile(prev_root, '**', 'TMSTrigger', pattern));
                if ~isempty(prev_files)
                    prev_files = prev_files([prev_files.bytes] > 10000);
                end
            else
                prev_files = dir(fullfile(prev_root, '**', 'GUMMarkers', pattern));
            end
        end

        % Pick first non-duplicate
        for k = 1:numel(files)
            fpath = fullfile(files(k).folder, files(k).name);
            is_dup = false;
            for p = 1:numel(prev_files)
                ppath = fullfile(prev_files(p).folder, prev_files(p).name);
                [s, ~] = system(sprintf('diff "%s" "%s"', fpath, ppath));
                if s == 0, is_dup = true; break; end
            end
            if is_dup, continue; end

            try
                localite    = readstruct(fpath);
                markertype  = mt;
                source_name = files(k).name;
                return;
            catch
                warning('neuronav_ingest_markers: failed to read %s', fpath);
            end
        end
    end
end


% ======================================================================
%  LOCAL HELPER: write InstrumentMarker XML
% ======================================================================
function write_instrument_marker_xml(xml_path, target_map, results)

    default_colors = {'#ffff00','#ff0000','#00ff00','#0000ff','#ff00ff','#00ffff','#ffc800'};

    fid = fopen(xml_path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('neuronav_ingest_markers: cannot open %s for writing', xml_path);
    end

    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<InstrumentMarkerList coordinateSpace="RAS">\n');
    fprintf(fid, '    <!--All positions are saved in the coordinate system of the corresponding medical data.\n');
    fprintf(fid, 'NIfTI image data are recommended using RAS system (x-axis increases from the left hand side to the right hand side of the patient,\n');
    fprintf(fid, 'y-axis increases from the posterior side to the anterior side of the patient and z-axis increases from the feet toward the head of the patient).-->\n');

    for t = 1:numel(target_map)
        tm  = target_map(t);
        res = results{tm.series_index};

        % Mean 4×4 matrix: results stores [1×4×4]; reshape to [4×4]
        M = reshape(squeeze(res.matrix4d_mean), [4, 4])';

        if isfield(tm, 'color') && ~isempty(tm.color)
            color = tm.color;
        else
            color = default_colors{mod(t-1, numel(default_colors)) + 1};
        end
        if iscell(color), color = color{1}; end

        fprintf(fid, '    <InstrumentMarker alwaysVisible="false" index="%d" selected="false">\n', t-1);
        fprintf(fid, '        <Marker additionalInformation="" color="%s"\n', color);
        fprintf(fid, '            description="%s" set="true">\n', tm.target_name);
        fprintf(fid, '            <Matrix4D');
        for row = 0:3
            for col = 0:3
                fprintf(fid, ' data%d%d="%.17g"', row, col, M(row+1, col+1));
            end
        end
        fprintf(fid, '/>\n');
        fprintf(fid, '        </Marker>\n');
        fprintf(fid, '    </InstrumentMarker>\n');
    end

    fprintf(fid, '</InstrumentMarkerList>\n');
    fclose(fid);
end


% ======================================================================
%  LOCAL HELPER: write JSON sidecar
% ======================================================================
function write_json_sidecar(json_path, target_map, results, markertype, source_file)

    entries = struct();
    for t = 1:numel(target_map)
        tm  = target_map(t);
        res = results{tm.series_index};

        s.series_index   = tm.series_index;
        s.transducer_id  = tm.transducer_id;
        s.markertype     = markertype;
        s.source_file    = source_file;
        s.N_pulses       = res.N_pulses;
        s.position_dev_mm = struct( ...
            'mean', res.position_dev_mm.mean, ...
            'std',  res.position_dev_mm.std);
        s.rotation_dev_rad = struct( ...
            'mean', res.rotation_dev_rad.mean, ...
            'std',  res.rotation_dev_rad.std);

        entries.(matlab.lang.makeValidName(tm.target_name)) = s;
    end

    % jsonencode is available since R2016b
    json_str = jsonencode(entries, 'PrettyPrint', true);
    fid = fopen(json_path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('neuronav_ingest_markers: cannot open %s for writing', json_path);
    end
    fprintf(fid, '%s\n', json_str);
    fclose(fid);
end


% ======================================================================
%  LOCAL HELPER: select file for a specific coil index
% ======================================================================
function [localite, source_name] = select_coil_file(raw_root, coil_idx, sub_id, session)
% Returns the localite struct for a given coil index, checking for
% all-zero positions and falling back to other timestamps if needed.

    localite    = [];
    source_name = '';

    pattern = sprintf('TriggerMarkers_Coil%d_*.xml', coil_idx);
    files   = dir(fullfile(raw_root, '**', 'TMSTrigger', pattern));
    files   = files([files.bytes] > 10000);

    if isempty(files)
        warning('neuronav_ingest_markers: no TriggerMarkers files found for Coil%d (%s %s)', ...
            coil_idx, sub_id, session);
        return;
    end

    % Sort newest first
    for i = 1:numel(files)
        tok = regexp(files(i).name, '_(\d{17})', 'tokens', 'once');
        if ~isempty(tok)
            try, files(i).dt = datetime(tok{1}, 'InputFormat', 'yyyyMMddHHmmssSSS');
            catch, files(i).dt = NaT; end
        else
            files(i).dt = NaT;
        end
    end
    valid = ~arrayfun(@(f) isnat(f.dt), files);
    files = files(valid);
    if isempty(files), return; end
    [~, sidx] = sort([files.dt], 'descend');
    files = files(sidx);

    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);
        try
            loc = readstruct(fpath);
        catch
            warning('neuronav_ingest_markers: failed to read %s', fpath);
            continue;
        end

        % Check for all-zero positions
        if isfield(loc, 'TriggerMarker')
            markers = loc.TriggerMarker;
            xs = arrayfun(@(m) m.Matrix4D.data03Attribute, markers);
            ys = arrayfun(@(m) m.Matrix4D.data13Attribute, markers);
            zs = arrayfun(@(m) m.Matrix4D.data23Attribute, markers);
            if all(xs == 0) && all(ys == 0) && all(zs == 0)
                warning(['neuronav_ingest_markers: Coil%d positions are all-zero in %s ' ...
                    '(%s %s) — trying next timestamp'], ...
                    coil_idx, files(k).name, sub_id, session);
                continue;
            end
        end

        localite    = loc;
        source_name = files(k).name;
        return;
    end

    warning('neuronav_ingest_markers: no non-zero file found for Coil%d (%s %s)', ...
        coil_idx, sub_id, session);
end


% ======================================================================
%  LOCAL HELPER: safe nested field access with default
% ======================================================================
function val = get_nested(s, fields, default)
    val = s;
    for i = 1:numel(fields)
        if isstruct(val) && isfield(val, fields{i})
            val = val.(fields{i});
        else
            val = default;
            return;
        end
    end
end
