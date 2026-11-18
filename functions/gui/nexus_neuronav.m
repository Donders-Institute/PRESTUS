function out = nexus_neuronav(args)
% NEXUS_NEURONAV  Localite marker ingestion for the PRESTUS Nexus (through NEXUS_BRIDGE)
%
% The Localite page of the Nexus (Placement tab) calls this with one action
% per step. Settings come from the saved study config, read with
% LOAD_PARAMETERS as a run reads them (path.localite_raw / _post,
% placement.localite.*, io.overwrite_files, neuronav.*), so the page has no
% settings of its own. The steps are the ones of NEURONAV_GUI:
% NEURONAV_DETECT_SERIES + NEURONAV_NAME_SERIES, NEURONAV_CHECK_TARGET_MAP +
% NEURONAV_INGEST_MARKERS, NEURONAV_SCAN_SESSIONS and NEURONAV_SUMMARY_TO_MNI.
% The batch itself (NEURONAV_INGEST_BATCH) runs as a script, not here.
%
% Use as:
%   out = nexus_neuronav(args)
%
% Input:
%   args - struct with .action and .config (path of the saved study config), and:
%     'sessions'  -                     subjects and sessions under the raw root
%     'detect'    .sub, .ses            series, QC, names, pulse positions, raw files
%     'write'     .sub, .ses, .target_map (list of .series_index, .transducer_id,
%                                       .target_name, .color)
%     'scan'      -                     session and label tables (batch)
%     'mni'       .csv, .seg, .template, .head   summary rows in MNI space
%                 (and the head outline when .head is true)
%     'file_series' .file               series of one Localite XML export and the
%                                       transducer/focus voxels each places
%                                       (Import trajectory, placement.localite.file)
%
% Output:
%   out - struct of plain data per action (see the local functions)
%
% See also: NEXUS_BRIDGE, NEURONAV_GUI, NEURONAV_DETECT_SERIES, NEURONAV_INGEST_MARKERS,
%           NEURONAV_LOCALITE_FILE_SERIES

    p = load_parameters(char(args.config));
    switch char(args.action)
        case 'sessions'
            out = sessions(p);
        case 'detect'
            out = detect(p, char(args.sub), char(args.ses));
        case 'write'
            out = write(p, char(args.sub), char(args.ses), get_nonempty_field(args, 'target_map', []));
        case 'scan'
            out = scan(p);
        case 'mni'
            out = mni(p, args);
        case 'file_series'
            out = file_series(p, char(args.file));
        otherwise
            error('nexus_neuronav:action', 'Unknown action ''%s''.', char(args.action));
    end
end

function out = file_series(p, file)
% The series of one Localite XML, each with the positions a run places with it:
% NEURONAV_LOCALITE_FILE_SERIES, then LOCALITE_MATRIX_TO_POSITIONS with the T1
% header PREPROC_TRANSDUCER_PLACEMENT uses (1-based PRESTUS voxels).
    seg_len = get_nested(p, {'placement', 'localite', 'expected_segment_length'}, []);
    [series, markertype] = neuronav_localite_file_series(file, seg_len);
    t1_header = [];
    note = '';
    t1_file = fullfile(p.path.anat, sprintf(p.path.t1_pattern, p.subject_id));
    if isfile(t1_file)
        t1_header = niftiinfo(t1_file);
    else
        note = sprintf('No T1 at %s: positions are resolved when the run starts.', t1_file);
    end
    rows = cell(1, numel(series));
    for k = 1:numel(series)
        s = series(k);
        r = struct('index', s.index, 'description', s.description, 'n_markers', s.n_markers, ...
            'pos_dev_mm', s.pos_dev_mm, 'is_pose', s.is_pose, 'tracked', s.tracked, ...
            'trans_pos', [], 'focus_pos', [], 'error', '');
        if s.is_pose && s.tracked && ~isempty(t1_header)
            try
                [r.trans_pos, r.focus_pos] = localite_matrix_to_positions(s.matrix, t1_header, p);
            catch err
                r.error = err.message;
            end
        end
        rows{k} = r;
    end
    out = struct('file', file, 'markertype', markertype, 'series', {rows}, 'note', note);
end

function out = sessions(p)
    paths = neuronav_resolve_paths(p);
    subs = list_dirs(paths.raw, 'sub-*');
    out = struct('raw', paths.raw, 'post', paths.post, 'subjects', {cell(1, numel(subs))});
    for k = 1:numel(subs)
        out.subjects{k} = struct('sub', subs{k}, ...
            'sessions', {list_dirs(fullfile(paths.raw, subs{k}), 'ses-*')});
    end
end

function out = detect(p, sub, ses)
    D = neuronav_detect_series(p, sub, ses);
    N = neuronav_name_series(D, p);
    suggest = logical(get_nested(p, {'neuronav', 'suggest_names'}, true));
    n = numel(D.results);
    series = cell(1, n);
    for k = 1:n
        r = D.results{k};
        if isempty(D.coil_ids)
            coil = ternary(strcmp(D.markertype, 'TriggerMarkers'), '0', '—');
        else
            coil = sprintf('%d', D.coil_ids(k));
        end
        s = struct('index', k, 'coil', coil, 'transducer_id', 0, 'empty', isempty(r), ...
            'name', N(k).name, 'rule', N(k).rule, 'label', N(k).label, 'tracked', N(k).tracked, ...
            'name_from', describe_name(N(k), suggest, ~isempty(D.candidates)), ...
            'dist_mm', N(k).dist_mm, 'depth_mm', N(k).depth_mm, 'n_pulses', 0, ...
            'pos_dev_mean', NaN, 'pos_dev_std', NaN, 'rot_dev_mean', NaN, 'rot_dev_std', NaN, ...
            'flag', 'no data', 'positions', {{}}, 'mean_pos', [], 'beam_dir', [], 'pos_dev', []);
        if ~isempty(D.coil_ids); s.transducer_id = D.coil_ids(k); end
        if ~isempty(r)
            P = reshape(r.matrix4d_all(:, 1:3, 4), [], 3);
            M = reshape(r.matrix4d_mean, 4, 4);
            s.n_pulses = r.N_pulses;
            s.pos_dev_mean = r.position_dev_mm.mean; s.pos_dev_std = r.position_dev_mm.std;
            s.rot_dev_mean = r.rotation_dev_rad.mean; s.rot_dev_std = r.rotation_dev_rad.std;
            s.flag = neuronav_series_qc_flag(r);
            s.positions = num2cell(P, 2)';          % one [x y z] per pulse
            s.mean_pos = M(1:3, 4)';
            s.beam_dir = M(1:3, 1)' / max(norm(M(1:3, 1)), eps);
            s.pos_dev = r.position_dev_mm.all(:)';
        end
        series{k} = s;
    end
    cand = cell(1, numel(D.candidates));
    for k = 1:numel(D.candidates)
        c = D.candidates(k);
        cand{k} = struct('description', c.description, 'xyz', double(c.xyz(:)'), 'source', c.source);
    end
    sources = D.sources(~cellfun(@isempty, D.sources));
    out = struct('sub', D.sub_id, 'ses', D.session, 'markertype', D.markertype, ...
        'sources', {unique(sources, 'stable')}, 'series', {series}, 'candidates', {cand}, ...
        'raw_files', {raw_files(p, sub, D.session, sources, D.candidates)}, 'log', D.log, ...
        'n_from_coil_map', sum(strcmp({N.rule}, 'coil_map')), ...
        'n_from_targets', sum(strcmp({N.rule}, 'heuristic')));
end

function out = write(p, sub, ses, tm)
    paths = neuronav_resolve_paths(p);
    assert(~isempty(paths.post), 'nexus_neuronav:post', ...
        'Set path.localite_post (where ingested markers are written).');
    target_map = struct('series_index', {}, 'transducer_id', {}, 'target_name', {}, 'color', {});
    if iscell(tm); tm = [tm{:}]; end
    for k = 1:numel(tm)
        target_map(k) = struct('series_index', double(tm(k).series_index), ...
            'transducer_id', double(tm(k).transducer_id), 'target_name', strtrim(char(tm(k).target_name)), ...
            'color', strtrim(char(get_nonempty_field(tm(k), 'color', ''))));
    end
    msg = neuronav_check_target_map(target_map);
    if ~isempty(msg); error('nexus_neuronav:targetMap', '%s', msg); end
    [xml_path, json_path, res] = neuronav_ingest_markers(p, sub, ses, target_map);
    assert(~isempty(xml_path), 'nexus_neuronav:nothing', 'Ingestion wrote nothing; see the log.');
    out = struct('xml', xml_path, 'json', json_path, 'kept_existing', isempty(res), ...
        'n_targets', numel(target_map), 'json_text', fileread(json_path));
end

function out = scan(p)
    [S, L] = neuronav_scan_sessions(p);
    out = struct('sessions', {table_rows(S)}, 'labels', {table_rows(L)});
end

function out = mni(p, args)
    paths = neuronav_resolve_paths(p);
    csv = char(get_nonempty_field(args, 'csv', fullfile(paths.post, 'neuronav_ingest_summary.csv')));
    assert(isfile(csv), 'nexus_neuronav:noSummary', 'No summary CSV: %s (run the batch first).', csv);
    seg = char(get_nonempty_field(args, 'seg', char(get_nested(p, {'path', 'seg'}, ''))));
    assert(isfolder(seg), 'nexus_neuronav:seg', 'Set the SimNIBS root (path.seg): %s', seg);
    p.path.seg = seg;
    T = neuronav_summary_to_mni(csv, p, 'progress', ...
        @(k, n, sub, ~) nexus_progress(sprintf('%d / %d subjects (%s)', k, n, sub)));
    T = T(~cellfun(@isempty, T.target), :);   % sessions without written targets
    out = struct('csv', csv, 'rows', {table_rows(T)}, 'template', '', 'faces', {{}}, 'vertices', {{}});
    if logical(get_nonempty_field(args, 'head', false))
        f = char(get_nonempty_field(args, 'template', neuronav_find_mni_template(p)));
        if isempty(f)
            ok = find(cellfun(@isempty, T.mni_status), 1);
            if ~isempty(ok)
                f = fullfile(seg, sprintf('m2m_%s', T.sub{ok}), 'toMNI', 'final_tissues_MNI.nii.gz');
            end
        end
        if ~isempty(f) && isfile(f)
            [F, V] = neuronav_head_surface(f);
            out.template = f;
            out.faces = num2cell(F, 2)';
            out.vertices = num2cell(V, 2)';
        end
    end
end

%% Helpers ------------------------------------------------------------------
function rows = raw_files(p, sub, ses, sources, candidates)
% Marker files of the session and which ones the detection used (as NEURONAV_GUI)
    paths = neuronav_resolve_paths(p);
    loc = char(get_nested(p, {'placement', 'localite', 'loc_folder'}, 'loc'));
    base = fullfile(paths.raw, sub, ses, loc);
    rows = {};
    if ~isfolder(base); return; end
    trig = dir(fullfile(base, '**', 'TMSTrigger', 'TriggerMarkers_Coil*.xml'));
    gum  = dir(fullfile(base, '**', 'GUMMarkers*.xml'));
    et   = dir(fullfile(base, '**', 'EntryTarget*.xml'));
    files = [trig; gum; et];
    types = [repmat({'TriggerMarkers'}, numel(trig), 1); ...
             repmat({'GUMMarkers'}, numel(gum), 1); repmat({'EntryTarget'}, numel(et), 1)];
    label_src = {};
    if ~isempty(candidates); label_src = unique({candidates.source}); end
    for k = 1:numel(files)
        f = files(k);
        rel = strrep(f.folder, [base filesep], '');
        if strcmp(f.folder, base); rel = '.'; end
        status = '';
        if any(strcmp(sources, f.name))
            status = 'used';
        elseif strcmp(types{k}, 'TriggerMarkers') && f.bytes <= 10000
            status = 'ignored (≤ 10 kB)';
        end
        if any(strcmp(label_src, [types{k} ' ' f.name]))
            status = strtrim([status '  target labels']);
        end
        rows{end + 1} = struct('type', types{k}, 'file', f.name, 'folder', rel, ...
            'kb', round(f.bytes / 1024, 1), 'status', status); %#ok<AGROW>
    end
end

function t = describe_name(g, suggest, has_candidates)
% Where a series name came from, in words (as the "Name from" column of NEURONAV_GUI)
    lab = g.label;
    stem = regexprep(g.name, '_(rep\d+|target)$', '');
    if ~isempty(lab) && ~strcmp(lab, stem)
        lab = sprintf('"%s" → %s', lab, stem);
    elseif ~isempty(lab)
        lab = sprintf('"%s"', lab);
    end
    switch g.rule
        case 'coil_map'
            t = 'coil table';
            if ~isempty(g.label); t = sprintf('coil table (nearest planned: %s)', g.label); end
        case 'heuristic'
            if isnan(g.dist_mm)
                t = sprintf('GUMMarkers pair %s', lab);
            else
                t = sprintf('%s %s · %.1f mm off axis, %.0f mm deep', strtok(g.source), lab, g.dist_mm, g.depth_mm);
            end
        otherwise
            if ~g.tracked
                t = 'untracked';
            elseif ~suggest
                t = 'default';
            elseif ~has_candidates
                t = 'default (no planned targets)';
            elseif ~isnan(g.dist_mm)
                t = sprintf('default (nearest %.1f mm off axis)', g.dist_mm);
            else
                t = 'default (no match)';
            end
    end
end

function rows = table_rows(T)
% Table -> cell of row structs (cellstr cells as char, logicals kept)
    rows = cell(1, height(T));
    names = T.Properties.VariableNames;
    for i = 1:height(T)
        r = struct();
        for j = 1:numel(names)
            v = T.(names{j})(i, :);
            if iscell(v); v = v{1}; end
            if isstring(v); v = char(v); end
            r.(names{j}) = v;
        end
        rows{i} = r;
    end
end

function names = list_dirs(parent, pattern)
    names = {};
    if isempty(parent) || ~isfolder(parent); return; end
    d = dir(fullfile(parent, pattern));
    names = sort({d([d.isdir]).name});
end

function val = get_nested(s, fields, default)
    val = s;
    for i = 1:numel(fields)
        if isstruct(val) && isscalar(val) && isfield(val, fields{i})
            val = val.(fields{i});
        else
            val = default; return;
        end
    end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end
