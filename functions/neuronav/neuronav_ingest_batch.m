function summary = neuronav_ingest_batch(parameters, selection, varargin)
% NEURONAV_INGEST_BATCH  Ingest and export Localite targets for many subjects/sessions.
%
% For each selected session: detect the series (NEURONAV_DETECT_SERIES), name
% them (NEURONAV_NAME_SERIES), write the ingested XML/JSON pair
% (NEURONAV_INGEST_MARKERS) and collect one summary row per target with its
% coordinates. neuronav_gui's Batch window calls this same function, so a
% script and the GUI given the same parameters produce the same outputs.
%
% Use as:
%   summary = neuronav_ingest_batch(parameters)                 % all sessions
%   summary = neuronav_ingest_batch(parameters, selection)
%   summary = neuronav_ingest_batch(..., 'summary_file', f, 'progress', fcn)
%
% Input:
%   parameters - config struct (e.g. load_parameters('study.yaml')); uses
%                path.localite_raw / localite_post / localite_t1,
%                placement.localite.* (ingestion), neuronav.* (coil_map,
%                naming, batch_write) and, for voxel/focus coordinates,
%                transducer(1) geometry + focal distance and placement.localite
%                offsets as in LOCALITE_MATRIX_TO_POSITIONS
%   selection  - table with columns sub, ses (e.g. rows of
%                NEURONAV_SCAN_SESSIONS), or an N×2 cell {sub, ses};
%                [] / omitted = every session found
%
% Name-value options:
%   'summary_file' - CSV to merge the rows into (rows of the processed
%                    sessions are replaced); default
%                    <path.localite_post>/neuronav_ingest_summary.csv, '' = none
%   'progress'     - @(k, n, sub, ses, status) called after each session
%
% Which series are written (neuronav.batch_write):
%   'all'   (default) every tracked series (named by the naming rule)
%   'named' only series named by coil_map or the heuristic
% Untracked series (all positions at 0) and coil dropouts are never written.
%
% Output (also written to summary_file), one row per written target, plus one
% row per session without targets (status explains why):
%   sub, ses, target, series_index, coil, rule, label, label_source,
%   axis_dist_mm, depth_mm, n_pulses, pos_dev_mm, rot_dev_rad, markertype,
%   source_file,
%   pose_x/y/z        mean instrument position (Localite planning-image mm)
%   beam_x/y/z        mean beam direction (unit vector, matrix column 1)
%   trans_ras_x/y/z, focus_ras_x/y/z   transducer (bowl rear centre) / focus [mm]
%   trans_vox_x/y/z, focus_vox_x/y/z   the same in planning-T1 voxels
%                     (NaN when the transducer config or planning T1 is missing)
%   xml, status
%
% See also: NEURONAV_SCAN_SESSIONS, NEURONAV_DETECT_SERIES, NEURONAV_NAME_SERIES

    ip = inputParser;
    ip.addParameter('summary_file', []);
    ip.addParameter('progress', []);
    ip.parse(varargin{:});
    opt = ip.Results;

    paths = neuronav_resolve_paths(parameters);
    if isempty(paths.raw) || isempty(paths.post)
        error('neuronav_ingest_batch:paths', ...
            'Set path.localite_raw (raw exports) and path.localite_post (outputs).');
    end
    parameters.path.localite_raw  = paths.raw;
    parameters.path.localite_post = paths.post;
    if isempty(opt.summary_file) && ~ischar(opt.summary_file) && ~isstring(opt.summary_file)
        opt.summary_file = fullfile(paths.post, 'neuronav_ingest_summary.csv');
    end

    if nargin < 2 || isempty(selection)
        selection = neuronav_scan_sessions(parameters);
    end
    if istable(selection)
        pairs = [cellstr(selection.sub), cellstr(selection.ses)];
    else
        pairs = cellstr(selection);
    end

    write_mode = 'all';
    if isfield(parameters, 'neuronav') && isfield(parameters.neuronav, 'batch_write') && ...
            ~isempty(parameters.neuronav.batch_write)
        write_mode = lower(char(parameters.neuronav.batch_write));
    end

    rows = {};
    for s = 1:size(pairs, 1)
        sub = pairs{s, 1}; ses = pairs{s, 2};
        try
            [new_rows, status] = process_session(parameters, sub, ses, write_mode);
        catch ME
            new_rows = {}; status = ['error: ' ME.message];
        end
        if isempty(new_rows)
            new_rows = empty_row(sub, ses, status);
        end
        rows = [rows; new_rows]; %#ok<AGROW>
        fprintf('[batch] %s %s: %s\n', sub, ses, status);
        if ~isempty(opt.progress), opt.progress(s, size(pairs, 1), sub, ses, status); end
    end

    summary = cell2table(rows, 'VariableNames', summary_vars());
    if ~isempty(opt.summary_file)
        [summary_all, ok] = merge_summary(opt.summary_file, summary);
        out_file = opt.summary_file;
        if ~ok
            % Never overwrite a summary that could not be merged
            [d, f] = fileparts(out_file);
            out_file = fullfile(d, sprintf('%s_%s.csv', f, datestr(now, 'yyyymmddTHHMMSS'))); %#ok<TNOW1,DATST>
            warning('neuronav_ingest_batch:merge', ...
                'Could not merge into %s (different columns?); wrote %s instead.', ...
                opt.summary_file, out_file);
        end
        if ~isfolder(fileparts(out_file)), mkdir(fileparts(out_file)); end
        writetable(summary_all, out_file);
        fprintf('[batch] summary: %s\n', out_file);
    end
end


% ======================================================================
function [rows, status] = process_session(parameters, sub, ses, write_mode)
    rows = {};
    D = neuronav_detect_series(parameters, sub, ses);
    if isempty(D.results) || all(cellfun(@isempty, D.results))
        status = 'no usable markers';
        return;
    end
    N = neuronav_name_series(D, parameters);

    keep = [N.tracked];
    if strcmp(write_mode, 'named'), keep = keep & ~strcmp({N.rule}, 'default'); end
    if ~any(keep)
        status = sprintf('no series to write (%d detected, batch_write = %s)', ...
            numel(N), write_mode);
        return;
    end

    idx = find(keep);
    coils = zeros(1, numel(D.results));
    if ~isempty(D.coil_ids), coils = D.coil_ids; end
    target_map = struct('series_index', num2cell(idx), 'transducer_id', num2cell(coils(idx)), ...
        'target_name', {N(idx).name}, 'color', '');

    p = parameters;
    if ~isempty(D.coil_map), p.neuronav.coil_map = D.coil_map; end
    evalc('[xml_path, ~, res] = neuronav_ingest_markers(p, sub, ses, target_map);');
    if isempty(xml_path)
        status = 'ingestion wrote nothing';
        return;
    end
    status = sprintf('%d target(s) written', numel(idx));
    if isempty(res), status = 'kept existing outputs (io.overwrite_files = never)'; end

    [hdr, t1_note] = planning_header(sub, p);
    for k = idx
        r = D.results{k};
        M = reshape(squeeze(r.matrix4d_mean), [4, 4]);
        beam = M(1:3, 1)' / max(norm(M(1:3, 1)), eps);
        tr = nan(1, 3); fr = tr; tv = tr; fv = tr; note = t1_note;
        if ~isempty(hdr)
            try
                [tv, fv, tr, fr] = localite_matrix_to_positions(M, hdr, p);
            catch ME
                note = ME.message;
            end
        end
        st = status;
        if ~isempty(note), st = sprintf('%s; coordinates: %s', status, note); end
        rows(end+1, :) = {sub, ses, N(k).name, k, coils(k), N(k).rule, N(k).label, ...
            N(k).source, N(k).dist_mm, N(k).depth_mm, r.N_pulses, ...
            r.position_dev_mm.mean, r.rotation_dev_rad.mean, D.markertype, D.sources{k}, ...
            M(1,4), M(2,4), M(3,4), beam(1), beam(2), beam(3), ...
            tr(1), tr(2), tr(3), fr(1), fr(2), fr(3), ...
            tv(1), tv(2), tv(3), fv(1), fv(2), fv(3), xml_path, st}; %#ok<AGROW>
    end
end

function [hdr, note] = planning_header(sub, p)
    hdr = []; note = '';
    try
        hdr = niftiinfo(neuronav_find_planning_t1(sub, p));
    catch ME
        note = ME.message;
    end
end

function v = summary_vars()
    v = {'sub', 'ses', 'target', 'series_index', 'coil', 'rule', 'label', 'label_source', ...
        'axis_dist_mm', 'depth_mm', 'n_pulses', 'pos_dev_mm', 'rot_dev_rad', 'markertype', ...
        'source_file', 'pose_x', 'pose_y', 'pose_z', 'beam_x', 'beam_y', 'beam_z', ...
        'trans_ras_x', 'trans_ras_y', 'trans_ras_z', 'focus_ras_x', 'focus_ras_y', 'focus_ras_z', ...
        'trans_vox_x', 'trans_vox_y', 'trans_vox_z', 'focus_vox_x', 'focus_vox_y', 'focus_vox_z', ...
        'xml', 'status'};
end

function row = empty_row(sub, ses, status)
% A session without written targets: text columns '', numbers NaN.
    v = summary_vars();
    row = repmat({NaN}, 1, numel(v));
    row(ismember(v, text_vars())) = {''};
    row([1 2 end]) = {sub, ses, status};
end

function v = text_vars()
    v = {'sub', 'ses', 'target', 'rule', 'label', 'label_source', 'markertype', ...
        'source_file', 'xml', 'status'};
end

function [T, ok] = merge_summary(file, new)
% Replace the rows of the processed sessions in an existing summary.
    T = new; ok = true;
    if ~isfile(file), return; end
    try
        opts = detectImportOptions(file);
        opts = setvartype(opts, intersect(text_vars(), opts.VariableNames), 'char');
        old = readtable(file, opts);
        if ~isequal(sort(old.Properties.VariableNames), sort(new.Properties.VariableNames))
            ok = false; return;
        end
        old = old(:, new.Properties.VariableNames);
        keep = ~ismember(strcat(old.sub, '|', old.ses), strcat(new.sub, '|', new.ses));
        T = sortrows([old(keep, :); new], {'sub', 'ses'});
    catch
        T = new; ok = false;
    end
end
