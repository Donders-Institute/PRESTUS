function out = nexus_calibration(args)
% NEXUS_CALIBRATION  Transducer calibration for the PRESTUS Nexus (through NEXUS_BRIDGE)
%
% The Calibration page of the Nexus (Transducer tab) calls this with one
% action per step; it keeps the state of each page (equipment, profiles,
% library entry, unsaved fit) between calls, keyed by args.session, and
% returns plain data that JSON can carry. All computation is the one
% CALIBRATION_GUI uses: CALIBRATION_DESIRED_PROFILES, CALIBRATION_EVALUATE,
% CALIBRATION_FIT, CALIBRATION_SAVE_ENTRY, LIBRARY_SETTINGS_TABLE.
%
% Use as:
%   out = nexus_calibration(args)
%
% Input:
%   args - struct with .action and .session, and per action:
%     'open'    .serial, .driving_system ('' = none), .library_path (user
%               folder, '' = none), optional .equip_path, .profile_root, .csv
%               ('' = the combo's char_data_path under profile_root)
%     'load'    .csv
%     'view'    .predicted, .depths (predicted), .settings (see below)
%     'fit'     .settings, .free_params, .ring_profile ('auto' | 'fit' | 'none'),
%               .calibrate_at ('measured' | 'desired'), .depths (desired) or
%               .measured_depths (selected measured depths), .predicted
%               (show the result at .depths)
%     'save'    .out_dir, .allow_shipped
%     'export_settings' .path, .settings;  'export_profiles' .path
%     'write_fit' .path (MAT file for a k-Wave drive check script)
%     'report'  .out_path ('' = <tempdir>/<combo>_report.html): the HTML report
%               of the saved library entry (GENERATE_CALIBRATION_REPORT)
%     'close'
%   settings - .forward_model, .near_field_skip_mm, .huber_delta, .n_starts,
%              .domain, .isppa_wcm2, .align, .max_extrapolation_mm, .tolerance,
%              .advanced (FIT_TRANSDUCER_MODEL option names → values) and
%              .ring_profile_options (FIT_RING_PROFILE options, incl. OffAxis: a
%              plane-scan file, and OffAxisOptions)
%
% Output:
%   out - struct; 'open' and 'load' return the combo (entry, its meta.fit,
%         measured depths, profile CSV, default desired depths); 'view' and
%         'fit' return the profiles, R² / focus and runtime settings per
%         depth of each model (geometric, library, fit) and summary messages
%
% See also: NEXUS_BRIDGE, CALIBRATION_GUI, CALIBRATION_FIT, CALIBRATION_EVALUATE

    persistent sessions
    if isempty(sessions); sessions = containers.Map(); end
    action = char(args.action);
    key = char(get_nonempty_field(args, 'session', 'default'));
    if strcmp(action, 'open')
        S = open_combo(args);
        sessions(key) = S;
        out = combo_info(S);
        return;
    end
    if ~isKey(sessions, key)
        error('nexus_calibration:noSession', 'No calibration session ''%s'': open a transducer first.', key);
    end
    S = sessions(key);
    switch action
        case 'load'
            S = load_profiles(S, char(args.csv));
            out = combo_info(S);
        case 'view'
            [S, out] = view(S, args, logical(get_nonempty_field(args, 'predicted', false)));
        case 'fit'
            [S, out] = do_fit(S, args);
        case 'save'
            out = do_save(S, args);
        case 'export_settings'
            out = export_settings(S, args);
        case 'export_profiles'
            out = export_profiles(S, args);
        case 'write_fit'
            out = write_fit(S, args);
        case 'report'
            out = do_report(S, args);
        case 'close'
            remove(sessions, key);
            out = struct('closed', key);
            return;
        otherwise
            error('nexus_calibration:action', 'Unknown action ''%s''.', action);
    end
    sessions(key) = S;
end

%% Actions ------------------------------------------------------------------
function S = open_combo(args)
    root = get_prestus_path();
    S.equip_path = char(get_nonempty_field(args, 'equip_path', fullfile(root, 'config', 'equipment')));
    S.lib_path   = fullfile(root, 'config', 'calibration');
    S.prof_root  = char(get_nonempty_field(args, 'profile_root', fullfile(root, 'local', 'profiles')));
    S.user_path  = char(get_nonempty_field(args, 'library_path', ''));
    S.user_libs  = [{S.user_path}, transducer_library_folders()];   % study folder, then PRESTUS_TRANSDUCER_LIBRARY
    S.user_libs  = S.user_libs(~cellfun(@isempty, S.user_libs));
    S.equip = load_equipment_config(S.equip_path, S.user_libs);   % + equipment entries in your own libraries
    S.serial = char(args.serial);
    if ~isfield(S.equip.trans, S.serial)
        error('nexus_calibration:serial', 'Transducer %s is not in %s.', S.serial, S.equip_path);
    end
    S.tran = S.equip.trans.(S.serial);
    S.ds = char(get_nonempty_field(args, 'driving_system', ''));
    S.combo = S.serial;
    if ~isempty(S.ds); S.combo = [S.serial '_' S.ds]; end
    S.matrix = ~isfield(S.tran.transducer, 'annular') && isfield(S.tran.transducer, 'matrix');
    S.lib = []; S.lib_file = '';
    for folder = [S.user_libs, {S.lib_path}]
        if isempty(folder{1}); continue; end
        [f, has] = library_entry_file(folder{1}, S.combo);
        if has
            S.lib = yaml.loadFile(f, 'ConvertToArray', true); S.lib_file = f; break;
        end
    end
    S.meas = []; S.csv = ''; S.fit = []; S.prof = []; S.predicted = false;
    csv = char(get_nonempty_field(args, 'csv', ''));
    if isempty(csv) && isfield(S.tran, 'combos') && isfield(S.tran.combos, S.ds) ...
            && isfield(S.tran.combos.(S.ds), 'char_data_path')
        csv = resolve_profile_file(char(S.tran.combos.(S.ds).char_data_path), S.equip, S.serial, ...
            S.user_libs, S.prof_root);
    end
    if ~isempty(csv) && ~S.matrix; S = load_profiles(S, csv); end
end

function S = load_profiles(S, csv)
    S.meas = read_axial_profiles(csv);
    S.csv = csv;
    S.fit = []; S.prof = []; S.predicted = false;
end

function out = combo_info(S)
    out = struct('combo', S.combo, 'serial', S.serial, 'driving_system', S.ds, 'matrix', S.matrix, ...
        'csv', S.csv, 'measured_depths', [], 'measured_span', [], 'library_file', S.lib_file, ...
        'calibration_type', '', 'last_updated', '', 'fit_record', [], 'default_depths', [], ...
        'dist_geom_mm', NaN, 'shipped_path', S.lib_path);
    if ~isempty(S.meas)
        out.measured_depths = S.meas.depths_ep_mm(:)';
        out.measured_span = [min(S.meas.depths_ep_mm), max(S.meas.depths_ep_mm)];
    end
    if ~isempty(S.lib) && isfield(S.lib, 'meta')
        out.calibration_type = char(string(get_nonempty_field(S.lib.meta, 'calibration_type', '')));
        out.last_updated = char(string(get_nonempty_field(S.lib.meta, 'last_updated', '')));
        out.fit_record = get_nonempty_field(S.lib.meta, 'fit', []);
    end
    if ~S.matrix
        g = S.tran.transducer.annular.dist_geom_ep_mm;
        out.dist_geom_mm = g;
        lo = 5 * ceil(get_nonempty_field(S.tran, 'min_foc', 0.5 * g) / 5);
        hi = 5 * floor(get_nonempty_field(S.tran, 'max_foc', 1.2 * g) / 5);
        if hi <= lo; hi = lo + 5; end
        out.default_depths = lo:5:hi;
    end
end

function [S, out] = view(S, args, predicted)
    assert(~S.matrix, 'nexus_calibration:matrix', ...
        '%s is a matrix transducer: prediction and fitting support annular arrays.', S.serial);
    st = settings(args);
    if predicted
        S.prof = calibration_desired_profiles(S.tran, S.meas, row(args.depths), st.profile);
    else
        assert(~isempty(S.meas), 'nexus_calibration:noProfiles', 'Load a profile CSV first.');
        S.prof = S.meas;
        S.prof.profile_type = repmat({'measured'}, 1, numel(S.meas.depths_ep_mm));
    end
    S.predicted = predicted;
    eo = struct('Predicted', predicted, 'Library', S.lib, 'Fit', S.fit, ...
        'FitOptions', st.fit, 'IsppaWcm2', st.isppa);
    if ~isempty(S.meas); eo.MeasuredZmaxMm = max(S.meas.z_ep_mm); end
    [res, msg] = calibration_evaluate(S.tran, S.prof, eo);
    out = payload(S, res, msg, st);
end

function [S, out] = do_fit(S, args)
    assert(~S.matrix, 'nexus_calibration:matrix', 'Fitting supports annular arrays only.');
    assert(~isempty(S.meas), 'nexus_calibration:noProfiles', 'Fitting needs measured profiles: load a CSV first.');
    st = settings(args);
    at_desired = strcmp(char(get_nonempty_field(args, 'calibrate_at', 'measured')), 'desired');
    if at_desired
        sel = calibration_desired_profiles(S.tran, S.meas, row(args.depths), st.profile);
    else
        d = row(get_nonempty_field(args, 'measured_depths', S.meas.depths_ep_mm));
        keep = arrayfun(@(x) any(abs(d - x) < 0.01), S.meas.depths_ep_mm);
        assert(any(keep), 'nexus_calibration:noDepths', 'Select at least one depth.');
        sel = interpolate_axial_profiles(S.meas, S.meas.depths_ep_mm(keep));
    end
    fo = struct('Mode', char(get_nonempty_field(args, 'mode', 'per_depth')), ...
        'FreeParams', {cellstr(get_nonempty_field(args, 'free_params', {'phase'; 'amp'}))'}, ...
        'FitOptions', st.fit, 'Library', S.lib, 'IsppaWcm2', st.isppa, 'TranSerial', S.serial, ...
        'DrivingSystem', S.ds, 'Source', S.csv, 'CalibrateAt', ternary(at_desired, 'desired', 'measured'), ...
        'Tool', 'nexus', 'RingProfile', char(get_nonempty_field(args, 'ring_profile', 'auto')), ...
        'RingProfileOptions', st.ring, 'Measured', S.meas, ...
        'Progress', @(k, n, d, pt) nexus_progress( ...
            sprintf('Depth %.1f mm, %s profile (%d/%d)', d, pt, k, n)));
    if at_desired
        fo.ProfileOptions = rmfield(st.profile, 'NearFieldSkipMm');
    end
    [S.fit, msg] = calibration_fit(S.tran, sel, fo);
    show_desired = at_desired || logical(get_nonempty_field(args, 'predicted', false));
    [S, out] = view(S, args, show_desired);
    out.messages = [{msg}, out.messages];
    out.fit_record = S.fit.record;
end

function out = do_save(S, args)
    assert(~isempty(S.fit), 'nexus_calibration:noFit', 'Nothing to save: run a fit first.');
    so = struct('ShippedLibraryPath', S.lib_path, ...
        'AllowShipped', logical(get_nonempty_field(args, 'allow_shipped', false)));
    [file, prof_csv] = calibration_save_entry(S.fit.entry, S.fit, char(args.out_dir), so);
    out = struct('file', file, 'profiles_csv', prof_csv, 'out_dir', char(args.out_dir));
end

function out = do_report(S, args)
    % The entry is looked up again: a 'save' since 'open' may have written one
    lib_file = ''; lib_folder = '';
    for folder = [S.user_libs, {S.lib_path}]
        if isempty(folder{1}); continue; end
        [f, has] = library_entry_file(folder{1}, S.combo);
        if has; lib_file = f; lib_folder = folder{1}; break; end
    end
    assert(~isempty(lib_file), 'nexus_calibration:noEntry', ...
        'No library entry for %s yet: save a fit first.', S.combo);
    out_path = char(get_nonempty_field(args, 'out_path', ''));
    if isempty(out_path); out_path = fullfile(tempdir, [S.combo '_report.html']); end
    file = generate_calibration_report(S.combo, lib_folder, out_path, S.equip);
    out = struct('file', file, 'library_file', lib_file);
end

function out = export_settings(S, args)
    assert(~isempty(S.prof), 'nexus_calibration:noView', 'Nothing on view: Predict or evaluate first.');
    st = settings(args);
    if ~isempty(S.fit)
        lib = S.fit.entry; what = 'fit';
    elseif ~isempty(S.lib)
        lib = S.lib; what = 'library';
    else
        lib = struct(); what = 'geometric';
    end
    so = struct('IsppaWcm2', st.isppa, 'CsvPath', char(args.path), 'ForwardModel', st.fit.ForwardModel, ...
        'NearFieldSkipMm', st.fit.NearFieldSkipMm, 'ProfileType', {S.prof.profile_type});
    T = library_settings_table(S.tran, lib, S.prof.depths_ep_mm, so);
    out = struct('path', char(args.path), 'what', what, 'n_depths', height(T));
end

function out = export_profiles(S, args)
    if ~isempty(S.fit)
        p = S.fit.prof; what = 'calibration';
    elseif S.predicted && any(~strcmp(S.prof.profile_type, 'unavailable'))
        p = S.prof; p.I = p.I(1:numel(S.meas.z_ep_mm), :); p.z_ep_mm = S.meas.z_ep_mm; what = 'desired';
    else
        error('nexus_calibration:noProfiles', ...
            'No profiles to export: load a CSV and Predict at desired depths, or fit.');
    end
    write_axial_profiles(p, char(args.path), struct('combo', S.combo));
    shown = any(isfinite(p.I), 1);
    out = struct('path', char(args.path), 'what', what, 'n_profiles', nnz(shown), ...
        'types', {unique(p.profile_type(shown))});
end

function out = write_fit(S, args)
% The unsaved per-depth fit, for a k-Wave drive check run as a script
    assert(~isempty(S.fit) && strcmp(S.fit.mode, 'per_depth'), 'nexus_calibration:noFit', ...
        'The k-Wave drive check needs a per-depth fit: run Fit first.');
    fit = S.fit; combo = S.combo; %#ok<NASGU>
    f = char(args.path);
    d = fileparts(f);
    if ~isempty(d) && ~isfolder(d); mkdir(d); end
    save(f, 'fit', 'combo');
    out = struct('path', f, 'depths_ep_mm', fit.prof.depths_ep_mm(:)');
end

%% Payload ------------------------------------------------------------------
function out = payload(S, res, msg, st)
    out = struct('combo', S.combo, 'predicted', S.predicted, ...
        'depths', S.prof.depths_ep_mm(:)', 'profile_type', {S.prof.profile_type(:)'}, ...
        'measured_span', [], 'dist_geom_mm', S.tran.transducer.annular.dist_geom_ep_mm, ...
        'messages', {msg(:)'}, 'results', struct(), 'measured_focus', [], 'params', [], ...
        'has_fit', ~isempty(S.fit), 'isppa_wcm2', st.isppa);
    if ~isempty(S.meas)
        out.measured_span = [min(S.meas.depths_ep_mm), max(S.meas.depths_ep_mm)];
    end
    for f = fieldnames(res)'
        r = res.(f{1});
        out.results.(f{1}) = struct('depths', r.depths(:)', 'z', r.z(:)', ...
            'Imod', {cols(r.Imod)}, 'Imeas', {cols(r.Imeas)}, 'r2', r.r2(:)', ...
            'focus', r.focus(:)', 'settings', {settings_rows(r.settings)});
    end
    if S.predicted && isfield(res, 'geometric')
        % Focus of the measured / derived profiles themselves (NaN without one)
        r = res.geometric;
        fm = nan(1, numel(r.depths));
        has = any(isfinite(r.Imeas), 1);
        if any(has)
            fm(has) = calibration_focus_beyond(r.Imeas(:, has), r.z, st.fit.NearFieldSkipMm);
        end
        out.measured_focus = fm;
    end
    if ~S.predicted && ~isempty(S.fit)
        m = S.fit.model;
        out.params = struct('kind', 'per_depth', 'depths', cellfun(@(x) x.fit_depths_ep_mm, m), ...
            'phase_deg', {cellfun(@(x) x.elem_phase_correction_deg(:)', m, 'UniformOutput', false)});
    end
end

function c = cols(M)
    if isempty(M); c = {}; else; c = num2cell(M, 1); end
end

function rows = settings_rows(s)
% Runtime settings per depth (LIBRARY_MODEL_AT_DEPTH) as plain rows
    rows = {};
    for k = 1:numel(s)
        lbl = char(s(k).phase_model);
        if strcmp(lbl, 'blend'); lbl = sprintf('blend w=%.2f', s(k).refinement_weight); end
        drive = NaN;
        if ~isempty(s(k).drive_pa); drive = s(k).drive_pa / 1e3; end
        rows{end + 1} = struct('phase_model', lbl, 'drive_kpa', drive, ...
            'elem_phase_deg', {num2cell(s(k).elem_phase_deg(:)')}); %#ok<AGROW>
    end
end

%% Helpers ------------------------------------------------------------------
function st = settings(args)
    s = get_nonempty_field(args, 'settings', struct());
    st.fit = struct('ForwardModel', char(get_nonempty_field(s, 'forward_model', 'rayleigh')), ...
        'NearFieldSkipMm', double(get_nonempty_field(s, 'near_field_skip_mm', 10)), ...
        'HuberDelta', double(get_nonempty_field(s, 'huber_delta', 0.8)), ...
        'NStarts', double(get_nonempty_field(s, 'n_starts', 16)), ...
        'Domain', char(get_nonempty_field(s, 'domain', 'intensity')));
    adv = get_nonempty_field(s, 'advanced', struct());
    if isstruct(adv)
        for f = fieldnames(adv)'
            v = adv.(f{1});
            if isnumeric(v); v = double(v(:)'); end
            st.fit.(f{1}) = v;
        end
    end
    st.isppa = double(get_nonempty_field(s, 'isppa_wcm2', 30));
    st.profile = struct('Align', char(get_nonempty_field(s, 'align', 'label')), ...
        'MaxExtrapolationMm', double(get_nonempty_field(s, 'max_extrapolation_mm', 0)), ...
        'NearFieldSkipMm', st.fit.NearFieldSkipMm, 'Tolerance', double(get_nonempty_field(s, 'tolerance', 0.01)));
    st.ring = plain(get_nonempty_field(s, 'ring_profile_options', struct()));    % FIT_RING_PROFILE options
end

function v = plain(v)
% JSON values as MATLAB options: numeric rows, char, nested structs
    if isstruct(v)
        for f = fieldnames(v)'; v.(f{1}) = plain(v.(f{1})); end
    elseif iscell(v) && all(cellfun(@isnumeric, v))
        v = cell2mat(v(:)');
    elseif isnumeric(v)
        v = double(v(:)');
    elseif isstring(v)
        v = char(v);
    end
end

function v = row(v)
    v = double(v(:)');
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end
