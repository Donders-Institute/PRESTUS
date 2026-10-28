function fig = calibration_gui(options)
% CALIBRATION_GUI  Standalone GUI: fit transducer library entries to hydrophone profiles
%
% Pick a transducer and driving system, load an axial-profile CSV, compare
% the measured profiles with geometric steering and the current library
% entry, fit per-depth refinements, optionally check them in free-water k-Wave (drive corrected to the simulated ISPPA;
% KWAVE_CALIBRATE_DRIVE), and save the result to a user library folder
% (transducer.library_path).
%
% Desired depths. "Predict" works at any depths: each gets a profile from
% the measurements (measured, or interpolated between the neighbouring
% measured depths, or — within "extrapolate ≤" — the nearest one shifted;
% INTERPOLATE_AXIAL_PROFILES), and geometric steering, the library entry and
% the unsaved fit are modelled there through the runtime library path
% (EVALUATE_LIBRARY_ENTRY: refinements interpolated between refined depths,
% blended back to the global model outside them). The view shows the
% profiles, R² against the (derived) profile, focus vs target, and per depth
% the phase model, drive and phases. With "Calibrate at: desired depths",
% Fit calibrates at the desired depths on those profiles; the result is
% shown in the same view, the profiles fitted to are written next to a
% saved entry, and both the profiles ("Export profiles…") and the per-depth
% settings ("Export settings…", LIBRARY_SETTINGS_TABLE) can be exported as CSV.
% Without a CSV only "Predict" is available (model only).
%
% "View geometry" draws the selected
% transducer; "New…" adds an annular or matrix transducer as {serial}.yaml
% in the equipment folder, with a live geometry preview, after
% CHECK_EQUIPMENT_CONFIG. Matrix transducers are geometry only: prediction
% and fitting support annular arrays. All computation goes through
% CALIBRATION_DESIRED_PROFILES, CALIBRATION_EVALUATE, CALIBRATION_FIT and
% CALIBRATION_SAVE_ENTRY (shared with the Calibration page of the PRESTUS
% Nexus) and through them FIT_TRANSDUCER_MODEL and the library helpers, so
% results match REFINE_TRANSDUCER_LIBRARY / BUILD_TRANSDUCER_LIBRARY exactly.
%
% A saved entry records how it was fitted (meta.fit: mode, free parameters,
% fit options, depths, profile CSV, reference ISPPA, k-Wave drive check).
% Loading an entry that has meta.fit restores those settings, so a refit
% starts where the last one ended.
%
% Use as:
%   calibration_gui()
%   calibration_gui(struct('ProfileRoot', '/path/to/Axial_profiles', ...
%                          'OutputPath',  '/path/to/my_library'))
%
% Input (options, all optional):
%   .EquipPath   - equipment folder (default config/equipment/); new
%                  transducers are saved here
%   .LibraryPath - shipped library (default config/calibration/)
%   .ProfileRoot - folder that char_data_path in the equipment YAML is relative to
%                  (default local/profiles/)
%   .OutputPath  - folder to save entries to (default: ask on save)
%   .SimConfig   - study/simulation YAML for the k-Wave drive calibration (grid, code_type,
%                  platform/hpc); default: config_default
%   .ExportFolder - write "Export profiles…" / "Export settings…" CSVs here
%                  under their default names, without a file dialog
%
% Output:
%   fig - uifigure handle
%
% See also: CALIBRATION_FIT, CALIBRATION_EVALUATE, CALIBRATION_SAVE_ENTRY,
%           FIT_TRANSDUCER_MODEL, EVALUATE_LIBRARY_ENTRY, INTERPOLATE_AXIAL_PROFILES,
%           LIBRARY_SETTINGS_TABLE, REFINE_TRANSDUCER_LIBRARY, LOAD_TRANSDUCER_FROM_LIBRARY,
%           CHECK_EQUIPMENT_CONFIG, KWAVE_CALIBRATE_DRIVE

    if nargin < 1; options = struct(); end
    root = get_prestus_path();
    st.equip_path = opt(options, 'EquipPath',   fullfile(root, 'config', 'equipment'));
    st.lib_path   = opt(options, 'LibraryPath', fullfile(root, 'config', 'calibration'));
    st.prof_root  = opt(options, 'ProfileRoot', fullfile(root, 'local', 'profiles'));
    st.out_path   = opt(options, 'OutputPath',  '');
    st.export_dir = opt(options, 'ExportFolder', '');
    st.equip      = load_equipment_config(st.equip_path);
    st.meas = [];           % measured profiles (READ_AXIAL_PROFILES), [] = none loaded
    st.prof = [];           % profiles on view: st.meas, or the desired depths of Predict with
                            % their measured/derived profiles (.profile_type)
    st.lib = []; st.res = struct();
    st.fit = [];  % current (unsaved) fit: CALIBRATION_FIT result (.mode, .model, .prof, .entry, .record)
    st.kw = [];   % k-Wave drive calibration of the current per-depth fit: .lib, .report, .sims, .settings
    st.predicted = false;   % true: st.prof is the desired-depth grid
    st.adv = struct();      % advanced FIT_TRANSDUCER_MODEL options ("Advanced…"); empty = defaults
    st.pext = struct();     % desired-profile options ("Advanced…": Tolerance)
    st.kwo = struct();      % k-Wave drive calibration options ("Advanced…": GridDims, ApplyCorrection)
    st.rpo = struct();      % FIT_RING_PROFILE / plane-scan options ("Ring profile" Options…)

    serials = transducer_serials(st.equip);

    %% Layout ---------------------------------------------------------------
    fig = uifigure('Name', 'PRESTUS transducer calibration', 'Position', [80 80 1350 860]);
    g   = uigridlayout(fig, [1 2], 'ColumnWidth', {340, '1x'});
    L   = uigridlayout(g, [31 2], 'RowHeight', repmat({24}, 1, 31), 'ColumnWidth', {130, '1x'}, ...
                       'RowSpacing', 4, 'Scrollable', 'on');
    R   = uigridlayout(g, [3 1], 'RowHeight', {'1.3x', '1x', 170});

    function h = at(h, r, c)
        h.Layout.Row = r; h.Layout.Column = c;
    end
    lab = @(txt, r) at(uilabel(L, 'Text', txt), r, 1);

    lab('Transducer', 1);
    ui.tran = at(uidropdown(L, 'Items', serials, 'ValueChangedFcn', @(~,~) on_transducer()), 1, 2);
    tb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 2, 2);
    uibutton(tb, 'Text', 'View geometry', 'ButtonPushedFcn', @(~,~) on_view_transducer());
    uibutton(tb, 'Text', 'New…', 'Tooltip', 'Add a transducer to the equipment folder', ...
             'ButtonPushedFcn', @(~,~) on_new_transducer());
    lab('Driving system', 3);
    ui.ds   = at(uidropdown(L, 'Items', {''}, 'ValueChangedFcn', @(~,~) on_combo()), 3, 2);
    lab('Profile CSV', 4);
    ui.csv  = at(uieditfield(L, 'text'), 4, 2);
    b = at(uigridlayout(L, [1 2], 'Padding', 0), 5, 2);
    uibutton(b, 'Text', 'Browse…', 'ButtonPushedFcn', @(~,~) on_browse_csv());
    uibutton(b, 'Text', 'Load', 'ButtonPushedFcn', @(~,~) on_load());

    lab('Fit depths (measured)', 6);
    ui.depths = at(uilistbox(L, 'Items', {}, 'Multiselect', 'on'), [6 9], 2);
    lab('Desired depths (mm)', 10);
    pb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 60}, 'ColumnSpacing', 4), 10, 2);
    ui.pred = uieditfield(pb, 'text', 'Tooltip', 'Depths from the exit plane, e.g. 20:5:60 or 30, 42.5, 50');
    uibutton(pb, 'Text', 'Predict', 'Tooltip', ...
             ['Profiles at these depths (measured, or interpolated from the measured ones) and the ' ...
              'geometric, library and fitted models there as a study resolves them'], ...
             'ButtonPushedFcn', @(~,~) on_predict());
    lab('Desired profiles', 11);
    ip = at(uigridlayout(L, [1 3], 'Padding', 0, 'ColumnWidth', {'1x', 78, 42}, 'ColumnSpacing', 4), 11, 2);
    ui.align = uidropdown(ip, 'Items', {'label', 'peak'}, 'Tooltip', ...
        ['Focus alignment when interpolating between measured depths: column focal-depth ' ...
         'labels (nominal setpoints, default) or intensity peaks beyond the near-field skip']);
    uilabel(ip, 'Text', 'extrapolate ≤', 'HorizontalAlignment', 'right');
    ui.extrap = uieditfield(ip, 'numeric', 'Value', 0, 'Limits', [0 50], 'Tooltip', ...
        ['Desired depths up to this far outside the measured range [mm] get the nearest measured ' ...
         'profile shifted to the desired depth; 0 = no extrapolation']);
    lab('Calibrate at', 12);
    ui.calat = at(uidropdown(L, 'Items', {'Fit depths (measured)', 'Desired depths'}, 'Tooltip', ...
        ['Fit depths: the measured depths selected above. Desired depths: the desired depths, ' ...
         'on their measured or interpolated/extrapolated profiles']), 12, 2);

    lab('Ring profile', 13);
    rpb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 80}, 'ColumnSpacing', 4), 13, 2);
    ui.rpo_btn = uibutton(rpb, 'Text', 'Options…', 'Tooltip', ['Sub-rings, smoothness, optimiser and ' ...
        'plane scans (off-axis data) for the ring-profile fit'], 'ButtonPushedFcn', @(~,~) on_ring_options());
    ui.ring = at(uidropdown(rpb, 'Items', {'keep the entry''s', 'fit (shared by all depths)', 'uniform rings'}, ...
        'ItemsData', {'auto', 'fit', 'none'}, 'Tooltip', ...
        ['How each ring vibrates across its width (FIT_RING_PROFILE). Fit: every ring split into ' ...
         'sub-rings with weights shared by all measured depths, then phases, gains and drive per depth; ' ...
         'the near field from 4 mm is then fitted too. Improves transducers whose near-field / focal ' ...
         'peak ratio uniform rings cannot reproduce (e.g. CTX). Axial profiles alone do not determine ' ...
         'the off-axis field: add plane scans under Options… where available.']), 1, 1);
    lab('Free parameters', 14);
    fp = at(uigridlayout(L, [2 2], 'Padding', 0, 'RowHeight', {20, 20}, 'RowSpacing', 4), [14 15], 2);
    ui.fp_phase = uicheckbox(fp, 'Text', 'phase', 'Value', true);
    ui.fp_amp   = uicheckbox(fp, 'Text', 'gain', 'Value', true);
    lab('Forward model', 16);
    ui.fwd     = at(uidropdown(L, 'Items', {'rayleigh', 'oneil'}), 16, 2);
    lab('Near-field skip (mm)', 17);
    ui.skip    = at(uieditfield(L, 'numeric', 'Value', 10, 'Limits', [0 50]), 17, 2);
    lab('Loss on', 18);
    ui.domain  = at(uidropdown(L, 'Items', {'intensity', 'amplitude'}), 18, 2);
    lab('Huber δ', 19);
    ui.huber   = at(uieditfield(L, 'numeric', 'Value', 0.8, 'Limits', [1e-3 100], 'Tooltip', ...
        'Pseudo-Huber scale on the peak-normalised profile; large (e.g. 10) = least squares'), 19, 2);
    lab('Starts', 20);
    nsb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 96}, 'ColumnSpacing', 4), 20, 2);
    ui.nstarts = uieditfield(nsb, 'numeric', 'Value', 16, 'Limits', [1 200], 'RoundFractionalValues', 'on');
    ui.adv_btn = uibutton(nsb, 'Text', 'Advanced…', 'Tooltip', ['Phase constraints, absolute-intensity fit, ' ...
        'global optimisers, focal weighting (PERFORM_GLOBAL_SEARCH behaviour)'], ...
        'ButtonPushedFcn', @(~,~) on_advanced());
    lab('Reference ISPPA (W/cm²)', 21);
    ui.isppa   = at(uieditfield(L, 'numeric', 'Value', 30, 'Limits', [0.1 1000]), 21, 2);

    lab('k-Wave drive', 22);
    kb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 95}, 'ColumnSpacing', 4), 22, 2);
    ui.kw      = uicheckbox(kb, 'Text', 'calibrate in free water, after fit', 'Value', false, 'Tooltip', ...
        ['Per-depth only: simulate each refined depth in free-water k-Wave through the ' ...
         'runtime library path, compare with the (measured or derived) profile and correct the drive ' ...
         'to the simulated ISPPA (minutes per depth).']);
    ui.kw_grid = uidropdown(kb, 'Items', {'3D', 'axisymmetric'}, ...
        'Tooltip', 'Free-water grid, unless the sim config sets grid.axisymmetric');
    lab('Sim config (k-Wave)', 23);
    cb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 60}, 'ColumnSpacing', 4), 23, 2);
    ui.simcfg = uieditfield(cb, 'text', 'Value', opt(options, 'SimConfig', ''), 'Tooltip', ...
        'Study/simulation YAML for grid, code_type and platform/hpc; empty = config_default');
    uibutton(cb, 'Text', 'Browse…', 'ButtonPushedFcn', @(~,~) on_browse_simcfg());

    ab = at(uigridlayout(L, [1 2], 'Padding', 0), 24, [1 2]);
    uibutton(ab, 'Text', 'Evaluate measured', 'Tooltip', ...
             'Geometric, library and fit at the measured depths (R² against the hydrophone profiles)', ...
             'ButtonPushedFcn', @(~,~) on_evaluate_measured());
    uibutton(ab, 'Text', 'Fit', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) on_fit());
    lab('Save to folder', 25);
    ui.out = at(uieditfield(L, 'text', 'Value', st.out_path), 25, 2);
    sb = at(uigridlayout(L, [1 2], 'Padding', 0), 26, [1 2]);
    uibutton(sb, 'Text', 'Browse…', 'ButtonPushedFcn', @(~,~) on_browse_out());
    uibutton(sb, 'Text', 'Save entry', 'ButtonPushedFcn', @(~,~) on_save());
    eb = at(uigridlayout(L, [1 2], 'Padding', 0), 27, [1 2]);
    uibutton(eb, 'Text', 'Export profiles…', 'Tooltip', ...
             ['Profile CSV + provenance sidecar: the profiles of the last fit, else those at the ' ...
              'desired depths (measured and interpolated/extrapolated)'], ...
             'ButtonPushedFcn', @(~,~) on_export_profiles());
    uibutton(eb, 'Text', 'Export settings…', 'Tooltip', ...
             ['CSV of the phases, drive and phase model per depth that the fit (else the library entry) ' ...
              'resolves to, at the depths on view'], ...
             'ButtonPushedFcn', @(~,~) on_export_settings());
    ui.status = at(uitextarea(L, 'Editable', 'off', 'Value', {'Select a transducer and load profiles.'}), [28 31], [1 2]);

    % Right: profile view, R² by depth, parameters
    pv = uigridlayout(R, [2 1], 'RowHeight', {24, '1x'}, 'Padding', 0);
    pvh = uigridlayout(pv, [1 2], 'ColumnWidth', {90, 140}, 'Padding', 0);
    uilabel(pvh, 'Text', 'Show depth');
    ui.show = uidropdown(pvh, 'Items', {''}, 'ValueChangedFcn', @(~,~) draw_profile());
    ui.ax_prof = uiaxes(pv);
    ui.ax_r2   = uiaxes(R);
    ui.params  = uitable(R, 'ColumnName', {'Element', 'Phase offset (°)', 'Gain'}, 'RowName', {});

    on_transducer();

    %% Callbacks --------------------------------------------------------------
    function on_browse_simcfg()
        [f, p] = uigetfile({'*.yaml;*.yml', 'YAML config'}, 'Simulation config for the k-Wave drive calibration');
        if ~isequal(f, 0); ui.simcfg.Value = fullfile(p, f); end
    end

    function on_transducer()
        tran = st.equip.trans.(ui.tran.Value);
        if isfield(tran, 'combos') && ~isempty(tran.combos)
            ui.ds.Items = fieldnames(tran.combos);
        else
            ui.ds.Items = {'(none)'};
        end
        ui.pred.Value = default_predict_depths(tran);
        on_combo();
    end

    function on_combo()
        tran = st.equip.trans.(ui.tran.Value);
        ds   = ui.ds.Value;
        ui.csv.Value = '';
        st.meas = []; st.prof = []; st.predicted = false; st.res = struct(); st.fit = []; st.kw = [];
        ui.depths.Items = {};
        if is_matrix(tran)
            st.lib = [];
            ui.show.Items = {''};
            draw();
            say(sprintf(['%s: matrix transducer, geometry only (View geometry). ' ...
                'Prediction and fitting support annular arrays.'], ui.tran.Value));
            return;
        end
        if isfield(tran, 'combos') && isfield(tran.combos, ds) && isfield(tran.combos.(ds), 'char_data_path')
            cand = char(tran.combos.(ds).char_data_path);
            if ~is_absolute(cand); cand = fullfile(st.prof_root, cand); end
            if isfile(cand); ui.csv.Value = cand; end
        end
        st.lib = load_entry();
        say(sprintf('%s: %s', combo_name(), ternary(isempty(st.lib), 'no library entry', ...
            sprintf('library entry (%s)', char(string(get_field(get_field(st.lib, 'meta', []), 'calibration_type', '?')))))));
        apply_fit_record(st.lib);
        if ~isempty(ui.csv.Value); on_load(); else; on_predict(); end
    end

    function on_browse_csv()
        [f, p] = uigetfile('*.csv', 'Axial profile CSV', st.prof_root);
        if isequal(f, 0); return; end
        ui.csv.Value = fullfile(p, f); on_load();
    end

    function on_browse_out()
        p = uigetdir(ui.out.Value, 'Folder for the library entry');
        if ~isequal(p, 0); ui.out.Value = p; end
    end

    function on_load()
        try
            st.meas = read_axial_profiles(ui.csv.Value);
        catch ME
            say(['Could not read profiles: ' ME.message]); return;
        end
        items = depth_items(st.meas.depths_ep_mm);
        ui.depths.Items = items; ui.depths.Value = items;
        st.fit = []; st.kw = [];
        say(sprintf('Loaded %d depths (%.1f–%.1f mm) from %s', numel(items), ...
            min(st.meas.depths_ep_mm), max(st.meas.depths_ep_mm), ui.csv.Value));
        show_measured();
        on_evaluate();
    end

    function show_measured()
        % View the measured depths (fit view)
        st.prof = st.meas; st.prof.profile_type = repmat({'measured'}, 1, numel(st.meas.depths_ep_mm));
        st.predicted = false;
        items = depth_items(st.meas.depths_ep_mm);
        prev = ui.show.Value;
        ui.show.Items = items;
        if any(strcmp(items, prev)); ui.show.Value = prev; else; ui.show.Value = items{1}; end
    end

    function on_evaluate_measured()
        if isempty(st.meas); on_predict(); return; end
        show_measured();
        on_evaluate();
    end

    function [dep, ok] = desired_depths()
        [dep, ok] = parse_depths(ui.pred.Value);
        ok = ok && ~isempty(dep);
        if ok; dep = unique(round(dep(:)', 2)); end
    end

    function o = profile_options()
        % INTERPOLATE_AXIAL_PROFILES options for profiles at desired depths
        o = struct('Align', ui.align.Value, 'MaxExtrapolationMm', ui.extrap.Value, ...
                   'NearFieldSkipMm', ui.skip.Value);
        if isfield(st.pext, 'Tolerance'); o.Tolerance = st.pext.Tolerance; end
    end

    function on_predict()
        % Desired depths. The modelled profile shape depends only on geometry,
        % frequency and the resolved phases/gains, so no measurement is
        % needed; where profiles are loaded, each desired depth also gets its
        % measured or derived profile to compare with (and calibrate on).
        tran = st.equip.trans.(ui.tran.Value);
        if is_matrix(tran); say('Prediction supports annular arrays only.'); return; end
        [dep, ok] = desired_depths();
        if ~ok; say('Desired depths: enter depths such as 20:5:60 or 30, 42.5, 50.'); return; end
        st.prof = calibration_desired_profiles(tran, st.meas, dep, profile_options());
        st.predicted = true;
        items = depth_items(st.prof.depths_ep_mm);
        ui.show.Items = [{'All depths'}, items]; ui.show.Value = 'All depths';
        on_evaluate();
    end

    function on_evaluate()
        % Rebuild all results on view from the current state (profiles,
        % library entry, fit, k-Wave drive calibration) and redraw
        if is_matrix(st.equip.trans.(ui.tran.Value)); say('Evaluation supports annular arrays only.'); return; end
        if isempty(st.prof); on_predict(); return; end
        tran = st.equip.trans.(ui.tran.Value);
        eo = struct('Predicted', st.predicted, 'Library', st.lib, 'Fit', st.fit, ...
            'FitOptions', base_options(), 'IsppaWcm2', ui.isppa.Value);
        if ~isempty(st.fit); eo.FitEntry = fitted_entry(); end
        if ~isempty(st.kw); eo.Kwave = st.kw; end
        if ~isempty(st.meas); eo.MeasuredZmaxMm = max(st.meas.z_ep_mm); end
        [st.res, msg] = calibration_evaluate(tran, st.prof, eo);
        cellfun(@say, msg);
        draw();
    end

    function tf = fit_on_derived()
        tf = ~isempty(st.fit) && any(~strcmp(st.fit.prof.profile_type, 'measured'));
    end

    function on_fit()
        if isempty(st.meas); say('Fitting needs measured profiles: load a CSV first.'); return; end
        tran = st.equip.trans.(ui.tran.Value);
        at_desired = strcmp(ui.calat.Value, 'Desired depths');
        if at_desired
            [dep, ok] = desired_depths();
            if ~ok; say('Desired depths: enter depths such as 20:5:60 or 30, 42.5, 50.'); return; end
            sel  = calibration_desired_profiles(tran, st.meas, dep, profile_options());
            none = strcmp(sel.profile_type, 'unavailable');
            if any(none)
                say(sprintf(['No profile at %s mm (outside the measured range %.1f–%.1f mm; ' ...
                    'raise "extrapolate ≤" to include them): skipped.'], mat2str(sel.depths_ep_mm(none)), ...
                    min(st.meas.depths_ep_mm), max(st.meas.depths_ep_mm)));
            end
            if all(none); say('None of the desired depths has a profile to calibrate on.'); return; end
        else
            sel = selected_profiles();
            if isempty(sel.depths_ep_mm); say('Select at least one depth.'); return; end
        end
        % No z shift: it is not applied at run time, so phases fitted with it
        % would match the measurement only after a shift the entry never gets
        free = {'phase', 'amp'};
        free = free([ui.fp_phase.Value, ui.fp_amp.Value]);
        dlg  = uiprogressdlg(fig, 'Title', 'Fitting', 'Indeterminate', 'on');
        cleanup = onCleanup(@() close(dlg));
        fo = struct('FreeParams', {free}, 'FitOptions', base_options(), ...
            'Library', st.lib, 'IsppaWcm2', ui.isppa.Value, 'TranSerial', ui.tran.Value, ...
            'DrivingSystem', ternary(strcmp(ui.ds.Value, '(none)'), '', ui.ds.Value), ...
            'Source', ui.csv.Value, 'CalibrateAt', ternary(at_desired, 'desired', 'measured'), ...
            'Tool', 'calibration_gui', 'RingProfile', ui.ring.Value, 'RingProfileOptions', ring_options(), ...
            'Measured', st.meas, ...
            'Progress', @(k, n, d, pt) set_msg(dlg, sprintf('Depth %.1f mm, %s profile (%d/%d)', d, pt, k, n)));
        if at_desired
            po = profile_options(); fo.ProfileOptions = rmfield(po, 'NearFieldSkipMm');
        end
        st.kw = [];
        try
            [st.fit, msg] = calibration_fit(tran, sel, fo);
        catch ME
            say(ME.message); return;
        end
        say(msg);
        if ~isempty(st.fit.ring_profile)
            plot_ring_profile(tran, st.fit.entry.global_model.ring_profile, st.fit.ring_profile);
        end
        if ui.kw.Value
            run_kwave_drive_calibration(st.fit.prof, dlg);
        end
        if at_desired || st.predicted
            on_predict();    % show the result on the desired depths, as a study resolves it
        else
            show_measured(); on_evaluate();
        end
    end

    function run_kwave_drive_calibration(sel, dlg)
        % Free-water k-Wave of the refined depths through the runtime library path
        dlg.Indeterminate = 'off'; dlg.Title = 'k-Wave drive calibration';
        ko = struct('EquipParam', st.equip, 'Axisymmetric', strcmp(ui.kw_grid.Value, 'axisymmetric'), ...
            'NearFieldSkipMm', ui.skip.Value, 'WorkDir', fullfile(tempdir, 'prestus_kwave_drive', combo_name()), ...
            'ProgressFcn', @(k, n, d) set(dlg, 'Value', (k - 1) / n, ...
                'Message', sprintf('Free-water simulation at %.1f mm (%d/%d)', d, k, n)));
        if ~isempty(strtrim(ui.simcfg.Value)); ko.SimConfig = strtrim(ui.simcfg.Value); end
        for f = fieldnames(st.kwo)'; ko.(f{1}) = st.kwo.(f{1}); end      % GridDims, ApplyCorrection
        t1 = tic;
        try
            [lib, rep, sims] = kwave_calibrate_drive(st.fit.entry, sel, ko);
        catch ME
            say(['k-Wave drive calibration failed: ' ME.message]); return;
        end
        st.kw = struct('lib', lib, 'report', rep, 'sims', sims, 'settings', struct( ...
            'grid', ui.kw_grid.Value, 'sim_config', strtrim(ui.simcfg.Value)));
        say(sprintf(['k-Wave drive calibration (%d depths, %.0f s): ISPPA %s W/cm² at the analytical drive ' ...
            '(target %g) → drive ×%s; median R² %.3f'], height(rep), toc(t1), ...
            strjoin(compose('%.1f', rep.isppa_kwave_wcm2'), ' '), ui.isppa.Value, ...
            strjoin(compose('%.2f', rep.amp_scale'), ' '), median(rep.r2_kwave)));
    end

    function on_save()
        if isempty(st.fit); say('Nothing to save — run a fit first.'); return; end
        out = ui.out.Value;
        if isempty(out); on_browse_out(); out = ui.out.Value; end
        if isempty(out); return; end
        so = struct('ShippedLibraryPath', st.lib_path);
        if strcmp(canonical(out), canonical(st.lib_path))
            sel = uiconfirm(fig, ['This overwrites the SHIPPED library entry. ' ...
                'Use a user folder and transducer.library_path unless you maintain the library.'], ...
                'Shipped library', 'Options', {'Overwrite shipped entry', 'Cancel'}, 'DefaultOption', 2);
            if ~strcmp(sel, 'Overwrite shipped entry'); return; end
            so.AllowShipped = true;
        end
        if ~isempty(st.kw); so.Kwave = st.kw.settings; end
        try
            [f, pcsv] = calibration_save_entry(fitted_entry(), st.fit, out, so);
        catch ME
            say(['Could not save: ' ME.message]); return;
        end
        if ~isempty(pcsv); say(sprintf('Calibration profiles (measured + derived): %s', pcsv)); end
        say(sprintf('Saved %s — set transducer.library_path: %s', f, out));
    end

    function on_export_profiles()
        % The profiles of the last fit, else the desired-depth profiles on view
        if ~isempty(st.fit)
            p = st.fit.prof; what = 'calibration';
        elseif st.predicted && any(~strcmp(st.prof.profile_type, 'unavailable'))
            p = st.prof; p.I = p.I(1:numel(st.meas.z_ep_mm), :); p.z_ep_mm = st.meas.z_ep_mm; what = 'desired';
        else
            say('No profiles to export: load a CSV and Predict at desired depths, or fit.'); return;
        end
        f = ask_path(sprintf('%s_%s_profiles.csv', combo_name(), what), 'Export axial profiles');
        if isempty(f); return; end
        write_axial_profiles(p, f, struct('combo', combo_name()));
        [d, stem] = fileparts(f);
        say(sprintf('Exported %d profiles (%s) to %s (+ %s.yaml)', nnz(any(isfinite(p.I), 1)), ...
            strjoin(unique(p.profile_type(any(isfinite(p.I), 1))), ', '), f, fullfile(d, stem)));
    end

    function on_export_settings()
        % Per-depth runtime settings of the fit (else the library entry) at the depths on view
        if isempty(st.prof); say('Nothing on view: load a CSV or Predict first.'); return; end
        tran = st.equip.trans.(ui.tran.Value);
        if ~isempty(st.fit)
            lib = fitted_entry(); what = 'fit';
        elseif ~isempty(st.lib)
            lib = st.lib; what = 'library';
        else
            lib = struct(); what = 'geometric';
        end
        f = ask_path(sprintf('%s_%s_settings.csv', combo_name(), what), 'Export settings per depth');
        if isempty(f); return; end
        o = base_options();
        so = struct('IsppaWcm2', ui.isppa.Value, 'CsvPath', f, 'ForwardModel', o.ForwardModel, ...
                    'NearFieldSkipMm', o.NearFieldSkipMm, 'ProfileType', {st.prof.profile_type});
        T = library_settings_table(tran, lib, st.prof.depths_ep_mm, so);
        say(sprintf('Exported settings of the %s at %d depths to %s', what, height(T), f));
    end

    function f = ask_path(name, ttl)
        if ~isempty(st.export_dir)                      % options.ExportFolder: no dialog
            f = fullfile(st.export_dir, name); return;
        end
        start = ui.out.Value;
        if isempty(start) || ~isfolder(start); start = pwd; end
        [fn, p] = uiputfile('*.csv', ttl, fullfile(start, name));
        figure(fig);
        if isequal(fn, 0); f = ''; else; f = fullfile(p, fn); end
    end

    function lib = fitted_entry()
        % The entry Save would write: the per-depth fit with k-Wave-corrected
        % drives if checked, else the analytical fit (CALIBRATION_FIT_ENTRY)
        if ~isempty(st.kw); lib = st.kw.lib; else; lib = st.fit.entry; end
    end

    function on_view_transducer()
        tran = st.equip.trans.(ui.tran.Value);
        v  = uifigure('Name', sprintf('Geometry: %s', ui.tran.Value), ...
                      'Position', [fig.Position(1:2) + [80 80], 760 720]);
        vg = uigridlayout(v, [3 1], 'RowHeight', {'1.4x', '1x', 60});
        ax1 = uiaxes(vg); ax2 = uiaxes(vg);
        info = draw_transducer_geometry(ax1, ax2, tran);
        uitextarea(vg, 'Editable', 'off', 'Value', info);
    end

    function on_new_transducer()
        % Form for a new annular or matrix transducer with a live geometry
        % preview; writes {serial}.yaml to st.equip_path
        d = uifigure('Name', 'New transducer', 'WindowStyle', 'modal', ...
                     'Position', [fig.Position(1:2) + [40 10], 1200 800]);
        top = uigridlayout(d, [1 2], 'ColumnWidth', {470, '1x'});
        G = uigridlayout(top, [24 2], 'RowHeight', [repmat({24}, 1, 23), {'1x'}], ...
                         'ColumnWidth', {180, '1x'}, 'RowSpacing', 4, 'Padding', 0);
        P = uigridlayout(top, [2 1], 'RowHeight', {'1.4x', '1x'}, 'Padding', 0);
        nt.ax_side = uiaxes(P); nt.ax_face = uiaxes(P);
        upd  = @(~,~) preview();
        dl   = @(par, txt, r) at(uilabel(par, 'Text', txt), r, 1);
        num  = @(par, v, lim) uieditfield(par, 'numeric', 'Value', v, 'Limits', lim, 'ValueChangedFcn', upd);
        inum = @(par, v, lim) uieditfield(par, 'numeric', 'Value', v, 'Limits', lim, ...
                                          'RoundFractionalValues', 'on', 'ValueChangedFcn', upd);
        sub  = @(par, n, r) at(uigridlayout(par, [1 n], 'Padding', 0, 'ColumnSpacing', 4), r, 2);

        % Common fields
        dl(G, 'Copy from', 1);
        nt.tmpl = at(uidropdown(G, 'Items', [{'(blank)'}; transducer_serials(st.equip)], ...
                    'ValueChangedFcn', @(~,~) fill_from_template()), 1, 2);
        dl(G, 'Type', 2);
        nt.type = at(uidropdown(G, 'Items', {'annular', 'matrix'}, 'ValueChangedFcn', @(~,~) on_type()), 2, 2);
        dl(G, 'Serial', 3);                  nt.serial = at(uieditfield(G, 'text', 'ValueChangedFcn', upd), 3, 2);
        dl(G, 'Name', 4);                    nt.name   = at(uieditfield(G, 'text'), 4, 2);
        dl(G, 'Manufacturer', 5);            nt.manu   = at(uieditfield(G, 'text'), 5, 2);
        dl(G, 'Frequency (kHz)', 6);         nt.freq   = at(num(G, 250, [1 1e4]), 6, 2);
        dl(G, 'Curvature radius (mm)', 7);   nt.rc     = at(num(G, 63, [1 1000]), 7, 2);
        dl(G, 'Exit plane → geom. focus (mm)', 8);
        nt.dep = at(num(G, 52, [0 1000]), 8, 2);
        nt.dep.Tooltip = 'dist_geom_ep_mm: distance from the exit plane to the centre of curvature';
        dl(G, 'Focal range (mm from EP)', 9);
        fr = sub(G, 2, 9);
        nt.fmin = num(fr, 10, [0 1000]); nt.fmax = num(fr, 60, [0 1000]);

        % Annular geometry
        A = uigridlayout(G, [7 2], 'RowHeight', repmat({24}, 1, 7), 'ColumnWidth', {180, '1x'}, ...
                         'RowSpacing', 4, 'Padding', 0);
        A.Layout.Row = [10 16]; A.Layout.Column = [1 2];
        dl(A, 'Physical channels', 1);       nt.nch = at(inum(A, 4, [1 1024]), 1, 2);
        dl(A, 'Ring inner Ø (mm)', 2);
        nt.id = at(uieditfield(A, 'text', 'Placeholder', 'e.g. 0, 30.2, 42.1, 51.1', 'ValueChangedFcn', upd), 2, 2);
        dl(A, 'Ring outer Ø (mm)', 3);
        nt.od = at(uieditfield(A, 'text', 'Placeholder', 'e.g. 29.6, 41.6, 50.6, 57.9', 'ValueChangedFcn', upd), 3, 2);
        at(uilabel(A, 'Text', 'One value per ring, comma- or space-separated; rings group into channels in order.', ...
           'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4]), [4 5], [1 2]);

        % Matrix geometry
        M = uigridlayout(G, [7 2], 'RowHeight', repmat({24}, 1, 7), 'ColumnWidth', {180, '1x'}, ...
                         'RowSpacing', 4, 'Padding', 0, 'Visible', 'off');
        M.Layout.Row = [10 16]; M.Layout.Column = [1 2];
        dl(M, 'Steering', 1);
        ms = sub(M, 2, 1);
        nt.m_steer  = uidropdown(ms, 'Items', {'3D', '1D'}, 'ValueChangedFcn', upd);
        nt.m_curved = uicheckbox(ms, 'Text', 'curved (spherical)', 'Value', true, 'ValueChangedFcn', upd);
        dl(M, 'Aperture Ø (mm)', 2);         nt.m_D = at(num(M, 60, [1 1000]), 2, 2);
        dl(M, 'Element shape', 3);
        nt.m_shape = at(uidropdown(M, 'Items', {'rect', 'disc', 'bowl'}, 'ValueChangedFcn', upd), 3, 2);
        dl(M, 'Element width × height (mm)', 4);
        mw = sub(M, 2, 4);
        nt.m_w = num(mw, 3, [0.01 100]); nt.m_h = num(mw, 3, [0.01 100]);
        dl(M, 'Grid layout', 5);
        nt.m_grid = at(uidropdown(M, 'Items', {'rect', 'fibonacci', 'fermat'}, ...
                       'ValueChangedFcn', @(~,~) on_grid()), 5, 2);
        GR = uigridlayout(M, [2 2], 'RowHeight', {24, 24}, 'ColumnWidth', {180, '1x'}, 'RowSpacing', 4, 'Padding', 0);
        GR.Layout.Row = [6 7]; GR.Layout.Column = [1 2];
        dl(GR, 'Rows × columns', 1);
        g1 = sub(GR, 2, 1);
        nt.r_nrow = inum(g1, 16, [1 256]); nt.r_ncol = inum(g1, 16, [1 256]);
        dl(GR, 'Spacing h, w (mm); sparsity', 2);
        g2 = sub(GR, 3, 2);
        nt.r_sh = num(g2, 0.5, [0 100]); nt.r_sw = num(g2, 0.5, [0 100]); nt.r_sp = num(g2, 1, [0.01 1]);
        set([nt.r_sh, nt.r_sw], 'LowerLimitInclusive', 'off');   % validate_matrix_transducer requires > 0
        GF = uigridlayout(M, [1 2], 'ColumnWidth', {180, '1x'}, 'Padding', 0, 'Visible', 'off');
        GF.Layout.Row = 6; GF.Layout.Column = [1 2];
        dl(GF, 'Elements; kerf (mm)', 1);
        g3 = sub(GF, 2, 1);
        nt.f_n = inum(g3, 256, [1 4096]); nt.f_kerf = num(g3, 0.5, [0 100]);
        GM = uigridlayout(M, [1 2], 'ColumnWidth', {180, '1x'}, 'Padding', 0, 'Visible', 'off');
        GM.Layout.Row = 6; GM.Layout.Column = [1 2];
        dl(GM, 'Elements', 1);               nt.fm_n = at(inum(GM, 256, [1 4096]), 1, 2);

        % Driving system
        at(uilabel(G, 'Text', 'Driving system (optional)', 'FontWeight', 'bold'), 17, [1 2]);
        dl(G, 'DS serial', 18);              nt.ds_serial = at(uieditfield(G, 'text', 'ValueChangedFcn', upd), 18, 2);
        dl(G, 'DS name', 19);                nt.ds_name   = at(uieditfield(G, 'text'), 19, 2);
        dl(G, 'DS manufacturer', 20);        nt.ds_manu   = at(uieditfield(G, 'text'), 20, 2);
        dl(G, 'DS channels', 21);            nt.ds_ch     = at(inum(G, 4, [1 4096]), 21, 2);
        dl(G, 'Profile CSV', 22);
        pc = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 70}, 'ColumnSpacing', 4), 22, 2);
        nt.csv = uieditfield(pc, 'text', 'ValueChangedFcn', upd);
        nt.csv_btn = uibutton(pc, 'Text', 'Browse…', 'ButtonPushedFcn', @(~,~) browse_new_csv());
        bb = at(uigridlayout(G, [1 2], 'Padding', 0), 23, [1 2]);
        uibutton(bb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) close(d));
        uibutton(bb, 'Text', 'Save transducer', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) save_new());
        nt.msg = at(uitextarea(G, 'Editable', 'off'), 24, [1 2]);
        preview();

        function on_type()
            is_m = strcmp(nt.type.Value, 'matrix');
            A.Visible = ~is_m; M.Visible = is_m;
            nt.csv.Enable = ~is_m; nt.csv_btn.Enable = ~is_m;
            preview();
        end

        function on_grid()
            GR.Visible = strcmp(nt.m_grid.Value, 'rect');
            GF.Visible = strcmp(nt.m_grid.Value, 'fibonacci');
            GM.Visible = strcmp(nt.m_grid.Value, 'fermat');
            preview();
        end

        function preview()
            [t, issues] = collect_new();
            try
                info = draw_transducer_geometry(nt.ax_side, nt.ax_face, t);
            catch ME
                clear_axes(nt.ax_side); clear_axes(nt.ax_face);
                info = {['Preview unavailable: ' ME.message]};
            end
            nt.msg.Value = [info(:); issues(:); {sprintf('Saves {serial}.yaml to %s', st.equip_path)}];
        end

        function fill_from_template()
            if strcmp(nt.tmpl.Value, '(blank)'); return; end
            t  = st.equip.trans.(nt.tmpl.Value);
            tr = t.transducer;
            nt.serial.Value = ''; nt.name.Value = '';
            nt.manu.Value   = char(string(get_field(t, 'manufact', '')));
            nt.freq.Value   = tr.freq_hz / 1e3;
            nt.fmin.Value   = get_field(t, 'min_foc', nt.fmin.Value);
            nt.fmax.Value   = get_field(t, 'max_foc', nt.fmax.Value);
            if is_matrix(t)
                m = tr.matrix;
                nt.type.Value = 'matrix';
                nt.rc.Value  = get_field(m, 'curv_radius_mm', nt.rc.Value);
                nt.dep.Value = get_field(m, 'dist_geom_ep_mm', nt.dep.Value);
                nt.m_steer.Value  = char(get_field(m, 'steering', nt.m_steer.Value));
                nt.m_curved.Value = logical(get_field(m, 'is_curved', true));
                nt.m_D.Value      = m.outer_diameter_mm;
                nt.m_shape.Value  = char(m.elem_shape);
                nt.m_w.Value = m.elem_width_mm; nt.m_h.Value = m.elem_height_mm;
                gs = m.matrix_shape.define_here.grid_shape;
                nt.m_grid.Value = char(gs.type);
                switch nt.m_grid.Value
                    case 'rect'
                        g = gs.rect;
                        nt.r_nrow.Value = g.elem_n_row; nt.r_ncol.Value = g.elem_n_col;
                        nt.r_sh.Value = g.elem_spacing_height_mm; nt.r_sw.Value = g.elem_spacing_width_mm;
                        nt.r_sp.Value = get_field(g, 'sparsity_factor', 1);
                    case 'fibonacci'
                        nt.f_n.Value = gs.fibonacci.elem_n; nt.f_kerf.Value = get_field(gs.fibonacci, 'kerf_mm', 0);
                    case 'fermat'
                        nt.fm_n.Value = gs.fermat.elem_n;
                end
                on_grid();
            else
                a = tr.annular;
                nt.type.Value = 'annular';
                nt.rc.Value   = a.curv_radius_mm;
                nt.dep.Value  = a.dist_geom_ep_mm;
                nt.nch.Value  = get_field(t, 'n_elem', a.elem_n);
                nt.id.Value   = strjoin(compose('%g', a.elem_id_mm(:)'), ', ');
                nt.od.Value   = strjoin(compose('%g', a.elem_od_mm(:)'), ', ');
            end
            on_type();
            nt.msg.Value = [{sprintf('Copied geometry from %s — enter a new serial.', nt.tmpl.Value)}; nt.msg.Value(:)];
        end

        function browse_new_csv()
            [f, p] = uigetfile('*.csv', 'Axial profile CSV', st.prof_root);
            figure(d);
            if isequal(f, 0); return; end
            nt.csv.Value = fullfile(p, f);
            preview();
        end

        function save_new()
            [t, issues] = collect_new();
            if ~isempty(issues)
                nt.msg.Value = [{'Not saved:'}; issues(:)]; return;
            end
            % CHECK_EQUIPMENT_CONFIG also flags some shipped entries (touching
            % sub-rings, max_foc beyond the curvature radius), so warn only
            warns = validate_equipment_entry(t);
            if ~isempty(warns)
                nt.msg.Value = [{'check_equipment_config:'}; warns(:)];
                sel = uiconfirm(d, strjoin(warns, newline), 'Equipment check', ...
                                'Options', {'Save anyway', 'Cancel'}, 'DefaultOption', 2, 'Icon', 'warning');
                if ~strcmp(sel, 'Save anyway'); return; end
            end
            f = fullfile(st.equip_path, [t.serial '.yaml']);
            if isfile(f)
                sel = uiconfirm(d, sprintf('%s already exists. Overwrite it?', f), 'Overwrite', ...
                                'Options', {'Overwrite', 'Cancel'}, 'DefaultOption', 2);
                if ~strcmp(sel, 'Overwrite'); return; end
            end
            yaml.dumpFile(f, t);
            st.equip = load_equipment_config(st.equip_path);
            ui.tran.Items = transducer_serials(st.equip);
            ui.tran.Value = t.serial;
            close(d);
            if is_matrix(t)
                say(sprintf('Added matrix transducer %s (%s).', t.serial, f));
            else
                say(sprintf('Added transducer %s (%s). No library entry yet: showing the geometric prediction.', ...
                    t.serial, f));
            end
            on_transducer();
        end

        function [t, issues] = collect_new()
            issues = {};
            serial = strtrim(nt.serial.Value);
            if ~isvarname(serial)
                issues{end+1} = 'Serial must start with a letter and contain only letters, digits and underscores.';
            end
            if nt.fmin.Value >= nt.fmax.Value
                issues{end+1} = 'Focal range: minimum must be below maximum.';
            end
            name = strtrim(nt.name.Value);
            if isempty(name); name = serial; end
            t = struct('type', 'transducer', 'serial', serial, 'name', name, ...
                       'manufact', strtrim(nt.manu.Value), 'n_elem', int32(1), ...
                       'min_foc', nt.fmin.Value, 'max_foc', nt.fmax.Value);
            if isempty(t.manufact); t = rmfield(t, 'manufact'); end

            if strcmp(nt.type.Value, 'annular')
                [id, ok_id] = parse_list(nt.id.Value);
                [od, ok_od] = parse_list(nt.od.Value);
                if ~ok_id || isempty(id); issues{end+1} = 'Ring inner Ø: enter one number per ring.'; end
                if ~ok_od || isempty(od); issues{end+1} = 'Ring outer Ø: enter one number per ring.'; end
                if ok_id && ok_od && numel(id) ~= numel(od)
                    issues{end+1} = sprintf('%d inner but %d outer diameters.', numel(id), numel(od));
                elseif any(od <= id)
                    issues{end+1} = 'Each outer diameter must exceed its inner diameter.';
                end
                if mod(numel(id), nt.nch.Value) ~= 0
                    issues{end+1} = sprintf(['Ring count (%d) must be a multiple of the physical channel count (%d); ' ...
                        'rings are grouped into channels in order.'], numel(id), nt.nch.Value);
                end
                t.n_elem = int32(nt.nch.Value);
                t.transducer.annular = struct('elem_n', int32(numel(id)), 'elem_id_mm', id, 'elem_od_mm', od, ...
                                              'curv_radius_mm', nt.rc.Value, 'dist_geom_ep_mm', nt.dep.Value);
            else
                gs.type = nt.m_grid.Value;
                switch gs.type
                    case 'rect'
                        gs.rect = struct('elem_n_row', int32(nt.r_nrow.Value), 'elem_n_col', int32(nt.r_ncol.Value), ...
                                         'elem_spacing_height_mm', nt.r_sh.Value, ...
                                         'elem_spacing_width_mm', nt.r_sw.Value, 'sparsity_factor', nt.r_sp.Value);
                    case 'fibonacci'
                        gs.fibonacci = struct('elem_n', int32(nt.f_n.Value), 'kerf_mm', nt.f_kerf.Value);
                    case 'fermat'
                        gs.fermat = struct('elem_n', int32(nt.fm_n.Value));
                end
                m = struct('steering', nt.m_steer.Value, 'elem_shape', nt.m_shape.Value, ...
                           'elem_width_mm', nt.m_w.Value, 'elem_height_mm', nt.m_h.Value, ...
                           'outer_diameter_mm', nt.m_D.Value, 'is_curved', nt.m_curved.Value, ...
                           'curv_radius_mm', nt.rc.Value, 'dist_geom_ep_mm', nt.dep.Value);
                m.matrix_shape = struct('type', 'define_here', 'define_here', struct('grid_shape', gs));
                if m.is_curved && m.curv_radius_mm <= m.outer_diameter_mm / 2
                    issues{end+1} = 'Curvature radius must exceed half the aperture diameter.';
                else
                    try
                        t.n_elem = int32(size(matrix_elem_pos_mm(m), 1));
                    catch ME
                        issues{end+1} = ME.message;
                    end
                end
                t.transducer.matrix = m;
            end
            t.transducer.freq_hz = nt.freq.Value * 1e3;

            ds_serial = strtrim(nt.ds_serial.Value);
            if ~isempty(ds_serial) && ~isvarname(ds_serial)
                issues{end+1} = 'DS serial must start with a letter and contain only letters, digits and underscores.';
            end
            csv = strtrim(nt.csv.Value);
            if ~isempty(ds_serial)
                ds_name = strtrim(nt.ds_name.Value);
                if isempty(ds_name); ds_name = ds_serial; end
                c = struct();
                c.ds = struct('serial', ds_serial, 'name', ds_name, ...
                              'manufact', strtrim(nt.ds_manu.Value), 'available_ch', int32(nt.ds_ch.Value));
                if isempty(c.ds.manufact); c.ds = rmfield(c.ds, 'manufact'); end
                if ~isempty(csv) && strcmp(nt.type.Value, 'annular')
                    % Store relative to ProfileRoot when possible, as the shipped entries do
                    if ~isempty(st.prof_root) && startsWith(canonical(csv), [canonical(st.prof_root) filesep])
                        csv = extractAfter(canonical(csv), [canonical(st.prof_root) filesep]);
                    end
                    c.char_data_path = csv;
                end
                t.combos.(ds_serial) = c;
            elseif ~isempty(csv) && strcmp(nt.type.Value, 'annular')
                issues{end+1} = 'A profile CSV belongs to a driving system: enter a DS serial.';
            end
        end
    end

    %% Drawing ----------------------------------------------------------------
    function draw()
        draw_profile(); draw_r2(); draw_params();
    end

    function draw_profile()
        ax = ui.ax_prof; clear_axes(ax); hold(ax, 'on');
        if isempty(st.prof) || ~isfield(st.res, 'geometric'); return; end
        show_all = strcmp(ui.show.Value, 'All depths');
        d = str2double(ui.show.Value);
        styles = result_styles();
        shown_meas = false;
        for k = 1:size(styles, 1)
            if ~isfield(st.res, styles{k, 1}); continue; end
            r = st.res.(styles{k, 1});
            if show_all; j = 1:numel(r.depths); else; j = find(abs(r.depths - d) < 0.01, 1); end
            if isempty(j); continue; end
            if ~shown_meas && ~isempty(r.Imeas)
                % Measured profiles solid; interpolated/extrapolated ones dashed
                jm = j(any(isfinite(r.Imeas(:, j)), 1));
                pt = arrayfun(@(x) profile_type_at(x), r.depths(jm), 'UniformOutput', false);
                for t = {'measured', '-', 'Hydrophone'; 'interpolated', '--', 'Interpolated profile'; ...
                         'extrapolated', ':', 'Extrapolated profile'}'
                    jt = jm(strcmp(pt, t{1}));
                    if isempty(jt); continue; end
                    h = plot(ax, r.z, r.Imeas(:, jt), t{2}, 'Color', 'k', 'LineWidth', 1.8, 'DisplayName', t{3});
                    set(h(2:end), 'HandleVisibility', 'off');
                end
                shown_meas = true;
            end
            name = styles{k, 3};
            if ~show_all && isfinite(r.r2(j)); name = sprintf('%s (R² %.2f)', name, r.r2(j)); end
            if ~show_all && ~isempty(r.settings); name = sprintf('%s [%s]', name, phase_label(r.settings(j))); end
            h = plot(ax, r.z, r.Imod(:, j), '-', 'Color', styles{k, 2}, 'LineWidth', 1.2, 'DisplayName', name);
            set(h(2:end), 'HandleVisibility', 'off');
        end
        if show_all
            xline(ax, st.prof.depths_ep_mm, ':', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            where = 'all depths';
        else
            xline(ax, d, ':', 'HandleVisibility', 'off');
            where = sprintf('%.1f mm', d);
        end
        xlabel(ax, 'Distance from exit plane (mm)'); ylabel(ax, 'Normalised intensity');
        title(ax, sprintf('%s — %s%s', strrep(combo_name(), '_', '\_'), where, ...
            ternary(st.predicted, ' (prediction, runtime library path)', '')));
        legend(ax, 'Location', 'northeast'); grid(ax, 'on'); box(ax, 'on');
    end

    function draw_r2()
        ax = ui.ax_r2; clear_axes(ax); hold(ax, 'on');
        styles = result_styles();
        for k = 1:size(styles, 1)
            if ~isfield(st.res, styles{k, 1}); continue; end
            r = st.res.(styles{k, 1});
            y = ternary(st.predicted, r.focus, r.r2);
            plot(ax, r.depths, y, '-o', 'Color', styles{k, 2}, 'MarkerFaceColor', styles{k, 2}, ...
                'LineWidth', 1.2, 'DisplayName', styles{k, 3});
        end
        tran = st.equip.trans.(ui.tran.Value);
        if ~is_matrix(tran)
            xline(ax, tran.transducer.annular.dist_geom_ep_mm, ':', 'natural focus', 'HandleVisibility', 'off');
        end
        if st.predicted
            lim = [min(st.prof.depths_ep_mm), max(st.prof.depths_ep_mm)];
            plot(ax, lim, lim, 'k--', 'DisplayName', 'Target');
            if isfield(st.res, 'geometric')
                % Focus of the measured (filled) and derived (open) profiles
                r  = st.res.geometric;
                jm = find(any(isfinite(r.Imeas), 1));
                pt = st.prof.profile_type(jm);
                for t = {'measured', 'k', 'Hydrophone'; 'interpolated', 'none', 'Interpolated profile'; ...
                         'extrapolated', 'none', 'Extrapolated profile'}'
                    jt = jm(strcmp(pt, t{1}));
                    if isempty(jt); continue; end
                    fm = calibration_focus_beyond(r.Imeas(:, jt), r.z, ui.skip.Value);
                    plot(ax, r.depths(jt), fm, 'kp', 'MarkerFaceColor', t{2}, 'MarkerSize', 9, 'DisplayName', t{3});
                end
            end
            if ~isempty(st.meas)
                % Shade the measured span: outside it the entry is extrapolated
                yl = ylim(ax); mr = [min(st.meas.depths_ep_mm), max(st.meas.depths_ep_mm)];
                h = patch(ax, mr([1 2 2 1]), yl([1 1 2 2]), [0.5 0.5 0.5], 'FaceAlpha', 0.08, ...
                    'EdgeColor', 'none', 'DisplayName', 'Measured range');
                uistack(h, 'bottom'); ylim(ax, yl);
            end
            xlabel(ax, 'Target depth (mm from EP)'); ylabel(ax, 'Focus, FLHM centre (mm)');
            legend(ax, 'Location', 'northwest');
        else
            ylim(ax, [max(-0.2, min([ylim(ax) 0])) 1]);
            xlabel(ax, 'Focal depth (mm from EP)'); ylabel(ax, 'R² (profile shape)');
            legend(ax, 'Location', 'southwest');
        end
        grid(ax, 'on'); box(ax, 'on');
    end

    function draw_params()
        if st.predicted
            draw_settings(); return;
        end
        if isempty(st.fit); ui.params.Data = {}; return; end
        m = st.fit.model;
        ph = cell2mat(cellfun(@(x) x.elem_phase_correction_deg(:), m, 'UniformOutput', false));
        ui.params.ColumnName = [{'Element'}, arrayfun(@(x) sprintf('%.1f mm', x.fit_depths_ep_mm), ...
            [m{:}], 'UniformOutput', false)];
        ui.params.Data = num2cell([(1:size(ph, 1))', round(ph, 1)]);
    end

    function draw_settings()
        % Runtime settings per desired depth of the fit (else the library entry)
        [src, label] = settings_source();
        if isempty(src)
            ui.params.ColumnName = {'Depth (mm)'}; ui.params.Data = {};
            return;
        end
        r = st.res.(src); s = r.settings;
        n_elem = numel(s(1).elem_phase_deg);
        % Drive: mean source pressure at the reference ISPPA (before per-ring gains)
        ui.params.ColumnName = [{'Depth (mm)', 'Profile', 'Entry', 'Phase model', 'Drive (kPa)', ...
            'Focus (mm)', 'R²'}, arrayfun(@(i) sprintf('Ring %d (°)', i), 1:n_elem, 'UniformOutput', false)];
        rows = cell(numel(s), 7 + n_elem);
        for k = 1:numel(s)
            drive = '—';
            if ~isempty(s(k).drive_pa); drive = sprintf('%.1f', s(k).drive_pa / 1e3); end
            r2 = '—';
            if isfinite(r.r2(k)); r2 = sprintf('%.3f', r.r2(k)); end
            pt = strrep(profile_type_at(r.depths(k)), 'unavailable', '—');
            rows(k, :) = [{sprintf('%g', r.depths(k)), pt, label, phase_label(s(k)), drive, ...
                sprintf('%.1f', r.focus(k)), r2}, compose('%.1f', s(k).elem_phase_deg)];
        end
        ui.params.Data = rows;
    end

    function t = profile_type_at(d)
        % Profile type of the depth on view ('measured' in the measured view)
        t = 'measured';
        if isfield(st.prof, 'profile_type')
            k = find(abs(st.prof.depths_ep_mm - d) < 0.01, 1);
            if ~isempty(k); t = st.prof.profile_type{k}; end
        end
    end

    function [src, label] = settings_source()
        src = ''; label = '';
        if isfield(st.res, 'fit') && ~isempty(st.res.fit.settings)
            src = 'fit'; label = 'Fit (unsaved)';
        elseif isfield(st.res, 'library') && ~isempty(st.res.library.settings)
            src = 'library'; label = 'Library';
        end
    end

    %% Helpers ----------------------------------------------------------------
    function o = base_options()
        o = struct('ForwardModel', ui.fwd.Value, 'NearFieldSkipMm', ui.skip.Value, ...
                   'HuberDelta', ui.huber.Value, 'NStarts', ui.nstarts.Value, 'Domain', ui.domain.Value);
        for f = fieldnames(st.adv)'; o.(f{1}) = st.adv.(f{1}); end
    end

    function on_advanced()
        % Advanced options: FIT_TRANSDUCER_MODEL (st.adv), the desired-profile
        % tolerance (st.pext) and the k-Wave drive calibration (st.kwo). Only
        % values that differ from the defaults are stored, so "Reset" restores
        % the default fit. "Start from preset" loads a CALIBRATION_FIT_OPTIONS
        % preset into this dialog and the main panel.
        D = advanced_defaults(); E = extra_defaults();
        cur = D; for fa = fieldnames(st.adv)'; cur.(fa{1}) = st.adv.(fa{1}); end
        ext = E; for fa = fieldnames(st.pext)'; ext.(fa{1}) = st.pext.(fa{1}); end
        for fa = fieldnames(st.kwo)'; ext.(fa{1}) = st.kwo.(fa{1}); end
        n_rows = 30;
        d = uifigure('Name', 'Advanced fit options', 'Position', [200 120 460 760], 'WindowStyle', 'modal');
        G = uigridlayout(d, [n_rows 2], 'RowHeight', repmat({22}, 1, n_rows), 'ColumnWidth', {210, '1x'}, ...
            'RowSpacing', 3, 'Scrollable', 'on');
        al  = @(txt, r, tip) at(uilabel(G, 'Text', txt, 'Tooltip', tip), r, 1);
        num = @(r) at(uieditfield(G, 'numeric'), r, 2);
        hdr = @(txt, r) at(uilabel(G, 'Text', txt, 'FontWeight', 'bold'), r, [1 2]);
        function [a, b] = pair_row(r, kind)
            g = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), r, 2);
            a = uieditfield(g, kind); b = uieditfield(g, kind);
        end
        al('Start from preset', 1, ['CALIBRATION_FIT_OPTIONS preset: sets the main panel (free parameters, ' ...
            'forward model, loss, near-field skip) and these options to the preset''s values']);
        w.pre = at(uidropdown(G, 'Items', {'(current settings)', 'rayleigh_phase_gain', 'oneil_phase_only'}, ...
            'ValueChangedFcn', @(src, ~) load_preset(src.Value)), 1, 2);
        hdr('Phase model', 2);
        al('Phase constraint', 3, 'none: free offsets per ring; linear / monotonic: PERFORM_GLOBAL_SEARCH constraints');
        w.pc = at(uidropdown(G, 'Items', {'none', 'linear', 'monotonic'}), 3, 2);
        al('Step bounds (rad; blank = default)', 4, 'PhaseStepBounds');
        [w.sb1, w.sb2] = pair_row(4, 'text');
        al('Phase prior λ', 5, 'Lambda: pull of the phase offsets towards geometric steering');
        w.lam = num(5);
        al('Gain prior λ', 6, 'LambdaAmp: pull of the log ring gains towards equal gains');
        w.lama = num(6);
        hdr('Gains', 7);
        al('Gain bounds (× prior)', 8, 'AmpBounds: bounds on each ring gain relative to the prior');
        [w.ab1, w.ab2] = pair_row(8, 'numeric');
        al('Plausible gain (× mean; flag only)', 9, 'PlausibleGain: gains outside are flagged as a possible hardware fault');
        [w.pg1, w.pg2] = pair_row(9, 'numeric');
        hdr('Amplitude and loss', 10);
        al('Amplitude fit', 11, 'shape: closed-form scale per depth; absolute: drive velocity fitted to the measured W/cm²');
        w.af = at(uidropdown(G, 'Items', {'shape', 'absolute'}), 11, 2);
        w.sv = at(uicheckbox(G, 'Text', 'one drive for all depths (absolute)'), 12, [1 2]);
        al('Focal weight (0 = uniform, ≥1)', 13, 'FocalWeight: Gaussian weights around the measured focus');
        w.fw = num(13);
        al('Fit window (mm EP; blank = all)', 14, 'FitWindowMm');
        [w.fw1, w.fw2] = pair_row(14, 'text');
        hdr('Near field', 15);
        al('Near-field penalty (0 = off)', 16, ['NearFieldPenalty: penalises model near-field intensity above the ' ...
            'measured near field (between the minimum and the near-field skip)']);
        w.nfp = num(16);
        al('Margin (× focal peak)', 17, 'NearFieldMargin: excess tolerated before the penalty applies');
        w.nfm = num(17);
        al('From (mm EP)', 18, 'NearFieldMinMm: start of the penalised near field');
        w.nfz = num(18);
        hdr('Medium and optimiser', 19);
        al('Sound speed (m/s) / density (kg/m³)', 20, 'CWater, Rho: free-water medium of the forward model');
        [w.cw, w.rho] = pair_row(20, 'numeric');
        al('Optimiser', 21, 'lsqnonlin (NStarts starts) | fexminimize | globalsearch');
        g5 = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4, 'ColumnWidth', {'1x', 70}), 21, 2);
        w.opt = uidropdown(g5, 'Items', {'lsqnonlin', 'fexminimize', 'globalsearch'});
        w.pol = uicheckbox(g5, 'Text', 'polish');
        al('FEX population / max evals', 22, 'PopSize, MaxFunEvals (fexminimize)');
        [w.pop, w.mfe] = pair_row(22, 'numeric');
        w.pop.Limits = [10 1e5]; w.pop.RoundFractionalValues = 'on';
        w.mfe.Limits = [100 1e8]; w.mfe.RoundFractionalValues = 'on';
        al('Random seed', 23, 'Seed of the extra starts');
        w.seed = num(23); w.seed.RoundFractionalValues = 'on';
        hdr('Desired profiles', 24);
        al('Same-depth tolerance (mm)', 25, ['Tolerance: a desired depth this close to a measured one uses the ' ...
            'measurement instead of an interpolated profile']);
        w.tol = num(25);
        hdr('k-Wave drive calibration', 26);
        al('Grid (blank = sim config)', 27, ['GridDims: ''auto'' (sized to the aperture and focus) or ' ...
            '[nx ny nz]; blank = the sim config''s grid']);
        w.grid = at(uieditfield(G, 'text'), 27, 2);
        w.apply = at(uicheckbox(G, 'Text', 'rescale the stored drives to the simulated ISPPA', ...
            'Tooltip', 'ApplyCorrection; off: report only'), 28, [1 2]);
        bb = at(uigridlayout(G, [1 3], 'Padding', 0), 30, [1 2]);
        uibutton(bb, 'Text', 'Reset', 'ButtonPushedFcn', @(~,~) apply(struct(), struct(), struct()));
        uibutton(bb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) close(d));
        uibutton(bb, 'Text', 'OK', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) ok());
        fill(cur, ext);

        function fill(c, e)
            w.pc.Value = char(c.PhaseConstraint);
            sb = c.PhaseStepBounds; if isempty(sb); sb = [NaN NaN]; end
            w.sb1.Value = num2str_blank(sb(1)); w.sb2.Value = num2str_blank(sb(2));
            w.lam.Value = c.Lambda; w.lama.Value = c.LambdaAmp;
            w.ab1.Value = c.AmpBounds(1); w.ab2.Value = c.AmpBounds(2);
            w.pg1.Value = c.PlausibleGain(1); w.pg2.Value = c.PlausibleGain(2);
            w.af.Value = char(c.AmplitudeFit); w.sv.Value = logical(c.SharedVelocity);
            w.fw.Value = c.FocalWeight;
            fw = c.FitWindowMm; if isempty(fw); fw = [NaN NaN]; end
            w.fw1.Value = num2str_blank(fw(1)); w.fw2.Value = num2str_blank(fw(2));
            w.nfp.Value = c.NearFieldPenalty; w.nfm.Value = c.NearFieldMargin; w.nfz.Value = c.NearFieldMinMm;
            w.cw.Value = c.CWater; w.rho.Value = c.Rho;
            w.opt.Value = char(c.Optimizer); w.pol.Value = logical(c.Polish);
            w.pop.Value = c.PopSize; w.mfe.Value = c.MaxFunEvals; w.seed.Value = c.Seed;
            w.tol.Value = e.Tolerance;
            g = e.GridDims; if isnumeric(g); g = strjoin(compose('%g', g), ' '); end
            w.grid.Value = char(g); w.apply.Value = logical(e.ApplyCorrection);
        end
        function load_preset(name)
            if strcmp(name, '(current settings)'); return; end
            p = calibration_fit_options(name);
            c = D;
            for fp_ = fieldnames(p)'; if isfield(c, fp_{1}); c.(fp_{1}) = p.(fp_{1}); end; end
            fill(c, ext);
            free = cellstr(p.FreeParams);
            ui.fp_phase.Value = ismember('phase', free); ui.fp_amp.Value = ismember('amp', free);
            ui.fwd.Value = p.ForwardModel; ui.domain.Value = p.Domain;
            ui.huber.Value = p.HuberDelta; ui.skip.Value = p.NearFieldSkipMm;
            say(sprintf('Preset %s loaded into the main panel and the advanced options (press OK to keep).', name));
        end
        function ok()
            v = struct('PhaseConstraint', w.pc.Value, 'PhaseStepBounds', pair_or_empty(w.sb1, w.sb2), ...
                'AmplitudeFit', w.af.Value, 'SharedVelocity', w.sv.Value, 'Optimizer', w.opt.Value, ...
                'PopSize', w.pop.Value, 'MaxFunEvals', w.mfe.Value, 'Polish', w.pol.Value, ...
                'FocalWeight', w.fw.Value, 'FitWindowMm', pair_or_empty(w.fw1, w.fw2), ...
                'Lambda', w.lam.Value, 'LambdaAmp', w.lama.Value, ...
                'AmpBounds', sort([w.ab1.Value w.ab2.Value]), 'PlausibleGain', sort([w.pg1.Value w.pg2.Value]), ...
                'NearFieldPenalty', w.nfp.Value, 'NearFieldMargin', w.nfm.Value, 'NearFieldMinMm', w.nfz.Value, ...
                'CWater', w.cw.Value, 'Rho', w.rho.Value, 'Seed', w.seed.Value);
            g = strtrim(w.grid.Value);
            gd = str2num(g); %#ok<ST2NM> numeric dims, else 'auto' or ''
            if ~isempty(gd); g = gd; end
            if ~(isempty(g) || isequal(g, 'auto') || (isnumeric(g) && any(numel(g) == [2 3])))
                uialert(d, 'Grid: enter auto, [nx nz] or [nx ny nz], or leave it blank.', 'Advanced'); return;
            end
            apply(non_default(v, D), non_default(struct('Tolerance', w.tol.Value), E), ...
                  non_default(struct('GridDims', g, 'ApplyCorrection', w.apply.Value), E));
        end
        function apply(v, pe, kw)
            st.adv = v; st.pext = pe; st.kwo = kw;
            update_advanced_button();
            names = [fieldnames(v); fieldnames(pe); fieldnames(kw)];
            if ~isempty(names)
                say(sprintf('Advanced options: %s', strjoin(names, ', ')));
            else
                say('Advanced options reset to defaults.');
            end
            close(d);
        end
    end

    function update_advanced_button()
        n = numel(fieldnames(st.adv)) + numel(fieldnames(st.pext)) + numel(fieldnames(st.kwo));
        ui.adv_btn.Text = ternary(n > 0, sprintf('Advanced (%d)', n), 'Advanced…');
        n = numel(fieldnames(st.rpo));
        ui.rpo_btn.Text = ternary(n > 0, sprintf('Options (%d)', n), 'Options…');
    end

    function on_ring_options()
        % FIT_RING_PROFILE options (st.rpo, non-defaults only), incl. a
        % plane-scan file (READ_OFFAXIS_SCANS) fitted with the axial profiles
        R0 = ring_defaults();
        cur = R0; for fr = fieldnames(st.rpo)'; cur.(fr{1}) = st.rpo.(fr{1}); end
        n_rows = 19;
        d = uifigure('Name', 'Ring profile options', 'Position', [220 160 470 560], 'WindowStyle', 'modal');
        G = uigridlayout(d, [n_rows 2], 'RowHeight', repmat({22}, 1, n_rows), 'ColumnWidth', {200, '1x'}, ...
            'RowSpacing', 3, 'Scrollable', 'on');
        al  = @(txt, r, tip) at(uilabel(G, 'Text', txt, 'Tooltip', tip), r, 1);
        num = @(v, r) at(uieditfield(G, 'numeric', 'Value', v), r, 2);
        hdr = @(txt, r) at(uilabel(G, 'Text', txt, 'FontWeight', 'bold'), r, [1 2]);
        hdr('Profile model (used when Ring profile = fit)', 1);
        al('Sub-rings per ring', 2, 'NSub');                       w.nsub = num(cur.NSub, 2);
        w.nsub.Limits = [1 32]; w.nsub.RoundFractionalValues = 'on';
        al('Smoothness (0 = none)', 3, 'Smooth: penalty between neighbouring sub-rings'); w.sm = num(cur.Smooth, 3);
        al('Near-field skip (mm EP)', 4, 'NearFieldSkipMm: axial samples fitted from here'); w.skip = num(cur.NearFieldSkipMm, 4);
        al('Weight bounds', 5, 'AmpBounds: bounds on sub-ring weights and channel gains');
        g = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 5, 2);
        w.ab1 = uieditfield(g, 'numeric', 'Value', cur.AmpBounds(1));
        w.ab2 = uieditfield(g, 'numeric', 'Value', cur.AmpBounds(2));
        hdr('Optimiser', 6);
        al('Alternations', 7, 'Iterations: drives ↔ sub-ring weights');  w.it = num(cur.Iterations, 7);
        al('Starts per depth', 8, 'NStarts');                            w.ns = num(cur.NStarts, 8);
        al('Random seed', 9, 'Seed');                                    w.seed = num(cur.Seed, 9);
        w.ho = at(uicheckbox(G, 'Text', 'leave-one-depth-out check (slower)', 'Value', logical(cur.HeldOut), ...
            'Tooltip', 'HeldOut: also score each depth with the profile fitted on the other depths'), 10, [1 2]);
        hdr('Plane / volume scans (off-axis data)', 11);
        al('File (blank = axial only)', 12, ['OffAxis: Sonic Concepts FUS Metrology export (.mat, folder or ' ...
            '.zip) or CSV focus_ep_mm, x_mm, y_mm, z_ep_mm, intensity_wcm2 (READ_OFFAXIS_SCANS)']);
        fb = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4, 'ColumnWidth', {'1x', 60}), 12, 2);
        w.file = uieditfield(fb, 'text', 'Value', char(cur.OffAxis));
        uibutton(fb, 'Text', 'Browse…', 'ButtonPushedFcn', @(~,~) browse());
        al('Scans used', 13, 'Kind: all | transverse (XY planes) | axial (XZ / YZ planes) | volume (3D scans)');
        w.kind = at(uidropdown(G, 'Items', {'all', 'transverse', 'axial', 'volume'}, 'Value', char(cur.Kind)), 13, 2);
        al('Weight per plane (× axial profile)', 14, 'Weight'); w.wt = num(cur.Weight, 14);
        al('Radial range / bin (mm)', 15, 'RMaxMm, BinMm');
        g = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 15, 2);
        w.rmax = uieditfield(g, 'numeric', 'Value', cur.RMaxMm);
        w.bin  = uieditfield(g, 'numeric', 'Value', cur.BinMm);
        al('Axial-plane / volume z bin (mm)', 16, 'ZBinMm: coarser = faster'); w.zbin = num(cur.ZBinMm, 16);
        al('Transverse / volume centre', 17, 'Centre: origin (scan axis at x = y = 0) | centroid of the top 50 %');
        w.cen = at(uidropdown(G, 'Items', {'origin', 'centroid'}, 'Value', char(cur.Centre)), 17, 2);
        bb = at(uigridlayout(G, [1 3], 'Padding', 0), 19, [1 2]);
        uibutton(bb, 'Text', 'Reset', 'ButtonPushedFcn', @(~,~) apply(struct()));
        uibutton(bb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) close(d));
        uibutton(bb, 'Text', 'OK', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) ok());

        function browse()
            [f, p] = uigetfile({'*.mat;*.zip;*.csv', 'Plane scans (*.mat, *.zip, *.csv)'}, 'Plane scans');
            figure(d);
            if ~isequal(f, 0); w.file.Value = fullfile(p, f); end
        end
        function ok()
            file = strtrim(w.file.Value);
            if ~isempty(file) && ~(isfile(file) || isfolder(file))
                uialert(d, sprintf('Not found: %s', file), 'Ring profile options'); return;
            end
            v = struct('NSub', w.nsub.Value, 'Smooth', w.sm.Value, 'NearFieldSkipMm', w.skip.Value, ...
                'AmpBounds', sort([w.ab1.Value w.ab2.Value]), 'Iterations', w.it.Value, 'NStarts', w.ns.Value, ...
                'Seed', w.seed.Value, 'HeldOut', w.ho.Value, 'OffAxis', file, 'Kind', w.kind.Value, ...
                'Weight', w.wt.Value, 'RMaxMm', w.rmax.Value, 'BinMm', w.bin.Value, 'ZBinMm', w.zbin.Value, ...
                'Centre', w.cen.Value);
            apply(non_default(v, R0));
        end
        function apply(v)
            st.rpo = v;
            update_advanced_button();
            if isfield(v, 'OffAxis') && ~strcmp(ui.ring.Value, 'fit')
                say('Plane scans are used when Ring profile is set to "fit".');
            end
            close(d);
        end
    end

    function ro = ring_options()
        % RingProfileOptions for CALIBRATION_FIT from st.rpo
        R0 = ring_defaults(); v = R0;
        for f = fieldnames(st.rpo)'; v.(f{1}) = st.rpo.(f{1}); end
        ro = struct('NSub', v.NSub, 'Smooth', v.Smooth, 'NearFieldSkipMm', v.NearFieldSkipMm, ...
            'AmpBounds', v.AmpBounds, 'Iterations', v.Iterations, 'NStarts', v.NStarts, 'Seed', v.Seed, ...
            'HeldOut', logical(v.HeldOut));
        if ~isempty(v.OffAxis)
            ro.OffAxis = char(v.OffAxis);
            ro.OffAxisOptions = struct('Kind', v.Kind, 'Weight', v.Weight, 'RMaxMm', v.RMaxMm, ...
                'BinMm', v.BinMm, 'ZBinMm', v.ZBinMm, 'Centre', v.Centre);
        end
    end

    function apply_fit_record(lib)
        % Start from the settings this entry was fitted with (meta.fit, written
        % by CALIBRATION_SAVE_ENTRY); entries without it keep the current settings
        rec = get_field(get_field(lib, 'meta', []), 'fit', []);
        if ~isstruct(rec); return; end
        fo = get_field(rec, 'fit_options', struct());
        ui.ring.Value = 'auto';   % refits keep the entry's ring profile unless changed
        free = cellstr(get_field(fo, 'FreeParams', {'phase', 'amp'}));
        ui.fp_phase.Value = ismember('phase', free);
        ui.fp_amp.Value   = ismember('amp', free);
        for c = {'ForwardModel', 'fwd'; 'Domain', 'domain'}'
            v = char(get_field(fo, c{1}, ''));
            if any(strcmp(ui.(c{2}).Items, v)); ui.(c{2}).Value = v; end
        end
        for c = {'NearFieldSkipMm', 'skip'; 'HuberDelta', 'huber'; 'NStarts', 'nstarts'}'
            if isfield(fo, c{1}); ui.(c{2}).Value = double(fo.(c{1})); end
        end
        ui.isppa.Value = double(get_field(rec, 'isppa_wcm2', ui.isppa.Value));
        po = get_field(rec, 'profile_options', struct());
        if any(strcmp(ui.align.Items, char(get_field(po, 'Align', '')))); ui.align.Value = char(po.Align); end
        if isfield(po, 'MaxExtrapolationMm'); ui.extrap.Value = double(po.MaxExtrapolationMm); end
        desired = strcmp(char(get_field(rec, 'calibrate_at', 'measured')), 'desired');
        ui.calat.Value = ui.calat.Items{1 + desired};
        if desired
            ui.pred.Value = strjoin(compose('%g', double(rec.depths_ep_mm(:)')), ', ');
        end
        D = advanced_defaults(); st.adv = struct();
        for f = fieldnames(D)'
            if isfield(fo, f{1}) && ~isequal(double_if_numeric(fo.(f{1})), D.(f{1}))
                st.adv.(f{1}) = fo.(f{1});
            end
        end
        st.pext = non_default(struct('Tolerance', double(get_field(po, 'Tolerance', 0.01))), extra_defaults());
        st.rpo = struct();
        rr = get_field(get_field(rec, 'ring_profile', struct()), 'options', struct());
        if isstruct(rr)
            oo = get_field(rr, 'OffAxisOptions', struct());
            for f = fieldnames(oo)'; rr.(f{1}) = oo.(f{1}); end
            R0 = ring_defaults(); v = struct();
            for f = fieldnames(R0)'
                if isfield(rr, f{1}) && ~strcmp(char(string(rr.(f{1}))), '(struct)'); v.(f{1}) = double_if_numeric(rr.(f{1})); end
            end
            st.rpo = non_default(v, R0);
        end
        update_advanced_button();
        say(sprintf('Fit settings taken from the entry (meta.fit: %s, %s, %s).', char(get_field(rec, 'mode', '')), ...
            char(get_field(rec, 'tool', '')), char(string(get_field(lib.meta, 'last_updated', '')))));
    end

    function p = selected_profiles()
        % Selected measured depths, in the profile form of INTERPOLATE_AXIAL_PROFILES
        d = str2double(ui.depths.Value);
        sel = arrayfun(@(x) any(abs(d - x) < 0.01), st.meas.depths_ep_mm);
        p = interpolate_axial_profiles(st.meas, st.meas.depths_ep_mm(sel));
    end

    function lib = load_entry()
        lib = [];
        for folder = {ui.out.Value, st.lib_path}
            f = fullfile(folder{1}, [combo_name() '.yaml']);
            if ~isempty(folder{1}) && isfile(f)
                lib = yaml.loadFile(f, 'ConvertToArray', true); return;
            end
        end
    end

    function n = combo_name()
        n = ui.tran.Value;
        if ~any(strcmp(ui.ds.Value, {'', '(none)'})); n = [n '_' ui.ds.Value]; end
    end

    function say(msg)
        ui.status.Value = [{sprintf('[%s] %s', datestr(now, 'HH:MM:SS'), msg)}; ui.status.Value(:)]; %#ok<TNOW1,DATST>
        drawnow limitrate;
    end
end

%% ==========================================================================
function s = result_styles()
% Result field, colour, legend name
    s = {'geometric', [0.6 0.6 0.6], 'Geometric'; 'library', [0.0 0.45 0.70], 'Library'; ...
         'fit', [0.85 0.37 0.01], 'Fit'; 'kwave', [0.20 0.63 0.17], 'k-Wave'};
end

function txt = phase_label(s)
% Runtime phase model of one LIBRARY_MODEL_AT_DEPTH result
    txt = s.phase_model;
    if strcmp(txt, 'blend'); txt = sprintf('blend w=%.2f', s.refinement_weight); end
end

function D = advanced_defaults()
% Advanced FIT_TRANSDUCER_MODEL options and their fitter defaults ("Advanced…").
% Not offered: the z shift and steering sound speed (not applied at run time),
% ChannelGroups (equipment), Prior (the entry) and EvaluateOnly ("Evaluate").
    D = struct('PhaseConstraint', 'none', 'PhaseStepBounds', [], 'AmplitudeFit', 'shape', ...
               'SharedVelocity', false, 'Optimizer', 'lsqnonlin', 'PopSize', 5000, ...
               'MaxFunEvals', 2e5, 'Polish', true, 'FocalWeight', 0, 'FitWindowMm', [], ...
               'Lambda', 1e-3, 'LambdaAmp', 1e-2, 'AmpBounds', [0.1 10], 'PlausibleGain', [0.5 2], ...
               'NearFieldPenalty', 1, 'NearFieldMargin', 0.1, 'NearFieldMinMm', 2, ...
               'CWater', 1500, 'Rho', 994, 'Seed', 1);
end

function E = extra_defaults()
% Desired-profile (INTERPOLATE_AXIAL_PROFILES) and k-Wave drive calibration
% (KWAVE_CALIBRATE_DRIVE) options of "Advanced…" and their defaults
    E = struct('Tolerance', 0.01, 'GridDims', '', 'ApplyCorrection', true);
end

function R = ring_defaults()
% FIT_RING_PROFILE and READ_OFFAXIS_SCANS options ("Ring profile options…")
    R = struct('NSub', 4, 'Smooth', 0.3, 'NearFieldSkipMm', 4, 'AmpBounds', [0.1 10], ...
               'Iterations', 6, 'NStarts', 12, 'Seed', 1, 'HeldOut', false, 'OffAxis', '', ...
               'Kind', 'all', 'Weight', 1, 'RMaxMm', 30, 'BinMm', 0.5, 'ZBinMm', 1, 'Centre', 'origin');
end

function v = non_default(v, D)
% Keep only the fields of v that differ from D
    for f = fieldnames(v)'
        if isfield(D, f{1}) && isequal(double_if_numeric(v.(f{1})), double_if_numeric(D.(f{1})))
            v = rmfield(v, f{1});
        end
    end
end

function v = double_if_numeric(v)
% Values loaded from YAML in the form of advanced_defaults (char, double rows)
    if isstring(v); v = char(v); elseif isnumeric(v); v = double(v(:)'); end
end

function items = depth_items(depths)
    items = arrayfun(@(d) sprintf('%.2f', d), depths, 'UniformOutput', false);
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function v = get_field(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function s = num2str_blank(x)
    if isfinite(x); s = num2str(x); else; s = ''; end
end

function v = pair_or_empty(e1, e2)
% [a b] from two text fields; [] when either is blank
    a = str2double(e1.Value); b = str2double(e2.Value);
    if isfinite(a) && isfinite(b); v = sort([a b]); else; v = []; end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end

function p = canonical(p)
    p = char(java.io.File(char(p)).getCanonicalPath());
end

function clear_axes(ax)
% CLA leaves objects with hidden handles (xline, HandleVisibility 'off') and
% keeps limits frozen by an earlier xlim/ylim call
    delete(allchild(ax));
    set(ax, 'XLimMode', 'auto', 'YLimMode', 'auto');
end

function tf = is_absolute(p)
    tf = java.io.File(char(p)).isAbsolute();
end

function s = transducer_serials(equip)
    s = fieldnames(equip.trans);
    s = s(cellfun(@(x) isfield(equip.trans.(x), 'transducer') && ...
                       any(isfield(equip.trans.(x).transducer, {'annular', 'matrix'})), s));
end

function tf = is_matrix(t)
    tf = ~isfield(t.transducer, 'annular') && isfield(t.transducer, 'matrix');
end

function s = default_predict_depths(t)
    if is_matrix(t); s = ''; return; end
    lo = get_field(t, 'min_foc', 0.5 * t.transducer.annular.dist_geom_ep_mm);
    hi = get_field(t, 'max_foc', 1.2 * t.transducer.annular.dist_geom_ep_mm);
    lo = 5 * ceil(lo / 5); hi = 5 * floor(hi / 5);
    if hi <= lo; hi = lo + 5; end
    s = sprintf('%g:5:%g', lo, hi);
end

function [v, ok] = parse_depths(txt)
% 'a:b:c', 'a:c' or a comma/space-separated list
    parts = strsplit(strtrim(char(txt)), ':');
    if numel(parts) == 1
        [v, ok] = parse_list(txt); return;
    end
    x  = str2double(parts);
    ok = any(numel(x) == [2 3]) && all(isfinite(x));
    if ~ok; v = []; return; end
    if numel(x) == 2; v = x(1):x(2); else; v = x(1):x(2):x(3); end
    ok = ~isempty(v) && numel(v) <= 200;
end

function pos = matrix_elem_pos_mm(m)
% Element centres [N x 3] in mm, apex at the origin, beam along +z, from
% CONVERT_TO_ELEMENT_POS so the preview matches the simulation
    [~, p] = evalc('convert_to_element_pos(struct(), struct(''matrix'', m), [0; 0; 0], [0; 0; 0])');
    pos = p' * 1e3;
end

function info = draw_transducer_geometry(axs, axf, t)
% Side (axial) and face views of an annular or matrix transducer entry
    clear_axes(axs); clear_axes(axf); hold(axs, 'on'); hold(axf, 'on');
    tr = t.transducer; info = {};
    blue = [0.0 0.45 0.70];
    if ~is_matrix(t)
        a  = tr.annular; R = a.curv_radius_mm; curved = true;
        id = a.elem_id_mm(:)'; od = a.elem_od_mm(:)'; n = numel(od);
        if numel(id) ~= n || n == 0; error('enter matching inner and outer diameters'); end
        nch = double(get_field(t, 'n_elem', n));
        if nch >= 1 && nch <= n && mod(n, nch) == 0
            grp = kron(1:nch, ones(1, n / nch));
        else
            grp = 1:n; nch = n;
        end
        col = lines(max(nch, 1));
        th  = linspace(0, 2*pi, 180);
        for i = 1:n
            r = linspace(id(i) / 2, od(i) / 2, 40);
            z = R - sqrt(max(R^2 - r.^2, 0));
            plot(axs, z, r, '-', 'Color', col(grp(i), :), 'LineWidth', 5, 'HandleVisibility', 'off');
            plot(axs, z, -r, '-', 'Color', col(grp(i), :), 'LineWidth', 5, 'HandleVisibility', 'off');
            px = od(i)/2 * cos(th); py = od(i)/2 * sin(th);
            if id(i) > 0                                  % annulus: trace the inner edge back
                px = [px, id(i)/2 * cos(fliplr(th))]; py = [py, id(i)/2 * sin(fliplr(th))];
            end
            patch(axf, px, py, col(grp(i), :), 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.3, 'FaceAlpha', 0.85);
        end
        rmax = max(od) / 2; dep = a.dist_geom_ep_mm;
        info{end+1} = sprintf('%d rings in %d channels, aperture Ø %.1f mm (colour = channel).', n, nch, max(od));
        bad = find(od(1:end-1) > id(2:end));
        if ~isempty(bad); info{end+1} = sprintf('Rings overlap after ring(s) %s.', num2str(bad)); end
        face_title = sprintf('Face view: %d rings, %d channels', n, nch);
    else
        m = tr.matrix; R = m.curv_radius_mm; curved = logical(m.is_curved);
        rmax = m.outer_diameter_mm / 2; dep = m.dist_geom_ep_mm;
        pos = matrix_elem_pos_mm(m); N = size(pos, 1);
        r = linspace(-rmax, rmax, 200);
        if curved; zc = R - sqrt(max(R^2 - r.^2, 0)); else; zc = zeros(size(r)); end
        plot(axs, zc, r, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1, 'HandleVisibility', 'off');
        scatter(axs, pos(:, 3), pos(:, 1), 10, blue, 'filled', 'HandleVisibility', 'off');
        % Footprints: rect is width (x) × height (y); disc/bowl use the
        % equal-area diameter, as in CREATE_MATRIX_KARRAY
        w = m.elem_width_mm; h = m.elem_height_mm;
        if strcmp(char(m.elem_shape), 'rect')
            fx = [-1 1 1 -1]' * w/2; fy = [-1 -1 1 1]' * h/2; esize = max(w, h);
        else
            dd = sqrt(4 * w * h / pi); th = linspace(0, 2*pi, 24)';
            fx = dd/2 * cos(th); fy = dd/2 * sin(th); esize = dd;
        end
        if N > 0
            patch(axf, pos(:, 1)' + fx, pos(:, 2)' + fy, blue, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end
        th = linspace(0, 2*pi, 200);
        plot(axf, rmax * cos(th), rmax * sin(th), 'k--');
        pitch = NaN;
        if N > 1
            nn = inf(N, 1);
            for i = 1:N
                dxy = hypot(pos(:, 1) - pos(i, 1), pos(:, 2) - pos(i, 2)); dxy(i) = inf;
                nn(i) = min(dxy);
            end
            pitch = median(nn);
        end
        gs = m.matrix_shape.define_here.grid_shape;
        info{end+1} = sprintf('%d elements (%s grid, %s elements), median pitch %.2f mm.', ...
            N, char(gs.type), char(m.elem_shape), pitch);
        if esize > pitch; info{end+1} = sprintf('Elements (%.2f mm) are larger than the pitch: they overlap.', esize); end
        if strcmp(char(gs.type), 'rect') && gs.rect.sparsity_factor < 1
            info{end+1} = 'sparsity_factor is validated but not applied by PRESTUS; all elements are shown.';
        end
        face_title = sprintf('Face view: %d elements', N);
    end

    % Exit plane, geometric focus, focal range (axial axis = distance from apex)
    z_ep = ternary(curved, R - dep, 0);
    yl = 1.15 * rmax;
    plot(axs, [z_ep z_ep], [-yl yl], 'k:', 'LineWidth', 1, 'DisplayName', 'Exit plane');
    if curved; plot(axs, R, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Geometric focus'); end
    fmin = get_field(t, 'min_foc', []); fmax = get_field(t, 'max_foc', []);
    zmax = max([R * curved, z_ep + 10]);
    if ~isempty(fmin) && ~isempty(fmax)
        plot(axs, z_ep + [fmin fmax], [0 0], '-', 'Color', [0.85 0.37 0.01], 'LineWidth', 4, ...
             'DisplayName', 'Steering range');
        zmax = max(zmax, z_ep + fmax);
    end
    plot(axs, [0 zmax + 5], [0 0], '-', 'Color', [0.85 0.85 0.85], 'HandleVisibility', 'off');
    daspect(axs, [1 1 1]); xlim(axs, [min(-2, z_ep - 5), zmax + 5]); ylim(axs, [-yl yl] * 1.1);
    xlabel(axs, 'Axial distance from apex (mm)'); ylabel(axs, 'Lateral (mm)');
    title(axs, 'Side view'); legend(axs, 'Location', 'southoutside', 'Orientation', 'horizontal');
    grid(axs, 'on'); box(axs, 'on');
    daspect(axf, [1 1 1]); xlim(axf, [-yl yl]); ylim(axf, [-yl yl]);
    xlabel(axf, 'x (mm)'); ylabel(axf, 'y (mm)'); title(axf, face_title); grid(axf, 'on'); box(axf, 'on');
end

function [v, ok] = parse_list(txt)
    tok = strsplit(strtrim(regexprep(char(txt), '[\[\],;\s]+', ' ')), ' ');
    v   = str2double(tok);
    ok  = ~any(isnan(v));
    if isempty(strtrim(char(txt))); v = []; ok = true; end
end

function issues = validate_equipment_entry(t)
% Run CHECK_EQUIPMENT_CONFIG on the entry alone, in a scratch folder
    tmp = tempname; mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, 's'));
    yaml.dumpFile(fullfile(tmp, [t.serial '.yaml']), t);
    [~, rep] = evalc('check_equipment_config(tmp)');
    issues = rep(1).issues;
end

function set_msg(dlg, msg)
% Progress text for CALIBRATION_FIT_DEPTHS (uiprogressdlg is a handle object)
    dlg.Message = msg;
end
