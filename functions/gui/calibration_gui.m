function fig = calibration_gui(options)
% CALIBRATION_GUI  Standalone GUI: fit transducer library entries to hydrophone profiles
%
% Pick a transducer and driving system, load an axial-profile CSV, compare
% the measured profiles with geometric steering and the current library
% entry, fit per-depth refinements (default) or — for debugging — the
% static model (all selected depths jointly), optionally check the per-depth
% refinements in free-water k-Wave (drive corrected to the simulated ISPPA;
% KWAVE_CHECK_REFINEMENT), and save the result to a user library folder
% (transducer.library_path).
%
% Desired depths. "Predict" works at any depths: each gets a profile from
% the measurements (measured, or interpolated between the neighbouring
% measured depths, or — within "extrapolate ≤" — the nearest one shifted;
% INTERPOLATE_AXIAL_PROFILES), and geometric steering, the library entry and
% the unsaved fit are modelled there through the runtime library path
% (EVALUATE_LIBRARY_ENTRY: refinements interpolated between refined depths,
% blended back to the static model outside them). The view shows the
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
% FIT_TRANSDUCER_MODEL and the library helpers, so results match
% REFINE_TRANSDUCER_LIBRARY / BUILD_TRANSDUCER_LIBRARY exactly.
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
%   .SimConfig   - study/simulation YAML for the k-Wave check (grid, code_type,
%                  platform/hpc); default: config_default
%   .ExportFolder - write "Export profiles…" / "Export settings…" CSVs here
%                  under their default names, without a file dialog
%
% Output:
%   fig - uifigure handle
%
% See also: FIT_TRANSDUCER_MODEL, EVALUATE_LIBRARY_ENTRY, INTERPOLATE_AXIAL_PROFILES,
%           LIBRARY_SETTINGS_TABLE, REFINE_TRANSDUCER_LIBRARY, LOAD_TRANSDUCER_FROM_LIBRARY,
%           CHECK_EQUIPMENT_CONFIG, KWAVE_CHECK_REFINEMENT

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
    st.lib = []; st.res = struct(); st.fit_model = []; st.fit_mode = '';
    st.fit_prof = [];       % profiles the current fit was fitted to (.profile_type per depth)
    st.kw = [];   % k-Wave check of the current per-depth fit: .lib, .report, .sims, .res
    st.predicted = false;   % true: st.prof is the desired-depth grid
    st.adv = struct();      % advanced FIT_TRANSDUCER_MODEL options ("Advanced…"); empty = defaults

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

    lab('Fit mode', 13);
    ui.mode = at(uidropdown(L, 'Items', {'Per-depth refinement', 'Static model (joint, debug)'}, ...
        'Tooltip', ['Per-depth: phases, gains and drive at each selected depth (default). ' ...
                    'Static: one depth-independent model over all selected depths, for debugging.'], ...
        'ValueChangedFcn', @(~,~) on_mode()), 13, 2);
    lab('Free parameters', 14);
    fp = at(uigridlayout(L, [2 2], 'Padding', 0, 'RowHeight', {20, 20}, 'RowSpacing', 4), [14 15], 2);
    ui.fp_phase = uicheckbox(fp, 'Text', 'phase', 'Value', true);
    ui.fp_amp   = uicheckbox(fp, 'Text', 'gain', 'Value', true);
    ui.fp_c     = uicheckbox(fp, 'Text', 'c steer', 'Value', false);
    lab('Forward model', 16);
    ui.fwd     = at(uidropdown(L, 'Items', {'rayleigh', 'oneil'}), 16, 2);
    lab('Near-field skip (mm)', 17);
    ui.skip    = at(uieditfield(L, 'numeric', 'Value', 10, 'Limits', [0 50]), 17, 2);
    lab('Loss on', 18);
    ui.domain  = at(uidropdown(L, 'Items', {'intensity', 'amplitude'}), 18, 2);
    lab('Huber δ', 19);
    ui.huber   = at(uieditfield(L, 'numeric', 'Value', 0.8, 'Limits', [1e-3 1]), 19, 2);
    lab('Starts', 20);
    nsb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 96}, 'ColumnSpacing', 4), 20, 2);
    ui.nstarts = uieditfield(nsb, 'numeric', 'Value', 16, 'Limits', [1 200], 'RoundFractionalValues', 'on');
    ui.adv_btn = uibutton(nsb, 'Text', 'Advanced…', 'Tooltip', ['Phase constraints, absolute-intensity fit, ' ...
        'global optimisers, focal weighting (PERFORM_GLOBAL_SEARCH behaviour)'], ...
        'ButtonPushedFcn', @(~,~) on_advanced());
    lab('Reference ISPPA (W/cm²)', 21);
    ui.isppa   = at(uieditfield(L, 'numeric', 'Value', 30, 'Limits', [0.1 1000]), 21, 2);

    lab('k-Wave check', 22);
    kb = at(uigridlayout(L, [1 2], 'Padding', 0, 'ColumnWidth', {'1x', 95}, 'ColumnSpacing', 4), 22, 2);
    ui.kw      = uicheckbox(kb, 'Text', 'free water, after fit', 'Value', false, 'Tooltip', ...
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

    on_mode();
    on_transducer();

    %% Callbacks --------------------------------------------------------------
    function on_mode()
        % c steer changes the steering law of the whole entry: static model only.
        % The k-Wave check applies to per-depth refinements.
        per_depth = ~startsWith(ui.mode.Value, 'Static');
        ui.fp_c.Enable = ternary(per_depth, 'off', 'on');
        ui.kw.Enable   = ternary(per_depth, 'on', 'off');
        ui.kw_grid.Enable = ui.kw.Enable;
    end

    function on_browse_simcfg()
        [f, p] = uigetfile({'*.yaml;*.yml', 'YAML config'}, 'Simulation config for the k-Wave check');
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
        st.meas = []; st.prof = []; st.predicted = false; st.res = struct(); st.fit_model = []; st.fit_prof = []; st.kw = [];
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
        st.fit_model = []; st.fit_prof = []; st.kw = [];
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

    function p = desired_profiles(dep)
        % Profiles at the desired depths from the measurements: measured,
        % interpolated or extrapolated (INTERPOLATE_AXIAL_PROFILES); NaN where
        % none can be derived
        p = interpolate_axial_profiles(st.meas, dep, struct('Align', ui.align.Value, ...
            'MaxExtrapolationMm', ui.extrap.Value, 'NearFieldSkipMm', ui.skip.Value));
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
        zmax = ceil(1.5 * max([dep, tran.transducer.annular.dist_geom_ep_mm]));
        if isempty(st.meas)
            z = (0:0.25:zmax)';
            p = struct('depths_ep_mm', dep, 'z_ep_mm', z, 'I', nan(numel(z), numel(dep)), 'source', '');
            p.profile_type = repmat({'unavailable'}, 1, numel(dep));
        else
            p  = desired_profiles(dep);
            z  = p.z_ep_mm; dz = median(diff(z));
            if z(end) < zmax                              % model beyond the scan; no profile there
                zx = (z(end) + dz:dz:zmax)';
                p.z_ep_mm = [z; zx]; p.I = [p.I; nan(numel(zx), numel(dep))];
            end
        end
        st.prof = p;
        st.predicted = true;
        items = depth_items(dep);
        ui.show.Items = [{'All depths'}, items]; ui.show.Value = 'All depths';
        on_evaluate();
    end

    function on_evaluate()
        % Rebuild all results on view from the current state (profiles,
        % library entry, fit, k-Wave check) and redraw
        if is_matrix(st.equip.trans.(ui.tran.Value)); say('Evaluation supports annular arrays only.'); return; end
        if isempty(st.prof); on_predict(); return; end
        tran = st.equip.trans.(ui.tran.Value);
        base = base_options();
        if st.predicted; base.NearFieldSkipMm = 0; end   % show the full modelled profile
        geo = base; geo.EvaluateOnly = true;
        rt  = base; rt.IsppaWcm2 = ui.isppa.Value;      % runtime path, drive at the reference ISPPA
        res = struct();
        has = any(isfinite(st.prof.I) & st.prof.I > 0, 1);   % depths with a (measured or derived) profile
        pg = st.prof; pg.I(:, ~has) = 1;                 % flat placeholder: the model scale is fitted to it
        res.geometric = pack(fit_transducer_model(tran, pg, geo));
        if ~isempty(st.lib)
            res.library = pack(evaluate_library_entry(tran, st.lib, st.prof, rt));
        end
        if ~isempty(st.fit_model)
            if st.predicted || fit_on_derived()
                % As a study resolves it (on the measured depths too, when the
                % fit was calibrated on desired depths)
                res.fit = pack(evaluate_library_entry(tran, fitted_entry(), st.prof, rt));
            elseif iscell(st.fit_model)
                res.fit = pack_many(st.fit_model);
                if ~isempty(st.kw); res.kwave = st.kw.res; end
            else
                res.fit = pack(st.fit_model);
            end
        end
        if st.predicted
            % The full profile is shown; R² and focus are taken beyond the near
            % field as in the fit (a near-field peak would otherwise count as
            % the focus), R² against the measured or derived profile
            zmax = inf; if ~isempty(st.meas); zmax = max(st.meas.z_ep_mm); end
            for f = fieldnames(res)'
                r = res.(f{1});
                r.Imeas(:, ~has) = NaN;
                r.focus = focus_beyond(r.Imod, r.z, ui.skip.Value);
                r.r2 = r2_beyond(r.Imod, r.Imeas, r.z, ui.skip.Value, zmax);
                res.(f{1}) = r;
            end
        end
        st.res = res;
        has_lib = isfield(st.res, 'library');
        has_fit = isfield(st.res, 'fit');
        if st.predicted
            ferr = @(r) median(abs(r.focus - r.depths));
            txt = sprintf('geometric %.1f mm', ferr(st.res.geometric));
            if has_lib; txt = [txt sprintf(' | library %.1f mm', ferr(st.res.library))]; end
            if has_fit; txt = [txt sprintf(' | fit (unsaved) %.1f mm', ferr(st.res.fit))]; end
            pt = st.prof.profile_type;
            say(sprintf(['Predicted %d depths (profiles: %d measured, %d interpolated, %d extrapolated, ' ...
                '%d none). Median |focus − target|: %s'], numel(st.prof.depths_ep_mm), ...
                nnz(strcmp(pt, 'measured')), nnz(strcmp(pt, 'interpolated')), ...
                nnz(strcmp(pt, 'extrapolated')), nnz(strcmp(pt, 'unavailable')), txt));
            [src, label] = settings_source();
            if ~isempty(src)
                pm = {st.res.(src).settings.phase_model};
                r2 = st.res.(src).r2;
                say(sprintf(['%s at the desired depths: %d refined (within the refined span), %d blended ' ...
                    'back to the static model, %d static/geometric; median R² %.3f on %d profiles ' ...
                    '(phases and drive at %g W/cm² in the table).'], label, ...
                    nnz(strcmp(pm, 'refined')), nnz(strcmp(pm, 'blend')), nnz(~ismember(pm, {'refined', 'blend'})), ...
                    median(r2, 'omitnan'), nnz(isfinite(r2)), ui.isppa.Value));
            end
        else
            lib_txt = '';
            if has_lib; lib_txt = sprintf(' | library %.3f', median(st.res.library.r2, 'omitnan')); end
            say(sprintf('Geometric median R² %.3f%s', median(st.res.geometric.r2, 'omitnan'), lib_txt));
        end
        draw();
    end

    function tf = fit_on_derived()
        tf = ~isempty(st.fit_prof) && any(~strcmp(st.fit_prof.profile_type, 'measured'));
    end

    function on_fit()
        if isempty(st.meas); say('Fitting needs measured profiles: load a CSV first.'); return; end
        tran = st.equip.trans.(ui.tran.Value);
        at_desired = strcmp(ui.calat.Value, 'Desired depths');
        if at_desired
            [dep, ok] = desired_depths();
            if ~ok; say('Desired depths: enter depths such as 20:5:60 or 30, 42.5, 50.'); return; end
            sel  = desired_profiles(dep);
            none = strcmp(sel.profile_type, 'unavailable');
            if any(none)
                say(sprintf(['No profile at %s mm (outside the measured range %.1f–%.1f mm; ' ...
                    'raise "extrapolate ≤" to include them): skipped.'], mat2str(sel.depths_ep_mm(none)), ...
                    min(st.meas.depths_ep_mm), max(st.meas.depths_ep_mm)));
                sel = profile_subset(sel, ~none);
            end
            if isempty(sel.depths_ep_mm); say('None of the desired depths has a profile to calibrate on.'); return; end
        else
            sel = selected_profiles();
            if isempty(sel.depths_ep_mm); say('Select at least one depth.'); return; end
        end
        fo   = base_options();
        if startsWith(ui.mode.Value, 'Static') && ~strcmp(get_field(fo, 'PhaseConstraint', 'none'), 'none')
            say('Phase constraints act on one depth at a time: use per-depth refinement (Advanced…).'); return;
        end
        % No z shift: it is not applied at run time, so phases fitted with it
        % would match the measurement only after a shift the entry never gets
        free = {'phase', 'amp', 'c_steer'};
        free = free([ui.fp_phase.Value, ui.fp_amp.Value, ui.fp_c.Value]);
        dlg  = uiprogressdlg(fig, 'Title', 'Fitting', 'Indeterminate', 'on');
        cleanup = onCleanup(@() close(dlg));
        t0 = tic;
        st.kw = [];
        if startsWith(ui.mode.Value, 'Static')
            fo.FreeParams = free;
            st.fit_model = fit_transducer_model(tran, sel, fo);
            st.fit_mode  = 'static';
            r2 = st.fit_model.metrics.r2;
        else
            % z shift stays at 0 (also for the prior). The prior only pins the
            % unobservable common phase (default Lambda): a stronger pull stops
            % the fit from reaching split-focus / near-field-peak solutions.
            fo.FreeParams = intersect(free, {'phase', 'amp'}, 'stable');
            if ~isempty(st.lib); fo.Prior = library_static_model(st.lib); end
            fo.Prior.z_shift_mm = 0;
            models = cell(1, numel(sel.depths_ep_mm));
            for k = 1:numel(sel.depths_ep_mm)
                dlg.Message = sprintf('Depth %.1f mm, %s profile (%d/%d)', sel.depths_ep_mm(k), ...
                    sel.profile_type{k}, k, numel(sel.depths_ep_mm));
                pk = sel; pk.depths_ep_mm = sel.depths_ep_mm(k); pk.I = sel.I(:, k);
                models{k} = fit_transducer_model(tran, pk, fo);
                models{k}.profile_type = sel.profile_type{k};   % stored as refinement.profile_type
            end
            st.fit_model = models;
            st.fit_mode  = 'per_depth';
            r2 = cellfun(@(m) m.metrics.r2, models);
        end
        st.fit_prof = sel;
        pt = sel.profile_type;
        say(sprintf('Fit (%s, %s) on %d profiles (%d measured, %d interpolated, %d extrapolated) in %.1f s: median R² %.3f, min %.3f', ...
            st.fit_mode, strjoin(fo.FreeParams, '+'), numel(pt), nnz(strcmp(pt, 'measured')), ...
            nnz(strcmp(pt, 'interpolated')), nnz(strcmp(pt, 'extrapolated')), toc(t0), median(r2), min(r2)));
        if strcmp(st.fit_mode, 'per_depth') && ui.kw.Value
            run_kwave_check(sel, dlg);
        end
        if at_desired || st.predicted
            on_predict();    % show the result on the desired depths, as a study resolves it
        else
            show_measured(); on_evaluate();
        end
    end

    function run_kwave_check(sel, dlg)
        % Free-water k-Wave of the refined depths through the runtime library path
        dlg.Indeterminate = 'off'; dlg.Title = 'k-Wave check';
        ko = struct('EquipParam', st.equip, 'Axisymmetric', strcmp(ui.kw_grid.Value, 'axisymmetric'), ...
            'NearFieldSkipMm', ui.skip.Value, 'WorkDir', fullfile(tempdir, 'prestus_kwave_check', combo_name()), ...
            'ProgressFcn', @(k, n, d) set(dlg, 'Value', (k - 1) / n, ...
                'Message', sprintf('Free-water simulation at %.1f mm (%d/%d)', d, k, n)));
        if ~isempty(strtrim(ui.simcfg.Value)); ko.SimConfig = strtrim(ui.simcfg.Value); end
        t1 = tic;
        try
            [lib, rep, sims] = kwave_check_refinement(entry_from_fit(), sel, ko);
        catch ME
            say(['k-Wave check failed: ' ME.message]); return;
        end
        st.kw = struct('lib', lib, 'report', rep, 'sims', sims, ...
            'res', pack_kwave(sims, rep, st.fit_model{1}.z_ep_mm, ui.skip.Value));
        say(sprintf(['k-Wave check (%d depths, %.0f s): ISPPA %s W/cm² at the analytical drive ' ...
            '(target %g) → drive ×%s; median R² %.3f'], height(rep), toc(t1), ...
            strjoin(compose('%.1f', rep.isppa_kwave_wcm2'), ' '), ui.isppa.Value, ...
            strjoin(compose('%.2f', rep.amp_scale'), ' '), median(rep.r2_kwave)));
    end

    function on_save()
        if isempty(st.fit_model); say('Nothing to save — run a fit first.'); return; end
        out = ui.out.Value;
        if isempty(out); on_browse_out(); out = ui.out.Value; end
        if isempty(out); return; end
        if strcmp(canonical(out), canonical(st.lib_path))
            sel = uiconfirm(fig, ['This overwrites the SHIPPED library entry. ' ...
                'Use a user folder and transducer.library_path unless you maintain the library.'], ...
                'Shipped library', 'Options', {'Overwrite shipped entry', 'Cancel'}, 'DefaultOption', 2);
            if ~strcmp(sel, 'Overwrite shipped entry'); return; end
        end
        lib = fitted_entry();
        lib.meta.calibration_type = 'user_empirical';
        lib.meta.refinement_mode  = st.fit_mode;
        lib.meta.last_updated     = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
        try; lib.meta.prestus_version = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end
        if ~isfolder(out); mkdir(out); end
        if fit_on_derived()
            % The profiles the entry was fitted to, next to it (as REFINE_TRANSDUCER_LIBRARY)
            pcsv = fullfile(out, [combo_name() '_calibration_profiles.csv']);
            write_axial_profiles(st.fit_prof, pcsv, struct('combo', combo_name()));
            if strcmp(st.fit_mode, 'per_depth'); lib.refinement.calibration_profiles_csv = pcsv; end
            say(sprintf('Calibration profiles (measured + derived): %s', pcsv));
        end
        f = fullfile(out, [combo_name() '.yaml']);
        yaml.dumpFile(f, lib);
        say(sprintf('Saved %s — set transducer.library_path: %s', f, out));
    end

    function on_export_profiles()
        % The profiles of the last fit, else the desired-depth profiles on view
        if ~isempty(st.fit_prof)
            p = st.fit_prof; what = 'calibration';
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
        if ~isempty(st.fit_model)
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
        % drives if checked, else the analytical fit
        if ~isempty(st.kw); lib = st.kw.lib; else; lib = entry_from_fit(); end
    end

    function lib = entry_from_fit()
        % Library entry from the current entry and fit (not yet written)
        lib = st.lib;
        if isempty(lib); lib = struct(); end
        lib.meta.tran_serial = ui.tran.Value;
        if ~strcmp(ui.ds.Value, '(none)'); lib.meta.driving_system = ui.ds.Value; end
        if strcmp(st.fit_mode, 'static')
            lib = library_entry_from_model(lib, st.fit_model, ui.isppa.Value);
        else
            lib = library_refinement_from_fits(lib, st.fit_model, ui.isppa.Value, [], ui.csv.Value);
        end
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
                    fm = focus_beyond(r.Imeas(:, jt), r.z, ui.skip.Value);
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
        m = st.fit_model;
        if isempty(m); ui.params.Data = {}; return; end
        if iscell(m)
            ph = cell2mat(cellfun(@(x) x.elem_phase_correction_deg(:), m, 'UniformOutput', false));
            ui.params.ColumnName = [{'Element'}, arrayfun(@(x) sprintf('%.1f mm', x.fit_depths_ep_mm), ...
                [m{:}], 'UniformOutput', false)];
            ui.params.Data = num2cell([(1:size(ph, 1))', round(ph, 1)]);
        else
            ui.params.ColumnName = {'Element', 'Phase offset (°)', 'Gain'};
            ui.params.Data = num2cell([(1:numel(m.elem_amp_weight))', ...
                round(m.elem_phase_correction_deg(:), 1), round(m.elem_amp_weight(:), 3)]);
            say(sprintf('c_steer %.0f m/s | cost %.3g (%d starts)', ...
                m.steering_sound_speed, m.cost, numel(m.start_costs)));
        end
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
        % Advanced FIT_TRANSDUCER_MODEL options. Only values that differ from the
        % fitter defaults are stored (st.adv), so "Reset" restores the default fit.
        D = struct('PhaseConstraint', 'none', 'PhaseStepBounds', [], 'AmplitudeFit', 'shape', ...
                   'SharedVelocity', false, 'Optimizer', 'lsqnonlin', 'PopSize', 5000, ...
                   'MaxFunEvals', 2e5, 'Polish', true, 'FocalWeight', 0, 'FitWindowMm', [], ...
                   'Lambda', 1e-3, 'LambdaAmp', 1e-2, 'ZShiftBounds', [-3 3]);
        cur = D; for f = fieldnames(st.adv)'; cur.(f{1}) = st.adv.(f{1}); end
        d = uifigure('Name', 'Advanced fit options', 'Position', [200 200 430 470], 'WindowStyle', 'modal');
        G = uigridlayout(d, [15 2], 'RowHeight', repmat({24}, 1, 15), 'ColumnWidth', {190, '1x'});
        al = @(txt, r) at(uilabel(G, 'Text', txt), r, 1);
        num = @(v, r) at(uieditfield(G, 'numeric', 'Value', v), r, 2);
        at(uilabel(G, 'Text', 'Phase model', 'FontWeight', 'bold'), 1, [1 2]);
        al('Phase constraint', 2);
        w.pc = at(uidropdown(G, 'Items', {'none', 'linear', 'monotonic'}, 'Value', cur.PhaseConstraint), 2, 2);
        sb = cur.PhaseStepBounds; if isempty(sb); sb = [NaN NaN]; end
        al('Step bounds (rad; blank = default)', 3);
        g2 = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 3, 2);
        w.sb1 = uieditfield(g2, 'text', 'Value', num2str_blank(sb(1)));
        w.sb2 = uieditfield(g2, 'text', 'Value', num2str_blank(sb(2)));
        al('Phase prior λ', 4);                 w.lam  = num(cur.Lambda, 4);
        al('Gain prior λ', 5);                  w.lama = num(cur.LambdaAmp, 5);
        at(uilabel(G, 'Text', 'Amplitude and loss', 'FontWeight', 'bold'), 6, [1 2]);
        al('Amplitude fit', 7);
        w.af = at(uidropdown(G, 'Items', {'shape', 'absolute'}, 'Value', cur.AmplitudeFit), 7, 2);
        w.sv = at(uicheckbox(G, 'Text', 'one drive for all depths (absolute)', 'Value', cur.SharedVelocity), 8, [1 2]);
        al('Focal weight (0 = uniform, ≥1)', 9); w.fw = num(cur.FocalWeight, 9);
        fw = cur.FitWindowMm; if isempty(fw); fw = [NaN NaN]; end
        al('Fit window (mm EP; blank = all)', 10);
        g3 = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 10, 2);
        w.fw1 = uieditfield(g3, 'text', 'Value', num2str_blank(fw(1)));
        w.fw2 = uieditfield(g3, 'text', 'Value', num2str_blank(fw(2)));
        al('z shift bounds (mm)', 11);
        g4 = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 11, 2);
        w.zb1 = uieditfield(g4, 'numeric', 'Value', cur.ZShiftBounds(1));
        w.zb2 = uieditfield(g4, 'numeric', 'Value', cur.ZShiftBounds(2));
        at(uilabel(G, 'Text', 'Optimiser', 'FontWeight', 'bold'), 12, [1 2]);
        g5 = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4, 'ColumnWidth', {'1x', 80}), 13, [1 2]);
        w.opt = uidropdown(g5, 'Items', {'lsqnonlin', 'fexminimize', 'globalsearch'}, 'Value', cur.Optimizer);
        w.pol = uicheckbox(g5, 'Text', 'polish', 'Value', cur.Polish);
        al('FEX population / max evals', 14);
        g6 = at(uigridlayout(G, [1 2], 'Padding', 0, 'ColumnSpacing', 4), 14, 2);
        w.pop = uieditfield(g6, 'numeric', 'Value', cur.PopSize, 'Limits', [10 1e5], 'RoundFractionalValues', 'on');
        w.mfe = uieditfield(g6, 'numeric', 'Value', cur.MaxFunEvals, 'Limits', [100 1e8], 'RoundFractionalValues', 'on');
        bb = at(uigridlayout(G, [1 3], 'Padding', 0), 15, [1 2]);
        uibutton(bb, 'Text', 'Reset', 'ButtonPushedFcn', @(~,~) apply(struct()));
        uibutton(bb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) close(d));
        uibutton(bb, 'Text', 'OK', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) apply(read()));

        function v = read()
            v = struct('PhaseConstraint', w.pc.Value, 'PhaseStepBounds', pair_or_empty(w.sb1, w.sb2), ...
                'AmplitudeFit', w.af.Value, 'SharedVelocity', w.sv.Value, 'Optimizer', w.opt.Value, ...
                'PopSize', w.pop.Value, 'MaxFunEvals', w.mfe.Value, 'Polish', w.pol.Value, ...
                'FocalWeight', w.fw.Value, 'FitWindowMm', pair_or_empty(w.fw1, w.fw2), ...
                'Lambda', w.lam.Value, 'LambdaAmp', w.lama.Value, 'ZShiftBounds', [w.zb1.Value w.zb2.Value]);
            for f = fieldnames(v)'   % keep only non-defaults
                if isequal(v.(f{1}), D.(f{1})); v = rmfield(v, f{1}); end
            end
        end
        function apply(v)
            st.adv = v;
            n = numel(fieldnames(v));
            ui.adv_btn.Text = ternary(n > 0, sprintf('Advanced (%d)', n), 'Advanced…');
            if n > 0
                say(sprintf('Advanced fit options: %s', strjoin(fieldnames(v), ', ')));
            else
                say('Advanced fit options reset to defaults.');
            end
            close(d);
        end
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
function r = pack(m)
    r = struct('depths', m.fit_depths_ep_mm(:)', 'z', m.z_ep_mm, 'Imeas', m.I_meas, ...
               'Imod', m.I_model, 'r2', m.metrics.r2(:)');
    r.focus = arrayfun(@(j) compute_flhm_center(m.I_model(:, j), m.z_ep_mm), 1:size(m.I_model, 2));
    r.settings = [];   % runtime settings per depth (EVALUATE_LIBRARY_ENTRY results only)
    if isfield(m, 'settings'); r.settings = m.settings; end
end

function r = pack_kwave(sims, rep, z, skip_mm)
% k-Wave on-axis profiles on the fit's z grid, peak-normalised beyond the near field
    r = struct('depths', [sims.depth_ep_mm], 'z', z, 'Imeas', [], 'Imod', zeros(numel(z), numel(sims)), ...
               'r2', rep.r2_kwave(:)');
    for k = 1:numel(sims)
        I = interp1(sims(k).z_ep_mm, sims(k).I_wcm2, z, 'linear', 0);
        r.Imod(:, k) = I / max(max(I(z >= skip_mm)), eps);
    end
    r.focus = arrayfun(@(j) compute_flhm_center(r.Imod(:, j), z), 1:numel(sims));
    r.settings = [];
end

function r = pack_many(models)
    r = pack(models{1});
    for k = 2:numel(models)
        q = pack(models{k});
        r.depths = [r.depths, q.depths]; r.Imeas = [r.Imeas, q.Imeas];
        r.Imod = [r.Imod, q.Imod]; r.r2 = [r.r2, q.r2]; r.focus = [r.focus, q.focus];
    end
end

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

function f = focus_beyond(I, z, skip_mm)
% FLHM centre of each column of I, ignoring z < skip_mm (near field)
    k = z >= skip_mm;
    f = arrayfun(@(j) compute_flhm_center(I(k, j), z(k)), 1:size(I, 2));
end

function r2 = r2_beyond(Imod, Imeas, z, skip_mm, zmax)
% Profile-shape R² (peak-normalised intensity) beyond the near field and
% within the measured scan, as FIT_TRANSDUCER_MODEL; NaN without a profile
    k  = z >= skip_mm & z <= zmax;
    r2 = nan(1, size(Imod, 2));
    for j = 1:size(Imod, 2)
        yt = Imeas(k, j); ym = Imod(k, j);
        ok = isfinite(yt);
        if ~any(ok) || max(yt(ok)) <= 0; continue; end
        yt = yt(ok) / max(yt(ok)); ym = ym(ok) / max(max(ym(ok)), eps);
        r2(j) = 1 - sum((yt - ym).^2) / max(sum((yt - mean(yt)).^2), eps);
    end
end

function p = profile_subset(p, keep)
    p.depths_ep_mm  = p.depths_ep_mm(keep);
    p.I             = p.I(:, keep);
    p.profile_type  = p.profile_type(keep);
    p.source_depths = p.source_depths(keep, :);
    p.weights       = p.weights(keep, :);
    p.shift_mm      = p.shift_mm(keep, :);
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
