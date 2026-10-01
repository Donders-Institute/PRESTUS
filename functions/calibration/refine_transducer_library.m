function [out_file, fits, report] = refine_transducer_library(combo_name, profile_csv, out_library_path, options)
% REFINE_TRANSDUCER_LIBRARY  Optional refinement of a library entry to your own hydrophone profiles
%
% Fits the transducer model (FIT_TRANSDUCER_MODEL) to measured axial profiles
% and writes the result to a library YAML in out_library_path — a user
% folder, so the shipped library in config/calibration/ is not modified.
% Use the result by setting transducer.library_path in the study config.
%
% Two modes:
%   'per_depth' — (default) fit phase offsets and gains at each given depth,
%                 regularised towards the entry's static model (geometric
%                 steering for a geometric entry), and store them with the
%                 drive that gives Intensity at that depth under
%                 'refinement'. load_transducer_from_library uses them at and
%                 between those depths and fades back to the static model
%                 outside (refinement.blend_mm).
%   'static'    — (debugging) refit the depth-independent model jointly over
%                 the given depths: per-element phase offsets, gains and
%                 (optionally) the steering sound speed. Replaces global_model.
%
% Optionally (KWave = true) each refined depth is then checked in a
% free-water k-Wave simulation of the runtime library path, and its drive
% corrected to the simulated ISPPA (KWAVE_CHECK_REFINEMENT).
%
% Use as:
%   out_file = refine_transducer_library(combo_name, profile_csv, out_library_path)
%   [out_file, fits, report] = refine_transducer_library(..., options)
%
% Typical use — measure axial profiles at the depths you will stimulate, then:
%   refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', ...
%       '/project/my_study/transducer_library', struct('Plot', true));
% and set transducer.library_path to that folder in the study config.
% See examples/demo_transducer_refinement.m.
%
% Input:
%   combo_name       - '{tran_serial}_{ds_serial}', e.g. 'CTX_250_001_SC_105_010'
%   profile_csv      - axial profile CSV (doc_transducer_library.md, Appendix)
%   out_library_path - folder to write {combo_name}.yaml into (created if needed).
%                      An existing entry there is updated; otherwise the shipped
%                      entry (config/calibration/) is copied first.
%   options          - struct, all optional:
%     .Mode           - 'per_depth' (default) | 'static'
%     .Depths         - depths [mm from EP] to calibrate (default: all in the CSV).
%                       A depth within 0.5 mm of a measured one uses that
%                       measurement; any other is interpolated from the
%                       neighbouring measurements (INTERPOLATE_AXIAL_PROFILES);
%                       depths that cannot be derived are skipped with a warning
%     .InterpolateOptions - passed to INTERPOLATE_AXIAL_PROFILES (Align,
%                       MaxExtrapolationMm; default no extrapolation)
%     .ExportProfiles - true | false | CSV path: write the profiles fitted to
%                       (measured and derived, with a provenance sidecar);
%                       default true when any depth is derived, to
%                       {combo_name}_calibration_profiles.csv next to the YAML
%     .ExportSettings - true | CSV path: write the per-depth settings the entry
%                       resolves to (LIBRARY_SETTINGS_TABLE) to
%                       {combo_name}_settings.csv (default false)
%     .SettingsDepths - depths for ExportSettings (default: the calibrated depths)
%     .Intensity      - reference ISPPA for the amplitude [W/cm²] (default 30)
%     .FitOptions     - passed to FIT_TRANSDUCER_MODEL (FreeParams, ForwardModel,
%                       NearFieldSkipMm, HuberDelta, Lambda, NStarts, ...).
%                       per_depth defaults: FreeParams {'phase','amp'}, z shift
%                       fixed at 0, Domain 'intensity', HuberDelta 0.8, Prior =
%                       the entry's static model. Adding 'z_shift' fits an axial
%                       offset as a diagnostic (refinement.fit_z_shift_mm); it is
%                       not applied at run time, so the default keeps it at 0.
%     .BlendMm        - refinement.blend_mm for per_depth (default 5)
%     .KWave          - true: k-Wave free-water check of each refined depth
%                       (per_depth only; default false)
%     .KWaveOptions   - options for KWAVE_CHECK_REFINEMENT (SimConfig, WorkDir, ...)
%     .Plot           - true: save {combo_name}_refinement.png next to the YAML
%                       (measured vs library vs refined, per depth)
%     .EquipParam     - struct from load_equipment_config (default: loaded)
%     .ShippedLibrary - source of the base entry (default config/calibration/)
%
% Output:
%   out_file - path of the written YAML
%   fits     - FIT_TRANSDUCER_MODEL result(s): one model ('static') or one per depth
%   report   - table per calibrated depth: profile type, R², |focal error| and
%              near-field excess of the entry before (library) and after refinement
%              (on derived profiles these compare with the derived profile)
%
% See also: FIT_TRANSDUCER_MODEL, LOAD_TRANSDUCER_FROM_LIBRARY, KWAVE_CHECK_REFINEMENT,
%           LIBRARY_ENTRY_FROM_MODEL, READ_AXIAL_PROFILES, CALIBRATION_GUI

    if nargin < 4; options = struct(); end
    mode      = opt(options, 'Mode', 'per_depth');
    intensity = opt(options, 'Intensity', 30);
    fit_opts  = opt(options, 'FitOptions', struct());
    blend_mm  = opt(options, 'BlendMm', 5);
    do_plot   = opt(options, 'Plot', false);
    equip     = opt(options, 'EquipParam', []);
    shipped   = opt(options, 'ShippedLibrary', fullfile(get_prestus_path(), 'config', 'calibration'));
    assert(ismember(mode, {'per_depth', 'static'}), ...
        'refine_transducer_library: Mode must be ''per_depth'' or ''static''.');
    if isempty(equip); equip = load_equipment_config(); end

    %% Base entry ----------------------------------------------------------
    if ~isfolder(out_library_path); mkdir(out_library_path); end
    out_file  = fullfile(out_library_path, [combo_name '.yaml']);
    base_file = fullfile(shipped, [combo_name '.yaml']);
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    elseif isfile(base_file)
        lib = yaml.loadFile(base_file, 'ConvertToArray', true);
    else
        error(['refine_transducer_library: no entry for ''%s'' in %s or %s.\n' ...
            'Create one with build_transducer_library first.'], combo_name, out_library_path, shipped);
    end
    tran = equip.trans.(char(lib.meta.tran_serial));

    %% Profiles ------------------------------------------------------------
    % A requested depth within 0.5 mm of a measured one uses that measurement;
    % any other is interpolated from the neighbouring measurements (or, within
    % MaxExtrapolationMm of the measured range, extrapolated)
    meas = read_axial_profiles(profile_csv);
    use  = reshape(opt(options, 'Depths', meas.depths_ep_mm), 1, []);
    near = arrayfun(@(u) find(abs(meas.depths_ep_mm - u) < 0.5, 1), use, 'UniformOutput', false);
    is_m = ~cellfun(@isempty, near);
    want = unique([meas.depths_ep_mm([near{is_m}]), use(~is_m)]);
    io   = opt(options, 'InterpolateOptions', struct());
    if ~isfield(io, 'NearFieldSkipMm'); io.NearFieldSkipMm = opt(fit_opts, 'NearFieldSkipMm', 10); end
    prof = interpolate_axial_profiles(meas, want, io);
    bad  = strcmp(prof.profile_type, 'unavailable');
    if any(bad)
        warning('refine_transducer_library:outsideProfiles', ...
            ['Depth(s) %s mm lie outside the measured range [%.1f, %.1f] mm (MaxExtrapolationMm = %g) ' ...
             'and are skipped.'], mat2str(prof.depths_ep_mm(bad)), min(meas.depths_ep_mm), ...
            max(meas.depths_ep_mm), opt(io, 'MaxExtrapolationMm', 0));
        prof = profile_subset(prof, ~bad);
    end
    if isempty(prof.depths_ep_mm)
        error('refine_transducer_library: none of the requested depths can be taken from %s (available: %s).', ...
            profile_csv, mat2str(meas.depths_ep_mm));
    end
    derived = ~strcmp(prof.profile_type, 'measured');
    if any(derived)
        fprintf('refine_transducer_library: %d depth(s) calibrated on derived profiles: %s\n', nnz(derived), ...
            strjoin(compose('%.2f mm (%s)', prof.depths_ep_mm(derived)', string(prof.profile_type(derived))'), ', '));
    end
    if isfield(tran, 'min_foc') && isfield(tran, 'max_foc') && ...
            any(prof.depths_ep_mm < tran.min_foc | prof.depths_ep_mm > tran.max_foc)
        warning('refine_transducer_library:outOfRange', ...
            'Some profile depths lie outside [min_foc, max_foc] = [%.1f, %.1f] mm.', tran.min_foc, tran.max_foc);
    end

    % Library entry as shipped/found, scored on these profiles as a study
    % resolves it (static model + any existing refinements; for the report)
    prior    = library_static_model(lib);
    eval_opt = struct();
    if isfield(fit_opts, 'ForwardModel'); eval_opt.ForwardModel = fit_opts.ForwardModel; end
    before   = evaluate_library_entry(tran, lib, prof, eval_opt);
    if ~isfield(fit_opts, 'FreeParams') && strcmp(mode, 'static')
        fit_opts.FreeParams = {'phase', 'amp'};   % the shipped library model
    end

    %% Fit -----------------------------------------------------------------
    switch mode
        case 'static'
            fits = fit_transducer_model(tran, prof, fit_opts);
            lib  = library_entry_from_model(lib, fits, intensity);

        case 'per_depth'
            % Per-depth defaults as in CALIBRATION_GUI: phase and gain,
            % intensity-domain loss, z shift fixed at 0 — the loader does not
            % apply a fitted shift, and phases fitted with one match the
            % hydrophone only after that shift
            prior.z_shift_mm = 0;
            po = struct('FreeParams', {{'phase', 'amp'}}, 'Domain', 'intensity', ...
                        'HuberDelta', 0.8, 'Prior', prior);
            fn = fieldnames(fit_opts);
            for k = 1:numel(fn); po.(fn{k}) = fit_opts.(fn{k}); end
            fits = cell(1, numel(prof.depths_ep_mm));
            for k = 1:numel(prof.depths_ep_mm)
                pk = prof; pk.depths_ep_mm = prof.depths_ep_mm(k); pk.I = prof.I(:, k);
                fits{k} = fit_transducer_model(tran, pk, po);
                fits{k}.profile_type = prof.profile_type{k};   % stored as refinement.profile_type
            end
            lib = library_refinement_from_fits(lib, fits, intensity, blend_mm, profile_csv);
            if opt(options, 'KWave', false)
                ko = opt(options, 'KWaveOptions', struct());
                ko.EquipParam = equip;
                [lib, kw] = kwave_check_refinement(lib, prof, ko);
                disp(kw);
            end
    end

    %% Report -----------------------------------------------------------------
    if strcmp(mode, 'static')
        after = fits.metrics; I_after = fits.I_model; z_after = fits.z_ep_mm;
    else
        after = vertcat_metrics(fits);
        I_after = cell2mat(cellfun(@(f) f.I_model, fits, 'UniformOutput', false)); z_after = fits{1}.z_ep_mm;
    end
    report = table(prof.depths_ep_mm(:), string(prof.profile_type(:)), before.metrics.r2, after.r2, ...
        abs(before.metrics.focal_err_mm), abs(after.focal_err_mm), after.nearfield_excess, ...
        'VariableNames', {'depth_ep_mm', 'profile_type', 'r2_library', 'r2_refined', 'abs_focal_err_library_mm', ...
        'abs_focal_err_refined_mm', 'nearfield_excess_refined'});
    disp(report);

    %% Exports ---------------------------------------------------------------
    % The profiles the entry was fitted to (always when any is derived), so
    % an interpolated calibration can be traced and reproduced
    prof_csv = opt(options, 'ExportProfiles', any(derived));
    if islogical(prof_csv) || isnumeric(prof_csv)
        if prof_csv; prof_csv = fullfile(out_library_path, [combo_name '_calibration_profiles.csv']); else; prof_csv = ''; end
    end
    if ~isempty(prof_csv)
        write_axial_profiles(prof, prof_csv, struct('combo', combo_name));
        fprintf('Calibration profiles: %s\n', prof_csv);
        if strcmp(mode, 'per_depth'); lib.refinement.calibration_profiles_csv = prof_csv; end
    end

    %% Provenance and write ------------------------------------------------
    lib.meta.calibration_type = 'user_empirical';
    lib.meta.last_updated     = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    try; [~, lib.meta.prestus_version] = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end
    lib.meta.refinement_mode  = mode;
    lib.meta.refinement_report = struct('depths_ep_mm', report.depth_ep_mm', ...
        'r2_library', report.r2_library', 'r2_refined', report.r2_refined');
    yaml.dumpFile(out_file, lib);
    fprintf('refine_transducer_library (%s): %s\n', mode, out_file);

    % Per-depth settings the written entry resolves to (optional)
    set_csv = opt(options, 'ExportSettings', false);
    if islogical(set_csv) || isnumeric(set_csv)
        if set_csv; set_csv = fullfile(out_library_path, [combo_name '_settings.csv']); else; set_csv = ''; end
    end
    if ~isempty(set_csv)
        sd = opt(options, 'SettingsDepths', prof.depths_ep_mm);
        so = struct('IsppaWcm2', intensity, 'CsvPath', set_csv, ...
                    'NearFieldSkipMm', opt(fit_opts, 'NearFieldSkipMm', 10));
        if isfield(fit_opts, 'ForwardModel'); so.ForwardModel = fit_opts.ForwardModel; end
        library_settings_table(tran, lib, sd, so);
        fprintf('Settings per depth: %s\n', set_csv);
    end

    if do_plot
        n  = numel(prof.depths_ep_mm);
        fh = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 min(1800, 330 * n) 320]);
        tl = tiledlayout(fh, 1, n, 'TileSpacing', 'compact');
        for k = 1:n
            ax = nexttile(tl); hold(ax, 'on');
            plot(ax, before.z_ep_mm, before.I_meas(:, k), 'k-', 'LineWidth', 1.6, 'DisplayName', 'Hydrophone');
            plot(ax, before.z_ep_mm, before.I_model(:, k), '-', 'Color', [0 0.45 0.70], ...
                'DisplayName', sprintf('Library (R² %.2f)', report.r2_library(k)));
            plot(ax, z_after, I_after(:, k), '-', 'Color', [0.85 0.37 0.01], 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Refined (R² %.2f)', report.r2_refined(k)));
            xline(ax, prof.depths_ep_mm(k), ':', 'HandleVisibility', 'off');
            ttl = sprintf('%.1f mm', prof.depths_ep_mm(k));
            if ~strcmp(prof.profile_type{k}, 'measured'); ttl = sprintf('%s (%s profile)', ttl, prof.profile_type{k}); end
            title(ax, ttl, 'FontWeight', 'normal');
            xlabel(ax, 'mm from exit plane'); grid(ax, 'on'); box(ax, 'on');
            if k == 1; ylabel(ax, 'Normalised intensity'); legend(ax, 'Location', 'northeast', 'FontSize', 7); end
        end
        title(tl, sprintf('%s — %s refinement', strrep(combo_name, '_', '\_'), mode));
        png = fullfile(out_library_path, [combo_name '_refinement.png']);
        exportgraphics(fh, png, 'Resolution', 150); close(fh);
        fprintf('Report figure: %s\n', png);
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

function M = vertcat_metrics(fits)
    M = fits{1}.metrics;
    for k = 2:numel(fits); M = [M; fits{k}.metrics]; end %#ok<AGROW>
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
