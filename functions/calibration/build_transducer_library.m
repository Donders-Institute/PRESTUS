%% build_transducer_library  Build or refresh the full transducer calibration library
%
% For each equipment combo in config/equipment/ this script either:
%   - Builds an empirical entry when axial characterisation data
%     (char_data_path) is available: one refinement per measured depth, or
%   - Computes a pure geometry-only entry from ring dimensions alone.
%
% Both paths produce a global_model section in config/calibration/{combo}.yaml
% that load_transducer_from_library can interpolate at any focal depth.
%
% QUICK START
%   1. Put your axial profile CSVs in local/profiles/ (default CHAR_DATA_ROOT;
%      layout in doc_calibration.md, Appendix), or set CHAR_DATA_ROOT to
%      another folder ('' if char_data_path entries are absolute).
%   2. Optionally restrict processing with COMBOS_TO_RUN.
%   3. Run the script. Results are written to config/calibration/ inside the repo.
%
% ADDING A NEW TRANSDUCER
%   Add a YAML to config/equipment/, set char_data_path if you have hydrophone
%   data, then re-run this script. The new combo is processed automatically.
%
% SINGLE COMBOS / NEW PROFILES
%   REFINE_TRANSDUCER_LIBRARY (or CALIBRATION_GUI) refits one combo's per-depth
%   refinements, optionally with the free-water k-Wave check
%   (KWAVE_CALIBRATE_DRIVE).
%
% See also: LOAD_TRANSDUCER_FROM_LIBRARY, REFINE_TRANSDUCER_LIBRARY,
%           FIT_TRANSDUCER_MODEL, VALIDATE_TRANSDUCER_LIBRARY, CALIBRATION_GUI

%% ── Paths ────────────────────────────────────────────────────────────────────
repo_root    = fullfile(fileparts(mfilename('fullpath')), '..', '..');
equip_path   = fullfile(repo_root, 'config', 'equipment');
if ~exist('LIBRARY_PATH', 'var'); LIBRARY_PATH = fullfile(repo_root, 'config', 'calibration'); end
library_path = LIBRARY_PATH;   % output folder (default: the shipped library)

% Only in a fresh session: once PRESTUS is on the path, the path is left as set.
if isempty(which('load_parameters'))
    addpath(fullfile(repo_root, 'functions', 'helper'));
    safe_addpath(fullfile(repo_root, 'functions'));
    safe_addpath(fullfile(repo_root, 'external'));
end

%% ── CONFIGURATION ────────────────────────────────────────────────────────────
% Variables are only set here when not already defined in the caller's workspace.
% Site-specific wrapper scripts (e.g. build_transducer_library_donders.m) can
% pre-set any of these before calling run() on this file to override defaults.

% Which combos to process.
%   'all'  — every combo found in config/equipment/
%   cell   — subset, e.g. {'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'}
if ~exist('COMBOS_TO_RUN',      'var'); COMBOS_TO_RUN      = 'all'; end

% Root folder for axial-profile CSVs referenced by char_data_path in the
% equipment YAML (default: local/profiles/ in the repo, git-ignored). Set to
% '' if those paths are already absolute.
if ~exist('CHAR_DATA_ROOT',     'var'); CHAR_DATA_ROOT     = fullfile(repo_root, 'local', 'profiles'); end

% Reference ISPPA [W/cm²] the stored drive is scaled to (empirical combos;
% only the first value is used).
if ~exist('DESIRED_INTENSITIES','var'); DESIRED_INTENSITIES = [30]; end

% Library model: one refinement per measured depth (phases, gains, drive;
% FIT_TRANSDUCER_MODEL with PER_DEPTH_FIT_OPTIONS) on a geometric global model,
% interpolated between measured depths. Quality: held-out interpolation (each
% interior depth predicted from its neighbours) vs geometric, in meta.quality.
% LIBRARY_MODE, if set, must be 'per_depth'.
if exist('LIBRARY_MODE', 'var') && ~strcmp(LIBRARY_MODE, 'per_depth')
    error('build_transducer_library:staticRemoved', ...
        'LIBRARY_MODE ''%s'' has been removed; only per-depth refinements are built.', LIBRARY_MODE);
end

% Per-depth fits (as CALIBRATION_GUI / REFINE_TRANSDUCER_LIBRARY): phase and gain,
% intensity-domain loss, z shift fixed at 0 (a fitted shift is not applied at run
% time). Settings tuned on held-out depths.
if ~exist('PER_DEPTH_FIT_OPTIONS', 'var'); PER_DEPTH_FIT_OPTIONS = calibration_fit_options('rayleigh_phase_gain'); end

% Within-ring velocity profile (FIT_RING_PROFILE): every ring split into
% sub-rings with weights shared by all depths, fitted before the per-depth
% refinements (which then include the near field from the profile's
% NearFieldSkipMm). 'none' (default) | 'fit' (every empirical combo) |
% cellstr of combo names to fit.
if ~exist('RING_PROFILE', 'var');         RING_PROFILE = 'none'; end
if ~exist('RING_PROFILE_OPTIONS', 'var'); RING_PROFILE_OPTIONS = struct('NSub', 4, 'Smooth', 0.3, 'NearFieldSkipMm', 4); end

% ── END CONFIGURATION ─────────────────────────────────────────────────────────

%% ── Load equipment config ────────────────────────────────────────────────────
equip_param = load_equipment_config(equip_path);
all_serials = fieldnames(equip_param.trans);

%% ── Build combo list ─────────────────────────────────────────────────────────
combo_list = {};   % {combo_name, tran_serial, ds_serial, has_char}

for si = 1:numel(all_serials)
    tran_serial = all_serials{si};
    tran        = equip_param.trans.(tran_serial);
    if isfield(tran, 'status') && strcmp(char(tran.status), 'excluded')
        fprintf('  %s: excluded (status: excluded), skipped\n', tran_serial);
        continue;
    end

    if ~isfield(tran, 'combos') || isempty(tran.combos)
        combo_name = tran_serial;
        if ~isequal(COMBOS_TO_RUN, 'all') && ~any(strcmp(COMBOS_TO_RUN, combo_name))
            continue;
        end
        combo_list(end+1,:) = {combo_name, tran_serial, '', false}; %#ok<AGROW>
        continue;
    end

    ds_keys = fieldnames(tran.combos);
    for di = 1:numel(ds_keys)
        ds_serial  = ds_keys{di};
        combo_name = [tran_serial '_' ds_serial];

        if ~isequal(COMBOS_TO_RUN, 'all') && ~any(strcmp(COMBOS_TO_RUN, combo_name))
            continue;
        end

        has_char = isfield(tran.combos.(ds_serial), 'char_data_path') && ...
                   ~isempty(tran.combos.(ds_serial).char_data_path);
        combo_list(end+1,:) = {combo_name, tran_serial, ds_serial, has_char}; %#ok<AGROW>
    end
end

if isempty(combo_list)
    fprintf('No combos matched COMBOS_TO_RUN — nothing to do.\n');
    return;
end

n_total     = size(combo_list, 1);
n_empirical = sum([combo_list{:,4}]);
n_geometric = n_total - n_empirical;

fprintf('\n=== build_transducer_library: %d combo(s) — %d empirical, %d geometric ===\n\n', ...
    n_total, n_empirical, n_geometric);

%% ── Process each combo ───────────────────────────────────────────────────────
results = struct('combo', {}, 'level', {}, 'status', {}, 'msg', {});

for ci = 1:n_total
    combo_name  = combo_list{ci,1};
    tran_serial = combo_list{ci,2};
    ds_serial   = combo_list{ci,3};
    has_char    = combo_list{ci,4};
    tran        = equip_param.trans.(tran_serial);

    level = 'geometric';
    if has_char; level = 'empirical'; end

    fprintf('[%d/%d] %s  (%s)\n', ci, n_total, combo_name, level);
    rec = struct('combo', combo_name, 'level', level, 'status', '', 'msg', '');

    try
        if strcmp(level, 'geometric')
            build_geometric_entry(tran_serial, ds_serial, library_path, tran);
        else
            combo_entry = tran.combos.(ds_serial);
            rp_opts = [];
            if isequal(RING_PROFILE, 'fit') || (iscell(RING_PROFILE) && any(strcmp(RING_PROFILE, combo_name)))
                rp_opts = RING_PROFILE_OPTIONS;
            end
            build_per_depth_entry(combo_name, tran_serial, ds_serial, tran, combo_entry, ...
                CHAR_DATA_ROOT, DESIRED_INTENSITIES(1), PER_DEPTH_FIT_OPTIONS, library_path, rp_opts);
        end
        rec.status = 'ok';
        fprintf('  [OK]\n');
    catch ME
        rec.status = 'error';
        rec.msg    = ME.message;
        fprintf('  [ERROR] %s  (%s line %d)\n', ME.message, ME.stack(1).file, ME.stack(1).line);
    end

    results(end+1) = rec; %#ok<AGROW>
end

%% ── Post-build validation ────────────────────────────────────────────────────
fprintf('\n--- Validating library ---\n');
try
    validate_transducer_library(library_path, equip_path);
catch ME
    fprintf('[WARN] Validation error: %s\n', ME.message);
end

%% ── Summary ──────────────────────────────────────────────────────────────────
n_ok  = sum(strcmp({results.status}, 'ok'));
n_err = sum(strcmp({results.status}, 'error'));

fprintf('\n=== Done: %d ok, %d error(s). Library at: %s ===\n', n_ok, n_err, library_path);
if n_err > 0
    fprintf('\nFailed combos:\n');
    failed = results(strcmp({results.status}, 'error'));
    for fi = 1:numel(failed)
        fprintf('  %-50s %s\n', failed(fi).combo, failed(fi).msg);
    end
end

%% ═══════════════════════ Local helpers ════════════════════════════════════════

function build_geometric_entry(tran_serial, ds_serial, library_path, tran)
% Write a geometry-only library entry (no empirical phase correction).
%
% Stores an empty global_model.elem_phase_correction_deg so that
% load_transducer_from_library reconstructs pure geometric steering phases
% from ring geometry at the requested depth. No amplitude calibration is
% available without characterisation data.

    if isempty(ds_serial)
        combo_key = tran_serial;
    else
        combo_key = [tran_serial '_' ds_serial];
    end

    out_file = library_entry_file(library_path, combo_key, tran_serial, ds_serial);
    if ~isfolder(fileparts(out_file)); mkdir(fileparts(out_file)); end
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    else
        lib = struct();
    end

    lib.meta.tran_serial = tran_serial;
    if ~isempty(ds_serial); lib.meta.driving_system = ds_serial; end
    if ~isfield(lib, 'meta') || ~isfield(lib.meta, 'created_at') || isempty(lib.meta.created_at)
        lib.meta.created_at = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    end
    lib.meta.last_updated     = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    lib.meta.calibration_type = 'geometric';
    lib.meta.geometry_fingerprint = equipment_geometry_fingerprint(tran);
    [lib.meta.prestus_version, lib.meta.prestus_commit] = library_provenance();

    % Drop any legacy per-depth phase tables; phases are computed at load time.
    if isfield(lib, 'global_model')
        for legacy = {'depths_ep_mm', 'elem_phase_deg'}
            if isfield(lib.global_model, legacy{1})
                lib.global_model = rmfield(lib.global_model, legacy{1});
            end
        end
    end
    lib.global_model.elem_phase_correction_deg = [];

    yaml.dumpFile(out_file, lib);
    fprintf('  Written (geometric): %s\n', out_file);
end


function build_per_depth_entry(combo_name, tran_serial, ds_serial, tran, combo_entry, ...
        char_data_root, ref_intensity, fit_options, library_path, rp_opts)
% One refinement per measured depth on a geometric global model. The global
% model (zero offsets, unit gains, measured drive per depth) applies outside
% the measured span; within it, LOAD_TRANSDUCER_FROM_LIBRARY interpolates the
% per-depth phases, gains and drives.

    LIBRARY_MODEL_VERSION = 3;   % 3: per-depth refinements with drive and gains

    char_rel = combo_entry.char_data_path;
    if isempty(char_data_root); char_path = char_rel; else; char_path = fullfile(char_data_root, char_rel); end
    if ~isfile(char_path)
        error('build_transducer_library:missingCharData', ...
            ['Characterisation data not found: %s\n' ...
             'Set CHAR_DATA_ROOT (paths in the equipment YAML are relative), or remove ' ...
             'char_data_path to build a geometric entry on purpose.'], char_path);
    end
    prof = read_axial_profiles(char_path);
    n_d  = numel(prof.depths_ep_mm);
    fprintf('  Depths: %s mm\n', mat2str(prof.depths_ep_mm));

    % Optional ring profile: the whole entry (global model, refinements) then
    % refers to the sub-ring geometry
    rp = [];
    if ~isempty(rp_opts)
        rp = fit_ring_profile(tran, prof, rp_opts);
        fprintf('  Ring profile (%d sub-rings): median R² %.3f → %.3f (z >= %g mm)\n', rp.n_sub, ...
            rp.r2_uniform_median, rp.r2_profile_median, rp.near_field_skip_mm);
        tran = ring_profile_transducer(tran, rp);
        if ~isfield(fit_options, 'NearFieldSkipMm') || isequal(fit_options.NearFieldSkipMm, 10)
            fit_options.NearFieldSkipMm = rp.near_field_skip_mm;
        end
    end

    ev  = struct('EvaluateOnly', true);
    for f = {'ForwardModel', 'Domain', 'NearFieldSkipMm'}
        if isfield(fit_options, f{1}); ev.(f{1}) = fit_options.(f{1}); end
    end
    geo  = fit_transducer_model(tran, prof, ev);
    fits = calibration_fit_depths(tran, prof, fit_options);
    r2_fit = cellfun(@(f) f.metrics.r2, fits);
    dz     = cellfun(@(f) f.z_shift_mm, fits);

    % Held-out interpolation: interior depth j predicted from j-1 and j+1
    % (phases on the unit circle, gains linearly), without z shift as at run time
    r2_int = nan(1, n_d);
    for j = 2:n_d - 1
        a = fits{j - 1}; b = fits{j + 1};
        t = (prof.depths_ep_mm(j) - a.fit_depths_ep_mm) / (b.fit_depths_ep_mm - a.fit_depths_ep_mm);
        u = unwrap([a.elem_phase_correction_deg; b.elem_phase_correction_deg] * pi/180, [], 1);
        pr = struct('elem_phase_correction_deg', mod(((1 - t) * u(1, :) + t * u(2, :)) * 180/pi, 360), ...
                    'elem_amp_weight', (1 - t) * a.elem_amp_weight + t * b.elem_amp_weight);
        pj = prof; pj.depths_ep_mm = prof.depths_ep_mm(j); pj.I = prof.I(:, j);
        sc = ev; sc.Prior = pr;
        r2_int(j) = fit_transducer_model(tran, pj, sc).metrics.r2;
    end
    r2_geo = geo.metrics.r2(:)';
    q = struct('mode', 'per_depth', 'n_depths', n_d, ...
        'geometric_r2_median', median(r2_geo), 'in_sample_r2_median', median(r2_fit), ...
        'in_sample_r2_min', min(r2_fit), ...
        'interp_r2_median', median(r2_int, 'omitnan'), ...
        'interp_worst_delta_r2', min(r2_int - r2_geo, [], 'omitnan'), ...
        'z_shift_mm', dz, 'flags', {flags_of(fits)});
    fprintf(['  Per-depth: in-sample median %.3f (min %.3f) | held-out interpolation %.3f ' ...
        '(worst vs geometric %+.3f) | geometric %.3f | dz %s mm\n'], q.in_sample_r2_median, ...
        q.in_sample_r2_min, q.interp_r2_median, q.interp_worst_delta_r2, q.geometric_r2_median, ...
        mat2str(round(dz, 1)));

    out_file = library_entry_file(library_path, combo_name, tran_serial, ds_serial);
    if ~isfolder(fileparts(out_file)); mkdir(fileparts(out_file)); end
    if isfile(out_file); lib = yaml.loadFile(out_file, 'ConvertToArray', true); else; lib = struct(); end
    for f = {'global_model', 'refinement'}
        if isfield(lib, f{1}); lib = rmfield(lib, f{1}); end
    end
    lib.meta.tran_serial = tran_serial;
    if ~isempty(ds_serial); lib.meta.driving_system = ds_serial; end
    if ~isfield(lib.meta, 'created_at') || isempty(lib.meta.created_at)
        lib.meta.created_at = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    end
    lib.meta.last_updated          = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    lib.meta.library_model_version = LIBRARY_MODEL_VERSION;
    lib.meta.calibration_type      = 'empirical';
    lib.meta.geometry_fingerprint  = equipment_geometry_fingerprint(tran);
    lib.meta.library_mode          = 'per_depth';
    lib.meta.char_data             = char_rel;
    lib.meta.quality               = q;
    [lib.meta.prestus_version, lib.meta.prestus_commit] = library_provenance();

    lib = library_entry_from_model(lib, geo, ref_intensity);   % geometric global model, measured drive
    lib.global_model.elem_phase_correction_deg = [];           % pure geometric steering outside the refinements
    lib.global_model = rmfield(lib.global_model, {'elem_amp_weight', 'free_params', 'fit_z_shift_mm'});
    lib = library_refinement_from_fits(lib, fits, ref_intensity, [], char_rel);
    if ~isempty(rp)
        lib.global_model.ring_profile = rp;
        lib.meta.ring_profile = true;
    elseif isfield(lib.meta, 'ring_profile')
        lib.meta = rmfield(lib.meta, 'ring_profile');
    end
    yaml.dumpFile(out_file, lib);
    fprintf('  Written: %s\n', out_file);
end

function fl = flags_of(fits)
    fl = cellfun(@(f) reshape(cellstr(f.flags), 1, []), fits, 'UniformOutput', false);
    fl = unique([fl{:}]);
end
