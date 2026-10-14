function [fit, msg] = calibration_fit(tran, prof, options)
% CALIBRATION_FIT  Interactive calibration fit: per-depth refinements
%
% The "Fit" step of CALIBRATION_GUI and the Calibration page of the PRESTUS
% Nexus. Depths without a profile ('unavailable', e.g. desired depths
% outside the measured range) are skipped. Per-depth mode fits each depth on
% its own (CALIBRATION_FIT_DEPTHS) with phase and gain free; the axial shift
% stays 0 because it is not applied at run time, and the prior (the current
% entry's global model, if any) only pins the unobservable common phase.
% The result carries the library entry Save would write
% (CALIBRATION_FIT_ENTRY) and a record of how it was fitted, which
% CALIBRATION_SAVE_ENTRY stores as meta.fit.
%
% Use as:
%   fit = calibration_fit(tran, prof)
%   [fit, msg] = calibration_fit(tran, prof, options)
%
% Input:
%   tran    - annular transducer struct (load_equipment_config)
%   prof    - profiles to fit: .depths_ep_mm, .z_ep_mm, .I, .profile_type
%             (CALIBRATION_DESIRED_PROFILES / INTERPOLATE_AXIAL_PROFILES);
%             with .n_z_scan only the first n_z_scan z samples are fitted
%   options - struct, all optional:
%     .Mode          - 'per_depth' (default; the former 'static' has been removed)
%     .FreeParams    - subset of {'phase', 'amp'} (default both)
%     .FitOptions    - further FIT_TRANSDUCER_MODEL options (ForwardModel,
%                      NearFieldSkipMm, HuberDelta, NStarts, Domain, advanced ones)
%     .Library       - current library entry (prior and base of the new entry) or []
%     .IsppaWcm2     - reference ISPPA for the stored drives (default 30)
%     .TranSerial    - transducer serial for the entry's meta (default '')
%     .DrivingSystem - driving-system serial ('' = none)
%     .Source        - provenance of the profiles, e.g. the profile CSV
%     .CalibrateAt   - 'measured' (default) | 'desired', recorded in meta.fit
%     .ProfileOptions - INTERPOLATE_AXIAL_PROFILES options used for derived
%                      profiles (Align, MaxExtrapolationMm), recorded in meta.fit
%     .Progress      - @(k, n, depth_ep_mm, profile_type) before each per-depth fit
%     .Tool          - name of the caller, recorded in meta.fit (default 'calibration_fit')
%     .RingProfile   - within-ring velocity profile: 'auto' (default: keep the
%                      Library entry's) | 'fit' (FIT_RING_PROFILE on the measured
%                      profiles) | 'none' | a profile struct
%     .RingProfileOptions - FIT_RING_PROFILE options (NSub, Smooth, NearFieldSkipMm)
%     .Measured      - measured profiles for the ring-profile fit (default: the
%                      measured columns of prof)
%
% Output:
%   fit - .mode ('per_depth'), .model ({1 x N} FIT_TRANSDUCER_MODEL results),
%         .prof (the profiles fitted),
%         .r2 [1 x N], .free_params, .entry (CALIBRATION_FIT_ENTRY), .record
%         (meta.fit), .seconds, .ring_profile (FIT_RING_PROFILE info, or [])
%   msg - one-line summary
%
% See also: CALIBRATION_FIT_DEPTHS, FIT_TRANSDUCER_MODEL, CALIBRATION_FIT_ENTRY,
%           CALIBRATION_SAVE_ENTRY, CALIBRATION_EVALUATE, CALIBRATION_GUI

    if nargin < 3; options = struct(); end
    mode = char(opt(options, 'Mode', 'per_depth'));
    assert(strcmp(mode, 'per_depth'), 'calibration_fit:mode', ...
        'Mode must be ''per_depth'', not ''%s'' (the static joint model has been removed).', mode);
    free = cellstr(opt(options, 'FreeParams', {'phase', 'amp'}));
    fo   = opt(options, 'FitOptions', struct());
    lib  = opt(options, 'Library', []);
    msg  = '';

    if ~isfield(prof, 'profile_type')
        prof.profile_type = repmat({'measured'}, 1, numel(prof.depths_ep_mm));
    end
    if isfield(prof, 'n_z_scan')            % drop the model-only range of desired profiles
        prof.z_ep_mm = prof.z_ep_mm(1:prof.n_z_scan); prof.I = prof.I(1:prof.n_z_scan, :);
        prof = rmfield(prof, 'n_z_scan');
    end
    none = strcmp(prof.profile_type, 'unavailable') | ~any(isfinite(prof.I), 1);
    if any(none)
        msg = sprintf('No profile at %s mm: skipped. ', mat2str(prof.depths_ep_mm(none)));
        prof = profile_subset(prof, ~none);
    end
    assert(~isempty(prof.depths_ep_mm), 'calibration_fit:noProfiles', ...
        'None of the depths has a profile to calibrate on.');
    t0 = tic;
    % Within-ring velocity profile: fitted on the measured profiles, shared by all depths
    rp_info = []; rp = [];
    rp_mode = opt(options, 'RingProfile', 'auto');
    if isstruct(rp_mode) || any(strcmp(rp_mode, {'fit', 'none'}))
        if isstruct(rp_mode)
            rp = rp_mode;
        elseif strcmp(rp_mode, 'fit')
            meas = opt(options, 'Measured', profile_subset(prof, strcmp(prof.profile_type, 'measured')));
            assert(numel(meas.depths_ep_mm) >= 2, 'calibration_fit:ringProfileDepths', ...
                'A ring profile needs measured profiles at two or more depths.');
            ro = opt(options, 'RingProfileOptions', struct());
            if ~isfield(ro, 'ForwardModel') && isfield(fo, 'ForwardModel'); ro.ForwardModel = fo.ForwardModel; end
            [rp, rp_info] = fit_ring_profile(tran, meas, ro);
        end
        if ~isempty(lib) || ~isempty(rp); lib = library_set_ring_profile(lib, rp); end
    elseif ~strcmp(rp_mode, 'auto')
        error('calibration_fit:ringProfile', 'RingProfile must be ''auto'', ''fit'', ''none'' or a struct.');
    end
    tran = library_transducer(tran, lib);
    if isfield(tran.transducer.annular, 'ring_profile')
        rp = tran.transducer.annular.ring_profile;
        if ~isfield(fo, 'NearFieldSkipMm'); fo.NearFieldSkipMm = rp.near_field_skip_mm; end
    end
    % z shift stays at 0 (also for the prior). The prior only pins the
    % unobservable common phase (default Lambda): a stronger pull stops
    % the fit from reaching split-focus / near-field-peak solutions.
    fo.FreeParams = intersect(free, {'phase', 'amp'}, 'stable');
    if ~isempty(lib); fo.Prior = library_static_model(lib); end
    fo.Prior.z_shift_mm = 0;
    model = calibration_fit_depths(tran, prof, fo, opt(options, 'Progress', []));
    r2 = cellfun(@(m) m.metrics.r2, model);

    fit = struct('mode', mode, 'model', {model}, 'prof', prof, 'r2', r2, ...
        'free_params', {fo.FreeParams}, 'seconds', toc(t0), 'ring_profile', rp_info);
    isppa = opt(options, 'IsppaWcm2', 30);
    source = char(opt(options, 'Source', ''));
    fit.entry = calibration_fit_entry(lib, fit, char(opt(options, 'TranSerial', '')), ...
        char(opt(options, 'DrivingSystem', '')), isppa, source);

    fo_rec = rmfield_if(fo, 'Prior');
    fit.record = struct('tool', char(opt(options, 'Tool', 'calibration_fit')), 'mode', mode, ...
        'calibrate_at', char(opt(options, 'CalibrateAt', 'measured')), ...
        'depths_ep_mm', prof.depths_ep_mm(:)', 'profile_type', {prof.profile_type}, ...
        'profile_csv', source, 'isppa_wcm2', isppa, 'fit_options', fo_rec, ...
        'profile_options', opt(options, 'ProfileOptions', struct()));
    if ~isempty(rp)
        fit.record.ring_profile = struct('n_sub', rp.n_sub, 'smooth', rp.smooth, ...
            'near_field_skip_mm', rp.near_field_skip_mm, 'refitted', ~isempty(rp_info));
    end

    pt = prof.profile_type;
    if ~isempty(rp_info)
        msg = [msg sprintf('Ring profile (%d sub-rings): median R² %.3f → %.3f. ', rp.n_sub, ...
            rp.r2_uniform_median, rp.r2_profile_median)];
    end
    msg = [msg sprintf(['Fit (%s, %s) on %d profiles (%d measured, %d interpolated, %d extrapolated) ' ...
        'in %.1f s: median R² %.3f, min %.3f'], mode, strjoin(fo.FreeParams, '+'), numel(pt), ...
        nnz(strcmp(pt, 'measured')), nnz(strcmp(pt, 'interpolated')), nnz(strcmp(pt, 'extrapolated')), ...
        fit.seconds, median(r2), min(r2))];
end

function p = profile_subset(p, keep)
    p.depths_ep_mm = p.depths_ep_mm(keep);
    p.I            = p.I(:, keep);
    p.profile_type = p.profile_type(keep);
    for f = {'source_depths', 'weights', 'shift_mm'}
        if isfield(p, f{1}); p.(f{1}) = p.(f{1})(keep, :); end
    end
end

function s = rmfield_if(s, f)
    if isfield(s, f); s = rmfield(s, f); end
end

function v = opt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
