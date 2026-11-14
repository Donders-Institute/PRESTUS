function [lib, fits, report, prof, diag] = calibrate_library_entry(lib, tran, meas, options)
% CALIBRATE_LIBRARY_ENTRY  Calibrate a transducer library entry to hydrophone profiles at chosen depths
%
% The calibration core behind REFINE_TRANSDUCER_LIBRARY and
% CALIBRATION_TRANSDUCER. It does not read or write files: the caller loads
% the base entry and the measured profiles and decides where the result goes.
%
%   1. Profiles: a requested depth within 0.5 mm of a measured one uses that
%      measurement; any other is interpolated from its two neighbouring
%      measurements (INTERPOLATE_AXIAL_PROFILES), or extrapolated within
%      MaxExtrapolationMm. Unavailable depths are skipped with a warning.
%      An interpolated depth whose neighbours are more than MaxGapMm apart is
%      flagged: neither profile nor offset interpolation is reliable there.
%   2. Fit (Mode):
%        'per_depth' - one fit per depth (CALIBRATION_FIT_DEPTHS) with the
%                      preset's options, regularised towards the entry's
%                      static model; stored as refinements
%                      (LIBRARY_REFINEMENT_FROM_FITS).
%        'geometric' - no fit: geometric steering (plus any static offset)
%                      scored per depth and stored as refinements with the
%                      measured drive.
%   Optional within-ring velocity profile (RingProfile, FIT_RING_PROFILE):
%      every ring is split into sub-rings with fixed weights shared by all
%      depths, fitted on all measured depths before the per-depth fits; the
%      per-depth fits then calibrate the channels of that geometry.
%   3. Optional free-water k-Wave drive calibration of each refined depth: the
%      stored drive is rescaled to the simulated ISPPA (KWAVE_CALIBRATE_DRIVE).
%   4. Report: the entry before and after, scored on the same profiles.
%
% Use as:
%   [lib, fits, report] = calibrate_library_entry(lib, tran, meas)
%   [lib, fits, report, prof, diag] = calibrate_library_entry(lib, tran, meas, options)
%
% Input:
%   lib     - base library struct (yaml.loadFile) or struct('meta', ...) for a new entry
%   tran    - transducer struct from load_equipment_config (annular)
%   meas    - measured profiles (READ_AXIAL_PROFILES): .depths_ep_mm, .z_ep_mm, .I
%   options - struct, all optional:
%     .Depths             - depths to calibrate [mm from EP] (default: all measured)
%     .Mode               - 'per_depth' (default) | 'geometric'
%     .Preset             - CALIBRATION_FIT_OPTIONS preset, 'rayleigh_phase_gain' |
%                           'oneil_phase_only' (default 'rayleigh_phase_gain')
%     .FitOptions         - FIT_TRANSDUCER_MODEL options on top of the preset
%     .Intensity          - reference ISPPA for the stored drives [W/cm²] (default 30)
%     .BlendMm            - refinement.blend_mm (default: keep existing / 5)
%     .Source             - provenance string stored with the refinements
%     .InterpolateOptions - options for INTERPOLATE_AXIAL_PROFILES (Align,
%                           MaxExtrapolationMm, ...)
%     .MaxGapMm           - warn above this neighbour spacing [mm] (default 10)
%     .KWaveDriveCalibration - run KWAVE_CALIBRATE_DRIVE (default false;
%                           former name KWave, still read)
%     .KWaveDriveOptions  - options for KWAVE_CALIBRATE_DRIVE (former name
%                           KWaveOptions, still read)
%     .EquipParam         - struct from load_equipment_config (passed to the k-Wave drive calibration)
%     .RingProfile        - within-ring velocity profile: 'auto' (default: keep the
%                           entry's profile, if any) | 'fit' (FIT_RING_PROFILE on all
%                           measured depths) | 'none' (uniform rings) | a profile struct.
%                           A changed profile changes the geometry: the entry's
%                           refinements and static offsets are dropped (warning)
%     .RingProfileOptions - options for FIT_RING_PROFILE (NSub, Smooth, NearFieldSkipMm, ...)
%     .ProgressFcn        - @(k, n, depth_ep_mm, profile_type) before each fit
%
% Output:
%   lib    - calibrated entry (meta is left to the caller)
%   fits   - {1 x N} FIT_TRANSDUCER_MODEL results, one per depth
%   report - table per calibrated depth: depth_ep_mm, profile_type, r2_library,
%            r2_refined, abs_focal_err_library_mm, abs_focal_err_refined_mm,
%            nearfield_excess_refined, source_gap_mm (NaN unless interpolated)
%   prof   - the profiles the entry was calibrated on (with profile_type)
%   diag   - .ring_profile (second output of FIT_RING_PROFILE, or []),
%            .before (EVALUATE_LIBRARY_ENTRY of the base entry), .z_after,
%            .I_after (modelled profiles after calibration), .kwave (k-Wave
%            report table or []), .fit_options (options actually used)
%
% See also: REFINE_TRANSDUCER_LIBRARY, CALIBRATION_TRANSDUCER, CALIBRATION_FIT_OPTIONS,
%           CALIBRATION_FIT_DEPTHS, INTERPOLATE_AXIAL_PROFILES, KWAVE_CALIBRATE_DRIVE

    if nargin < 4; options = struct(); end
    mode      = opt(options, 'Mode', 'per_depth');
    preset    = opt(options, 'Preset', 'rayleigh_phase_gain');
    fit_opts  = opt(options, 'FitOptions', struct());
    intensity = opt(options, 'Intensity', 30);
    blend_mm  = opt(options, 'BlendMm', []);
    source    = opt(options, 'Source', '');
    max_gap   = opt(options, 'MaxGapMm', 10);
    if strcmp(mode, 'static')
        error('calibrate_library_entry:staticRemoved', ...
            'Mode ''static'' (joint depth-independent model) has been removed; use ''per_depth''.');
    end
    assert(ismember(mode, {'per_depth', 'geometric'}), 'calibrate_library_entry:mode', ...
        'calibrate_library_entry: Mode must be ''per_depth'' or ''geometric''.');

    %% Profiles at the requested depths ------------------------------------
    use  = reshape(opt(options, 'Depths', meas.depths_ep_mm), 1, []);
    near = arrayfun(@(u) find(abs(meas.depths_ep_mm - u) < 0.5, 1), use, 'UniformOutput', false);
    is_m = ~cellfun(@isempty, near);
    want = unique([meas.depths_ep_mm([near{is_m}]), use(~is_m)]);
    io   = opt(options, 'InterpolateOptions', struct());
    if ~isfield(io, 'NearFieldSkipMm'); io.NearFieldSkipMm = opt(fit_opts, 'NearFieldSkipMm', 10); end
    prof = interpolate_axial_profiles(meas, want, io);
    bad  = strcmp(prof.profile_type, 'unavailable');
    if any(bad)
        warning('calibrate_library_entry:outsideProfiles', ...
            ['Depth(s) %s mm lie outside the measured range [%.1f, %.1f] mm (MaxExtrapolationMm = %g) ' ...
             'and are skipped.'], mat2str(prof.depths_ep_mm(bad)), min(meas.depths_ep_mm), ...
            max(meas.depths_ep_mm), opt(io, 'MaxExtrapolationMm', 0));
        prof = profile_subset(prof, ~bad);
    end
    if isempty(prof.depths_ep_mm)
        error('calibrate_library_entry:noProfiles', ...
            'None of the requested depths %s mm can be taken from the measured profiles (available: %s).', ...
            mat2str(use), mat2str(meas.depths_ep_mm));
    end
    derived = ~strcmp(prof.profile_type, 'measured');
    if any(derived)
        fprintf('calibrate_library_entry: %d depth(s) calibrated on derived profiles: %s\n', nnz(derived), ...
            strjoin(compose('%.2f mm (%s)', prof.depths_ep_mm(derived)', string(prof.profile_type(derived))'), ', '));
    end
    gap = nan(1, numel(prof.depths_ep_mm));
    interp = strcmp(prof.profile_type, 'interpolated');
    gap(interp) = abs(diff(prof.source_depths(interp, :), 1, 2))';
    if any(gap > max_gap)
        warning('calibrate_library_entry:largeGap', ...
            ['Depth(s) %s mm are interpolated between measurements %s mm apart (> %g mm). ' ...
             'The focal behaviour can change within such a gap; measure closer to these depths.'], ...
            mat2str(prof.depths_ep_mm(gap > max_gap)), mat2str(gap(gap > max_gap)), max_gap);
    end
    if isfield(tran, 'min_foc') && isfield(tran, 'max_foc') && ...
            any(prof.depths_ep_mm < tran.min_foc | prof.depths_ep_mm > tran.max_foc)
        warning('calibrate_library_entry:outOfRange', ...
            'Some depths lie outside [min_foc, max_foc] = [%.1f, %.1f] mm.', tran.min_foc, tran.max_foc);
    end

    %% Base entry as a study resolves it (for the report) ------------------
    fo = calibration_fit_options(preset, fit_opts);
    before = evaluate_library_entry(tran, lib, prof, struct('ForwardModel', fo.ForwardModel));

    %% Ring profile -----------------------------------------------------------
    rp_info = [];
    rp_mode = opt(options, 'RingProfile', 'auto');
    if isstruct(rp_mode) || any(strcmp(rp_mode, {'fit', 'none'}))
        if isstruct(rp_mode)
            rp = rp_mode;
        elseif strcmp(rp_mode, 'fit')
            ro = opt(options, 'RingProfileOptions', struct());
            if ~isfield(ro, 'ForwardModel'); ro.ForwardModel = fo.ForwardModel; end
            [rp, rp_info] = fit_ring_profile(tran, meas, ro);
            fprintf('calibrate_library_entry: ring profile (%d sub-rings) median R² %.3f → %.3f (z >= %g mm)\n', ...
                rp.n_sub, rp.r2_uniform_median, rp.r2_profile_median, rp.near_field_skip_mm);
        else
            rp = [];
        end
        lib = library_set_ring_profile(lib, rp);
    elseif ~strcmp(rp_mode, 'auto')
        error('calibrate_library_entry:ringProfile', 'RingProfile must be ''auto'', ''fit'', ''none'' or a struct.');
    end
    tran = library_transducer(tran, lib);
    if isfield(tran.transducer.annular, 'ring_profile') && ~isfield(fit_opts, 'NearFieldSkipMm')
        % The profile is fitted with the near field; calibrate the channels on the same range
        fo.NearFieldSkipMm = tran.transducer.annular.ring_profile.near_field_skip_mm;
    end

    %% Fit -------------------------------------------------------------------
    % Regularise towards the entry's global model (static offset of a legacy
    % entry, else geometric steering). z shift stays 0: the loader does not
    % apply a fitted shift, and phases fitted with one match the hydrophone
    % only after that shift
    if ~isfield(fo, 'Prior') || isempty(fo.Prior); fo.Prior = library_static_model(lib); end
    fo.Prior.z_shift_mm = 0;
    if strcmp(mode, 'geometric')
        fo.EvaluateOnly = true;
    end
    fits = calibration_fit_depths(tran, prof, fo, opt(options, 'ProgressFcn', []));
    lib  = library_refinement_from_fits(lib, fits, intensity, blend_mm, source);

    %% Optional k-Wave drive calibration --------------------------------------
    kw = [];
    if opt(options, 'KWaveDriveCalibration', opt(options, 'KWave', false))
        ko = opt(options, 'KWaveDriveOptions', opt(options, 'KWaveOptions', struct()));
        if isfield(options, 'EquipParam') && ~isfield(ko, 'EquipParam'); ko.EquipParam = options.EquipParam; end
        [lib, kw] = kwave_calibrate_drive(lib, prof, ko);
        disp(kw);
    end

    lib.meta.geometry_fingerprint = equipment_geometry_fingerprint(tran);   % checked on load

    %% Report ------------------------------------------------------------------
    after = vertcat_metrics(fits);
    I_after = cell2mat(cellfun(@(f) f.I_model, fits, 'UniformOutput', false)); z_after = fits{1}.z_ep_mm;
    report = table(prof.depths_ep_mm(:), string(prof.profile_type(:)), before.metrics.r2, after.r2, ...
        abs(before.metrics.focal_err_mm), abs(after.focal_err_mm), after.nearfield_excess, gap(:), ...
        'VariableNames', {'depth_ep_mm', 'profile_type', 'r2_library', 'r2_refined', 'abs_focal_err_library_mm', ...
        'abs_focal_err_refined_mm', 'nearfield_excess_refined', 'source_gap_mm'});
    diag = struct('before', before, 'z_after', z_after, 'I_after', I_after, 'kwave', kw, 'fit_options', fo, ...
        'ring_profile', rp_info);
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
