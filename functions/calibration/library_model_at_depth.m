function s = library_model_at_depth(lib, tran, focal_distance_ep, desired_intensity, source_label)
% LIBRARY_MODEL_AT_DEPTH  Phases, gains and drive of a library entry at any focal depth
%
% The depth logic of LOAD_TRANSDUCER_FROM_LIBRARY, on an in-memory entry
% (no file lookup), so callers can evaluate an entry that is not yet saved
% (e.g. a fit in CALIBRATION_GUI) exactly as a study would resolve it.
%
%   phase(d) = set_real_phases(d; c_steer) + offset(d)
%
%   offset(d) = static offset (global_model.elem_phase_correction_deg)
%               blended towards the refinement offsets
%               (refinement.elem_phase_correction_deg at refinement.depths_ep_mm)
%
% Within the span of refined depths the refinement offsets are interpolated
% per element on the unit circle (weight w = 1). Outside it the nearest
% refinement fades back to the static offset over refinement.blend_mm
% (default 5 mm): w = max(0, 1 - distance / blend_mm). Drives and per-ring
% gains of the refinements are interpolated in depth (clamped outside the
% span) and blended with the global ones by the same weight.
%
% Legacy per-depth forms are read as refinements: global_model.depths_ep_mm +
% elem_phase_deg (absolute phases) and calibration.focal_depths.f* entries
% (linear/monotonic parametric phases).
%
% Use as:
%   s = library_model_at_depth(lib, tran, focal_distance_ep, desired_intensity)
%   s = library_model_at_depth(lib, tran, focal_distance_ep, desired_intensity, source_label)
%
% Input:
%   lib               - library struct (yaml.loadFile(..., 'ConvertToArray', true))
%   tran              - transducer struct from load_equipment_config (annular)
%   focal_distance_ep - focal depth from the exit plane [mm]
%   desired_intensity - target ISPPA [W/cm²]; NaN = drive at the stored reference
%   source_label      - (optional) entry name for error messages, e.g. its YAML path
%
% Output (s):
%   .offset_deg            - [1 x N_elem] offset relative to geometric steering [deg]
%   .elem_phase_deg/_rad   - [1 x N_elem] absolute element phases [0, 360) / [0, 2π)
%   .steering_sound_speed  - sound speed of the steering law [m/s]
%   .drive_pa              - mean source pressure for desired_intensity [Pa]
%                            ([] = no amplitude calibration)
%   .elem_amp_weight       - [1 x N_elem] per-ring gains (mean 1)
%   .elem_p_amp            - [1 x N_elem] round(drive_pa * elem_amp_weight) ([] = none)
%   .phase_model           - 'geometric' | 'static' | 'refined' | 'blend'
%   .refinement_weight     - w, 0 (static) … 1 (refined)
%   .amplitude_source      - 'global' | 'refinement' | 'blend'
%   .amp_kwave_calibrated    - true iff the resolved drive comes wholly from a
%                            refinement whose amplitude was corrected against
%                            a free-water k-Wave simulation (KWAVE_CALIBRATE_DRIVE,
%                            refinement.amp_source == 'kwave'). False for the
%                            analytical-only global model, an unchecked
%                            refinement, or any blend between them.
%   .amp_kwave_nearest_depth_mm/.amp_kwave_r2/.amp_kwave_focal_err_mm
%                          - nearest refinement depth used for the check above
%                            and its k-Wave-vs-measured agreement, for context
%                            even when amp_kwave_calibrated is false (NaN if no
%                            refinement drive exists at all)
%
% See also: LOAD_TRANSDUCER_FROM_LIBRARY, LIBRARY_STATIC_MODEL, LIBRARY_ADD_REFINEMENTS,
%           KWAVE_CALIBRATE_DRIVE

    if nargin < 5 || isempty(source_label); source_label = '(in-memory entry)'; end
    n_elem = tran.transducer.annular.elem_n;
    tran.n_elem = n_elem;   % one phase per ring, as LOAD_TRANSDUCER_FROM_LIBRARY

    gm = struct();
    if isfield(lib, 'global_model') && isstruct(lib.global_model)
        gm = lib.global_model;
    end

    % Steering law: geometric delays at the fitted steering sound speed
    c_steer = 1500;
    if isfield(gm, 'steering_sound_speed') && ~isempty(gm.steering_sound_speed)
        c_steer = double(gm.steering_sound_speed);
    end
    geo_params = struct();
    geo_params.medium_properties.water.sound_speed = c_steer;
    geo_params.transducer        = tran.transducer;
    geo_params.transducer.n_elem = n_elem;
    geo_at = @(d) reshape(set_real_phases(tran, d, geo_params), 1, []);

    %% Static offset -------------------------------------------------------
    static_deg = zeros(1, n_elem);
    has_static = false;
    if isfield(gm, 'elem_phase_correction_deg') && ~isempty(gm.elem_phase_correction_deg)
        static_deg = check_length(gm.elem_phase_correction_deg, n_elem, ...
            'global_model.elem_phase_correction_deg', source_label);
        has_static = true;
    elseif isfield(tran, 'elem_phase_correction') && ~isempty(tran.elem_phase_correction)
        % Legacy: correction stored in the equipment YAML (save_elem_correction)
        static_deg = check_length(tran.elem_phase_correction.deg, n_elem, ...
            'equipment elem_phase_correction.deg', source_label);
        has_static = true;
    end

    %% Refinement offsets --------------------------------------------------
    [ref_depths, ref_deg, ref_amp, blend_mm, ref_pamp, ref_isppa, ref_gain, ref_amp_src, ref_kw_r2, ref_kw_fe] = ...
        collect_refinements(lib, n_elem, geo_at, source_label);
    [ref_offset_deg, w] = refinement_at(ref_depths, ref_deg, focal_distance_ep, blend_mm);

    % Circular blend of static and refinement offsets
    if w > 0
        z = (1 - w) * exp(1i * static_deg * pi/180) + w * exp(1i * ref_offset_deg * pi/180);
        offset_deg = angle(z) * 180/pi;
    else
        offset_deg = static_deg;
    end
    elem_phase_rad = mod((geo_at(focal_distance_ep) + offset_deg) * pi/180, 2*pi);

    if w >= 1
        phase_model = 'refined';
    elseif w > 0
        phase_model = 'blend';
    elseif has_static
        phase_model = 'static';
    else
        phase_model = 'geometric';
    end

    %% Amplitude -----------------------------------------------------------
    amp_map = [];
    if isfield(gm, 'amplitude_scaling') && ~isempty(gm.amplitude_scaling)
        amp_map = gm.amplitude_scaling;
    elseif ~isempty(ref_amp)
        % Legacy per-depth entries only: amplitude of the nearest refinement
        [~, k]  = min(abs(ref_depths - focal_distance_ep));
        amp_map = ref_amp{k};
    end
    drive = [];
    if isfield(gm, 'amplitude_by_depth') && isstruct(gm.amplitude_by_depth) && ...
            ~isempty(gm.amplitude_by_depth.depths_ep_mm)
        % Water-calibrated drive per depth (hydrophone profiles at matched
        % ISPPA): interpolate, clamp outside the range, scale as sqrt(I)
        abd = gm.amplitude_by_depth;
        dd  = double(abd.depths_ep_mm(:)'); aa = double(abd.elem_p_amp(:)');
        dq  = min(max(focal_distance_ep, dd(1)), dd(end));
        if numel(dd) == 1; a_ref = aa; else; a_ref = interp1(dd, aa, dq, 'linear'); end
        drive = a_ref;
        if ~isnan(desired_intensity)
            drive = round(a_ref * sqrt(desired_intensity / double(abd.isppa_wcm2)));
        end
    elseif ~isempty(amp_map)
        drive = lookup_amplitude(amp_map, desired_intensity);
    end
    % Per-depth drive of the refinements (fitted with the refined phases, or
    % k-Wave-corrected): interpolated like the phase offsets and blended
    % with the global amplitude by the same weight w. Without a global
    % amplitude the nearest refined drive is used at any depth.
    amp_source = 'global';
    amp_kwave_calibrated = false;
    amp_kwave_nearest_depth_mm = NaN;
    amp_kwave_r2 = NaN;
    amp_kwave_focal_err_mm = NaN;
    ok = isfinite(ref_pamp);
    if any(ok)
        dk = ref_depths(ok); ak = ref_pamp(ok); src_k = ref_amp_src(ok);
        r2_k = ref_kw_r2(ok); fe_k = ref_kw_fe(ok);
        dq = min(max(focal_distance_ep, dk(1)), dk(end));
        if numel(dk) == 1; a_r = ak; else; a_r = interp1(dk, ak, dq, 'linear'); end
        if ~isnan(desired_intensity); a_r = a_r * sqrt(desired_intensity / ref_isppa); end
        % Nearest refinement point actually driving this depth's amplitude, for
        % k-Wave drive-calibration provenance (categorical, so nearest-neighbour
        % rather than interpolated like the numeric drive above).
        [~, nearest_i] = min(abs(dk - dq));
        nearest_kwave = strcmpi(src_k{nearest_i}, 'kwave');
        amp_kwave_nearest_depth_mm = dk(nearest_i);
        amp_kwave_r2 = r2_k(nearest_i);
        amp_kwave_focal_err_mm = fe_k(nearest_i);
        if isempty(drive)
            drive = round(a_r); amp_source = 'refinement';
            amp_kwave_calibrated = nearest_kwave;
        elseif w > 0
            drive = round((1 - w) * drive + w * a_r);
            amp_source = ternary(w >= 1, 'refinement', 'blend');
            amp_kwave_calibrated = (w >= 1) && nearest_kwave;
        end
    end

    gains = ones(1, n_elem);
    if isfield(gm, 'elem_amp_weight') && ~isempty(gm.elem_amp_weight)
        gains = check_length(gm.elem_amp_weight, n_elem, 'global_model.elem_amp_weight', source_label);
    end
    % Per-depth gains of the refinements, interpolated and blended like the drive
    okg = all(isfinite(ref_gain), 2)';
    if any(okg) && w > 0
        dk = ref_depths(okg); gk = ref_gain(okg, :);
        dq = min(max(focal_distance_ep, dk(1)), dk(end));
        if numel(dk) == 1; g_r = gk; else; g_r = interp1(dk(:), gk, dq, 'linear'); end
        gains = (1 - w) * gains + w * g_r;
        gains = gains / mean(gains);
    end

    s.offset_deg           = offset_deg;
    s.elem_phase_rad       = elem_phase_rad;
    s.elem_phase_deg       = elem_phase_rad / pi * 180;
    s.steering_sound_speed = c_steer;
    s.drive_pa             = drive;
    s.elem_amp_weight      = gains;
    s.elem_p_amp           = [];
    if ~isempty(drive); s.elem_p_amp = round(drive * gains); end
    s.phase_model          = phase_model;
    s.refinement_weight    = w;
    s.amplitude_source     = amp_source;
    s.amp_kwave_calibrated   = amp_kwave_calibrated;
    s.amp_kwave_nearest_depth_mm = amp_kwave_nearest_depth_mm;
    s.amp_kwave_r2         = amp_kwave_r2;
    s.amp_kwave_focal_err_mm = amp_kwave_focal_err_mm;
end

%% ======================================================================
function v = check_length(v, n_elem, what, source_label)
    v = double(v(:)');
    if numel(v) ~= n_elem
        error('load_transducer_from_library:correctionLength', ...
            ['%s has %d elements but the transducer geometry has %d.\n' ...
             'The entry was built for a different geometry — rebuild it ' ...
             '(build_transducer_library) or remove the correction:\n  %s'], ...
            what, numel(v), n_elem, source_label);
    end
end

function [depths, offs, amps, blend_mm, pamp, isppa, gain, amp_src, kw_r2, kw_fe] = ...
    collect_refinements(lib, n_elem, geo_at, source_label)
% Gather refinement offsets [deg] per depth from the current and legacy forms.
% pamp: per-depth drive [Pa] at isppa [W/cm²]; gain: [K x N_elem] per-ring
% gains (current form only; NaN = none). amp_src: 'kwave' | 'analytical' | ''
% per depth (KWAVE_CALIBRATE_DRIVE provenance); kw_r2/kw_fe: the k-Wave
% check's profile R² / focal error [mm] at that depth, NaN if unchecked.
    depths = []; offs = zeros(0, n_elem); amps = {}; blend_mm = 5; pamp = []; isppa = NaN;
    gain = zeros(0, n_elem); amp_src = {}; kw_r2 = []; kw_fe = [];

    if isfield(lib, 'refinement') && isstruct(lib.refinement) && ...
            isfield(lib.refinement, 'depths_ep_mm') && ~isempty(lib.refinement.depths_ep_mm)
        r = lib.refinement;
        d = double(r.depths_ep_mm(:)');
        c = rows_of(r.elem_phase_correction_deg, numel(d));
        for k = 1:numel(d)
            depths(end+1) = d(k); %#ok<AGROW>
            offs(end+1, :) = check_length(c{k}, n_elem, ...
                sprintf('refinement.elem_phase_correction_deg{%d}', k), source_label); %#ok<AGROW>
            amps{end+1} = []; %#ok<AGROW>
        end
        pamp = nan(1, numel(d));
        if isfield(r, 'elem_p_amp') && numel(r.elem_p_amp) == numel(d) && ...
                isfield(r, 'isppa_wcm2') && ~isempty(r.isppa_wcm2)
            pamp  = double(r.elem_p_amp(:)');
            isppa = double(r.isppa_wcm2);
        end
        amp_src = repmat({''}, 1, numel(d));
        if isfield(r, 'amp_source') && numel(r.amp_source) == numel(d)
            amp_src = cellstr(string(r.amp_source(:)'));
        end
        kw_r2 = nan(1, numel(d));
        if isfield(r, 'kwave_r2') && numel(r.kwave_r2) == numel(d)
            kw_r2 = double(r.kwave_r2(:)');
        end
        kw_fe = nan(1, numel(d));
        if isfield(r, 'kwave_focal_err_mm') && numel(r.kwave_focal_err_mm) == numel(d)
            kw_fe = double(r.kwave_focal_err_mm(:)');
        end
        gain = nan(numel(d), n_elem);
        if isfield(r, 'elem_amp_weight')
            gc = rows_of(r.elem_amp_weight, numel(d));
            for k = 1:numel(d)
                if all(isfinite(gc{k})); gain(k, :) = check_length(gc{k}, n_elem, ...
                        sprintf('refinement.elem_amp_weight{%d}', k), source_label); end
            end
        end
        if isfield(r, 'blend_mm') && ~isempty(r.blend_mm); blend_mm = double(r.blend_mm); end
    end

    % Legacy (b): absolute per-depth phase tables in global_model
    if isfield(lib, 'global_model') && isstruct(lib.global_model) && ...
            isfield(lib.global_model, 'depths_ep_mm') && ~isempty(lib.global_model.depths_ep_mm) && ...
            isfield(lib.global_model, 'elem_phase_deg') && ~isempty(lib.global_model.elem_phase_deg)
        d = double(lib.global_model.depths_ep_mm(:)');
        c = rows_of(lib.global_model.elem_phase_deg, numel(d));
        for k = 1:numel(d)
            abs_deg = check_length(c{k}, n_elem, ...
                sprintf('global_model.elem_phase_deg{%d}', k), source_label);
            [depths, offs, amps] = add_knot(depths, offs, amps, d(k), abs_deg - geo_at(d(k)), []);
        end
    end

    % Legacy (c): parametric per-depth entries in calibration.focal_depths.
    % Without a global model every entry is used; otherwise only is_refinement ones.
    if isfield(lib, 'calibration') && isstruct(lib.calibration) && ...
            isfield(lib.calibration, 'focal_depths') && isstruct(lib.calibration.focal_depths)
        has_gm = isfield(lib, 'global_model') && isstruct(lib.global_model) && ...
            ~isempty(fieldnames(lib.global_model));
        fd   = lib.calibration.focal_depths;
        keys = fieldnames(fd);
        for k = 1:numel(keys)
            e = fd.(keys{k});
            if has_gm && ~(isfield(e, 'is_refinement') && e.is_refinement); continue; end
            d_k = str2double(strrep(strrep(keys{k}, 'f', ''), 'p', '.'));
            abs_deg = parametric_phases(e, n_elem, source_label) * 180/pi;
            amp_k = [];
            if isfield(e, 'amplitude_scaling'); amp_k = e.amplitude_scaling; end
            [depths, offs, amps] = add_knot(depths, offs, amps, d_k, abs_deg - geo_at(d_k), amp_k);
        end
    end

    pamp(end+1:numel(depths)) = NaN;   % legacy knots carry no per-depth drive or gains
    gain(end+1:numel(depths), :) = NaN;
    amp_src(end+1:numel(depths)) = {''};   % legacy knots were never k-Wave-calibrated
    kw_r2(end+1:numel(depths)) = NaN;
    kw_fe(end+1:numel(depths)) = NaN;
    [depths, order] = sort(depths);
    offs = offs(order, :);
    amps = amps(order);
    pamp = pamp(order);
    gain = gain(order, :);
    amp_src = amp_src(order);
    kw_r2 = kw_r2(order);
    kw_fe = kw_fe(order);
end

function [depths, offs, amps] = add_knot(depths, offs, amps, d, off, amp)
    % Later sources do not override an explicit refinement at the same depth
    if any(abs(depths - d) < 0.5); return; end
    depths(end+1)  = d;
    offs(end+1, :) = off;
    amps{end+1}    = amp;
end

function c = rows_of(v, n_rows)
    if iscell(v)
        c = v;
    elseif n_rows == 1
        c = {v};
    else
        c = num2cell(v, 2);
    end
end

function ph = parametric_phases(e, n_elem, source_label)
    start = e.phase_start_deg * pi/180;
    switch char(e.precession_mode)
        case 'linear'
            ph = mod(start + (0:n_elem-1) * e.phase_step_deg * pi/180, 2*pi);
        case 'monotonic'
            ph = mod(start + [0, cumsum(e.phase_increments_deg(:)' * pi/180)], 2*pi);
        otherwise
            error('load_transducer_from_library: unknown precession_mode ''%s'' in\n  %s', ...
                char(e.precession_mode), source_label);
    end
end

function [off, w] = refinement_at(depths, offs, d, blend_mm)
% Refinement offset at depth d and its weight against the static offset.
    off = zeros(1, size(offs, 2)); w = 0;
    if isempty(depths); return; end
    if d >= depths(1) && d <= depths(end)
        if numel(depths) == 1
            off = offs(1, :);
        else
            % per-element interpolation on the unit circle via unwrapping over depth
            u   = unwrap(offs * pi/180, [], 1);
            off = interp1(depths(:), u, d, 'linear') * 180/pi;
        end
        w = 1;
    else
        [dist, k] = min(abs(depths - d));
        off = offs(k, :);
        w   = max(0, 1 - dist / max(blend_mm, eps));
    end
end

function drive = lookup_amplitude(amp_map, desired_intensity)
    int_keys   = fieldnames(amp_map);
    int_values = cellfun(@(k) str2double(strrep(strrep(k, 'i', ''), 'p', '.')), int_keys);
    [int_dist, idx] = min(abs(int_values - desired_intensity));
    base_amp = double(amp_map.(int_keys{idx}));
    if isnan(desired_intensity) || int_dist / desired_intensity <= 0.05
        drive = base_amp;
    else
        % I ∝ amp² in the linear regime — exact scaling, no warning needed
        drive = round(base_amp * sqrt(desired_intensity / int_values(idx)));
    end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end
