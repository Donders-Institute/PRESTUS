function [tic_val, info] = compute_tic(parameters, transducer, Isppa_ref_Wcm2, scale_p)
% COMPUTE_TIC  Cranial Thermal Index (TIC) after IEC 62359 / AIUM-NEMA ODS.
%
%   TIC = W0 / (C_TIC * D_eq)
%       C_TIC = 40 mW/cm   (cranial-bone constant, IEC 62359)
%       W0    = time-averaged emitted acoustic power [mW]
%       D_eq  = equivalent aperture diameter [cm]
%
% TIC is reported as an INFORMATIONAL index: ITRUSST (Aubry et al., 2025)
% bases non-significant-risk on temperature rise, absolute temperature and
% CEM43 rather than the output-display thermal indices, so there is no hard
% TIC limit. The value is provided for cross-referencing with device output.
%
% W0 MODEL (isolated in the local function estimate_W0_mW below):
%   PRESTUS does not store emitted acoustic power, so W0 is estimated from
%   the source definition:
%       W0 = sum_i (scale_p * elem_p_amp_i)^2 / (2 rho c) * A_i * ramp * duty_cycle
%   i.e. the plane-wave source intensity of each element (water rho and c
%   from parameters.medium_properties.water) times its projected area.
%   scale_p is the pressure factor of ISPPA targeting (APPLY_ISPPA_SCALING),
%   so W0 belongs to the reported, scaled field. Annular arrays use the ring
%   areas from elem_od_mm / elem_id_mm; other types use the full aperture
%   area with the RMS amplitude (an upper bound for sparse arrays).
%   ramp is the pulse-average / peak intensity factor of a ramped pulse
%   (timing.pulse_ramp_dur, COMPUTE_RAMP_FACTOR); 1 for rectangular pulses.
%   FALLBACK without elem_p_amp: W0 = Isppa_ref * A_aperture * duty_cycle
%   (Isppa_ref is a pulse average, so the ramp is already in it).
%   Isppa_ref is a focal (spatial-peak) intensity, so the fallback
%   overestimates W0 by about the focusing gain and is only an upper bound.
%   info.note records which model ran.
%
% DUTY CYCLE: pulse duty (timing.dc, or pd/pri) times pulse-train duty
% (ptd/ptri) when both are defined; 1 (continuous) when nothing is defined.
%
% Use as:
%   [tic_val, info] = compute_tic(parameters, parameters.transducer(1), results.Isppa)
%   [tic_val, info] = compute_tic(parameters, parameters.transducer(1), results.Isppa, scale_p)
%
% Input:
%   parameters     - PRESTUS parameters (used for the duty cycle)
%   transducer     - a single transducer struct, e.g. parameters.transducer(1)
%   Isppa_ref_Wcm2 - reference spatial-peak pulse-average intensity [W/cm^2]
%   scale_p        - pressure factor of ISPPA targeting (sensor_data.isppa_scale_p);
%                    default 1 (the field as simulated)
%
% Output:
%   tic_val - cranial thermal index [-], or NaN if inputs are unavailable
%   info    - struct with fields W0_mW, Deq_cm, duty_cycle, note (provenance)
%
% See also: GENERATE_SIMULATION_REPORT, GET_RISK_LIMITS, ACOUSTIC_ANALYSIS

    arguments
        parameters       (1,1) struct
        transducer       (1,1) struct
        Isppa_ref_Wcm2   (1,1) double
        scale_p          (1,1) double = 1
    end

    C_TIC   = 40;  % mW/cm  (IEC 62359 cranial constant)
    tic_val = NaN;
    info    = struct('W0_mW', NaN, 'Deq_cm', NaN, 'duty_cycle', NaN, 'note', '');

    % --- equivalent aperture diameter D_eq ---
    Deq_mm = local_aperture_mm(transducer);
    if isnan(Deq_mm) || Deq_mm <= 0
        info.note = 'TIC unavailable: transducer aperture diameter not found';
        return
    end
    Deq_cm = Deq_mm / 10;
    info.Deq_cm = Deq_cm;

    % --- duty cycle (time-average factor) ---
    dc = local_duty_cycle(parameters);
    info.duty_cycle = dc;

    % --- emitted time-averaged acoustic power W0 (see W0 MODEL above) ---
    [W0_mW, note] = estimate_W0_mW(parameters, transducer, Isppa_ref_Wcm2, Deq_cm, dc, scale_p);
    info.W0_mW = W0_mW;
    info.note  = note;
    if isnan(W0_mW)
        return
    end

    tic_val = W0_mW / (C_TIC * Deq_cm);
end

% ------------------------------------------------------------------------
function d_mm = local_aperture_mm(tr)
% Equivalent aperture diameter [mm] for annular or matrix transducers.
    d_mm = NaN;
    if ~isfield(tr, 'type') || (~ischar(tr.type) && ~isstring(tr.type)), return; end
    t = char(tr.type);
    if ~isfield(tr, t), return; end
    sub = tr.(t);
    if strcmp(t, 'annular') && isfield(sub, 'elem_od_mm') && ~isempty(sub.elem_od_mm)
        v = sub.elem_od_mm(:);
        d_mm = max(v);                       % outer diameter of the largest ring
    elseif isfield(sub, 'outer_diameter_mm')
        d_mm = sub.outer_diameter_mm;        % matrix / clover arrays
    elseif isfield(sub, 'aperture_diameter_mm')
        d_mm = sub.aperture_diameter_mm;
    end
end

% ------------------------------------------------------------------------
function f = local_ramp_factor(p)
% Pulse-average / peak intensity factor of a ramped pulse (COMPUTE_RAMP_FACTOR,
% as for Ipa in ACOUSTIC_ANALYSIS); 1 without a ramp or a usable timing block.
    f = 1;
    if ~isfield(p, 'timing') || ~isfield(p.timing, 'pulse_ramp_dur') || ~isfield(p.timing, 'pd')
        return
    end
    r = p.timing.pulse_ramp_dur;
    if ~isnumeric(r) || ~isscalar(r) || ~(r > 0)
        return
    end
    shape = "linear";
    if isfield(p.timing, 'pulse_ramp_shape') && ~isempty(p.timing.pulse_ramp_shape)
        shape = string(p.timing.pulse_ramp_shape);
    end
    try
        f = compute_ramp_factor(p.timing.pd, r, shape);
    catch
        f = 1;
    end
end

% ------------------------------------------------------------------------
function dc = local_duty_cycle(p)
% Time-average factor in [0,1]: pulse duty (dc, or pd/pri) times pulse-train
% duty (ptd/ptri) from the same struct. Defaults to 1 (continuous) when no
% pulsing is defined.
    dc = 1;
    srcs = {};
    if isfield(p, 'timing')  && isstruct(p.timing),  srcs{end+1} = p.timing;  end
    if isfield(p, 'thermal') && isstruct(p.thermal), srcs{end+1} = p.thermal; end
    for i = 1:numel(srcs)
        s = srcs{i};
        pulse = local_ratio(s, 'pd', 'pri');
        if isfield(s, 'dc') && isnumeric(s.dc) && isscalar(s.dc) && s.dc > 0 && s.dc <= 1
            pulse = s.dc;
        end
        if isnan(pulse), continue; end
        train = local_ratio(s, 'ptd', 'ptri');
        if isnan(train), train = 1; end
        dc = pulse * train;
        return
    end
end

% ------------------------------------------------------------------------
function r = local_ratio(s, num, den)
% s.(num)/s.(den) when both are positive scalars and the ratio is in (0,1].
    r = NaN;
    if isfield(s, num) && isfield(s, den) && isnumeric(s.(num)) && isnumeric(s.(den)) ...
            && isscalar(s.(num)) && isscalar(s.(den)) && s.(den) > 0
        cand = s.(num) / s.(den);
        if cand > 0 && cand <= 1, r = cand; end
    end
end

% ------------------------------------------------------------------------
function [W0_mW, note] = estimate_W0_mW(parameters, tr, Isppa_Wcm2, Deq_cm, dc, scale_p)
% Estimate time-averaged emitted acoustic power [mW]. See W0 MODEL in the
% header: this is the single place to refine the power model.
    W0_mW = NaN;
    if isnan(Deq_cm) || isnan(dc)
        note = 'insufficient inputs for W0 estimate';
        return
    end
    P_pulse_W = local_source_power_W(parameters, tr, Deq_cm, scale_p);
    if ~isnan(P_pulse_W)
        % elem_p_amp is the peak amplitude: a ramp lowers the pulse-average power
        W0_mW = P_pulse_W * local_ramp_factor(parameters) * dc * 1000;
        note  = 'W0 = source intensity (elem_p_amp^2/2rhoc, ISPPA-scaled) * element area * ramp * duty_cycle (informational)';
        return
    end
    if isnan(Isppa_Wcm2)
        note = 'insufficient inputs for W0 estimate';
        return
    end
    area_cm2   = pi * (Deq_cm / 2)^2;             % geometric aperture area [cm^2]
    W0_pulse_W = Isppa_Wcm2 * area_cm2;           % pulse-average power [W] (upper bound)
    W0_mW      = W0_pulse_W * dc * 1000;          % time-averaged power [mW]
    note = 'W0 ~ Isppa_peak * aperture_area * duty_cycle (fallback: no elem_p_amp; upper bound; informational)';
end

% ------------------------------------------------------------------------
function P_W = local_source_power_W(parameters, tr, Deq_cm, scale_p)
% Pulse-average emitted power [W] from the source pressure amplitude(s):
% plane-wave intensity p^2/(2 rho c) per element times its projected area.
% NaN when elem_p_amp is unavailable.
    P_W = NaN;
    t = char(tr.type);
    if ~isfield(tr, t) || ~isstruct(tr.(t))
        return
    end
    g = tr.(t);
    if ~isfield(g, 'elem_p_amp') || ~isnumeric(g.elem_p_amp) || isempty(g.elem_p_amp) ...
            || any(~isfinite(g.elem_p_amp(:)))
        return
    end
    amp = double(g.elem_p_amp(:)) * scale_p;

    rho = 994; c = 1500;                          % config_default water
    if isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'water')
        w = parameters.medium_properties.water;
        if isfield(w, 'density')     && isnumeric(w.density)     && isscalar(w.density),     rho = w.density;     end
        if isfield(w, 'sound_speed') && isnumeric(w.sound_speed) && isscalar(w.sound_speed), c   = w.sound_speed; end
    end

    area_m2 = [];
    if strcmp(t, 'annular') && isfield(g, 'elem_od_mm') && ~isempty(g.elem_od_mm)
        od = double(g.elem_od_mm(:));
        id = zeros(size(od));
        if isfield(g, 'elem_id_mm') && numel(g.elem_id_mm) == numel(od)
            id = double(g.elem_id_mm(:));
        end
        area_m2 = pi/4 * max(od.^2 - id.^2, 0) * 1e-6;    % ring areas [m^2]
    end
    if isempty(area_m2) || (numel(amp) ~= 1 && numel(amp) ~= numel(area_m2))
        % whole aperture, mean-square amplitude
        area_m2 = pi * (Deq_cm / 2)^2 * 1e-4;
        amp = sqrt(mean(amp.^2));
    end
    P_W = sum(amp.^2 / (2 * rho * c) .* area_m2);
end
