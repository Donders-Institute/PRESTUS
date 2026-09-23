function [tic_val, info] = compute_tic(parameters, transducer, Isppa_ref_Wcm2)
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
% W0 ASSUMPTION (the one modelling choice here, deliberately isolated in the
% local function estimate_W0_mW below so it is easy to review/replace):
%   PRESTUS does not store emitted acoustic power. We estimate it as
%       W0 = Isppa_ref * A_aperture * duty_cycle
%   i.e. a reference spatial-peak pulse-average intensity times the geometric
%   aperture area times the duty cycle. Because Isppa_ref is a spatial PEAK
%   (focal) intensity rather than the aperture-averaged source intensity, this
%   OVERESTIMATES the true emitted power by ~the focusing gain and should be
%   read as an upper bound. To make TIC physically exact, replace
%   estimate_W0_mW with either (a) the source intensity (elem_amp^2/2 rho c)
%   times the active element area, or (b) the integral of the temporal-average
%   intensity over a transverse plane of a free-field/water reference run.
%
% Use as:
%   [tic_val, info] = compute_tic(parameters, parameters.transducer(1), results.Isppa)
%
% Input:
%   parameters     - PRESTUS parameters (used for the duty cycle)
%   transducer     - a single transducer struct, e.g. parameters.transducer(1)
%   Isppa_ref_Wcm2 - reference spatial-peak pulse-average intensity [W/cm^2]
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

    % --- emitted time-averaged acoustic power W0 (see ASSUMPTION above) ---
    [W0_mW, note] = estimate_W0_mW(Isppa_ref_Wcm2, Deq_cm, dc);
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
    if ~isfield(tr, 'type') || ~ischar(tr.type) && ~isstring(tr.type), return; end
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
function dc = local_duty_cycle(p)
% Duty cycle in [0,1]; defaults to 1 (continuous) when no pulsing is defined.
    dc = 1;
    srcs = {};
    if isfield(p, 'timing')  && isstruct(p.timing),  srcs{end+1} = p.timing;  end
    if isfield(p, 'thermal') && isstruct(p.thermal), srcs{end+1} = p.thermal; end
    for i = 1:numel(srcs)
        s = srcs{i};
        if isfield(s, 'dc') && isnumeric(s.dc) && isscalar(s.dc) && s.dc > 0 && s.dc <= 1
            dc = s.dc; return
        end
        if isfield(s, 'pd') && isfield(s, 'pri') && isnumeric(s.pd) && isnumeric(s.pri) ...
                && isscalar(s.pd) && isscalar(s.pri) && s.pri > 0
            cand = s.pd / s.pri;
            if cand > 0 && cand <= 1, dc = cand; return; end
        end
    end
end

% ------------------------------------------------------------------------
function [W0_mW, note] = estimate_W0_mW(Isppa_Wcm2, Deq_cm, dc)
% Estimate time-averaged emitted acoustic power [mW]. See ASSUMPTION in the
% header — this is the single place to refine the power model.
    W0_mW = NaN;
    note  = '';
    if isnan(Isppa_Wcm2) || isnan(Deq_cm) || isnan(dc)
        note = 'insufficient inputs for W0 estimate';
        return
    end
    area_cm2   = pi * (Deq_cm / 2)^2;             % geometric aperture area [cm^2]
    W0_pulse_W = Isppa_Wcm2 * area_cm2;           % pulse-average power [W] (upper bound)
    W0_mW      = W0_pulse_W * dc * 1000;          % time-averaged power [mW]
    note = 'W0 ~ Isppa_peak * aperture_area * duty_cycle (upper bound; informational)';
end
