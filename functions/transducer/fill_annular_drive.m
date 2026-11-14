function tr = fill_annular_drive(tr, c_water, t_i)
% FILL_ANNULAR_DRIVE  Steering phases and drive of an annular transducer from its depth setting
%
% A transducer defined by its geometry (no serial, or a library entry that left
% them empty) can be driven by its focal depth setting and the target free-water
% ISPPA alone, as a library device is:
%   - no element phases: geometric steering phases to transducer.focal_distance_ep
%     (or focal_distance_bowl), from the ring geometry (compute_phases, the
%     library's geometric baseline);
%   - no element amplitude and a target ISPPA: a nominal uniform drive
%     (NOMINAL_P_AMP_PA); the results are scaled to transducer.target_isppa_wcm2
%     (ISPPA targeting). Without a target ISPPA the amplitude stays required.
% Phases and amplitudes that are set are kept: they come from calibration.
% Flat rings (curv_radius_mm = Inf) are left as they are.
%
% Use as:
%   tr = fill_annular_drive(tr, c_water, t_i)
%
% Input:
%   tr      - one transducer struct (type 'annular')
%   c_water - sound speed in water for the steering phases [m/s]
%   t_i     - transducer index (messages)
%
% Output:
%   tr - with annular.elem_phase_deg / elem_p_amp filled in where they were empty
%
% See also: COMPUTE_PHASES, VALIDATE_ANNULAR_TRANSDUCER, RESOLVE_TRANSDUCER_FROM_SERIAL

    NOMINAL_P_AMP_PA = 100e3;   % any value: results are rescaled to target_isppa_wcm2

    if ~isfield(tr, 'annular') || ~isstruct(tr.annular)
        return;
    end
    a = tr.annular;
    has = @(s, f) isfield(s, f) && ~isempty(s.(f));

    if ~has(a, 'elem_phase_deg') && ~has(a, 'elem_phase_rad') && has(a, 'curv_radius_mm') ...
            && isfinite(a.curv_radius_mm) && has(a, 'elem_id_mm') && has(a, 'elem_od_mm') ...
            && has(tr, 'freq_hz')
        if ~has(a, 'dist_geom_ep_mm')                  % as validate_annular_transducer
            a.dist_geom_ep_mm = sqrt(a.curv_radius_mm^2 - (max(a.elem_od_mm) / 2)^2);
        end
        ep = [];
        if has(tr, 'focal_distance_ep')
            ep = tr.focal_distance_ep;
        elseif has(tr, 'focal_distance_bowl')
            ep = tr.focal_distance_bowl - (a.curv_radius_mm - a.dist_geom_ep_mm);
        end
        if ~isempty(ep)
            tran = struct('transducer', struct('freq_hz', tr.freq_hz, 'annular', a));
            a.elem_phase_deg = compute_phases(c_water, tran, ep);
            a.elem_phase_rad = deg2rad(a.elem_phase_deg);
            fprintf('[transducer %d] Geometric steering phases for a focal depth setting of %.1f mm.\n', ...
                t_i, ep);
        end
    end

    % Only with a target ISPPA: without one the amplitude sets the result's scale.
    if ~has(a, 'elem_p_amp') && ~has(a, 'elem_amp') && has(a, 'elem_n') && has(tr, 'target_isppa_wcm2')
        a.elem_p_amp = repmat(NOMINAL_P_AMP_PA, 1, a.elem_n);
        fprintf(['[transducer %d] Nominal drive (%g kPa per ring); results are scaled to ' ...
                 'target_isppa_wcm2.\n'], t_i, NOMINAL_P_AMP_PA / 1e3);
    end
    tr.annular = a;
end
