function profile_sim = extract_analytical_profile(parameters)
% EXTRACT_ANALYTICAL_PROFILE  Compute an axial intensity profile analytically (no k-Wave simulation)
%
% Mirrors extract_simulated_profile but derives all quantities from a
% forward model instead of a k-Wave result.  Velocity and the axial
% distance grid follow the same formulas used in extract_simulated_profile
% so the two profiles are directly comparable.
%
% The forward model is selected by parameters.calibration.forward_model:
%   'oneil'    (default) — O'Neil closed-form via focusedAnnulusONeil
%   'rayleigh'           — Rayleigh–Sommerfeld via rayleigh_axial_intensity
% This mirrors the same flag used in phase_optimization_annulus_full_curve.
%
% When this profile is passed to compute_analytical_solution the resulting
% simulated_analytical_scaling equals 1 by construction, so amplitude
% calibration reduces to the pure I ∝ v² analytical relationship.
%
% Use as:
%   profile_sim = extract_analytical_profile(parameters)
%
% Input:
%   parameters - PRESTUS config with transducer.annular geometry
%                (elem_amp, elem_phase_rad, elem_id_mm, elem_od_mm,
%                curv_radius_mm, elem_n, freq_hz),
%                medium_properties.water (density, sound_speed),
%                calibration.desired_focal_distance_ep [mm],
%                calibration.forward_model ('oneil'|'rayleigh', optional),
%                and optionally grid.resolution_mm [mm]
%
% Output:
%   profile_sim - struct with axial_intensity [W/cm²], axial_distance_bowl [mm],
%                 velocity [m/s]  — same fields as extract_simulated_profile
%
% See also: EXTRACT_SIMULATED_PROFILE, CALIBRATION_TRANSDUCER, COMPUTE_ANALYTICAL_SOLUTION,
%           RAYLEIGH_AXIAL_INTENSITY, PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE

    rho = parameters.medium_properties.water.density;
    c   = parameters.medium_properties.water.sound_speed;

    if isfield(parameters.grid, 'resolution_mm') && ~isempty(parameters.grid.resolution_mm)
        dx = parameters.grid.resolution_mm;
    else
        dx = 0.5;
    end

    if isfield(parameters.calibration, 'forward_model')
        forward_model = parameters.calibration.forward_model;
    else
        forward_model = 'oneil';
    end

    profile_sim.velocity            = parameters.transducer.annular.elem_amp(1) / (rho * c);
    profile_sim.axial_distance_bowl = 0 : dx : parameters.calibration.desired_focal_distance_ep * 2;

    if strcmp(forward_model, 'rayleigh')
        profile_sim.axial_intensity = rayleigh_axial_intensity( ...
            parameters.transducer.annular.elem_phase_rad, profile_sim.velocity, ...
            parameters, profile_sim.axial_distance_bowl);
    elseif strcmp(forward_model, 'oneil')
        p_axial = focusedAnnulusONeil( ...
            parameters.transducer.annular.curv_radius_mm / 1e3, ...
            [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
            repmat(profile_sim.velocity, 1, parameters.transducer.annular.elem_n), ...
            parameters.transducer.annular.elem_phase_rad, ...
            parameters.transducer.freq_hz, ...
            c, rho, ...
            (profile_sim.axial_distance_bowl - 0.5) * 1e-3);
        profile_sim.axial_intensity = p_axial .^ 2 / (2 * c * rho) * 1e-4;
    else
        error('extract_analytical_profile: unknown forward_model ''%s''; expected ''oneil'' or ''rayleigh''.', forward_model);
    end

end
