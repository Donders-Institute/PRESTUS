function [phases_rad, velocity_out, fit_info] = calibrate_babelbrain_phase_profile(tr, parameters, opts)
% CALIBRATE_BABELBRAIN_PHASE_PROFILE  Whole-profile Rayleigh fit of annular drive phases for BabelBrain
%
% Optimises the per-element drive phases so that the discrete-ring
% Rayleigh-Sommerfeld axial field (RAYLEIGH_AXIAL_INTENSITY) reproduces the
% continuous-bowl O'Neil profile (focusedAnnulusONeil) across the whole focal
% region, rather than matching a single landmark as CALIBRATE_FOCAL_STEERING
% does. Target and forward model differ only by the ring-gap discretisation,
% so the fit is non-circular and needs no external (k-Wave or hydrophone)
% data. Used when babelbrain.profile_calibration is enabled.
%
% The result is BabelBrain-specific: the returned phases feed only the weights
% delta (EXPORT_BABELBRAIN_WEIGHTS). The shared transducer-library phases used
% by k-Wave and jwave are untouched.
%
% Use as:
%   [phases_rad, velocity_out, fit_info] = calibrate_babelbrain_phase_profile(tr, parameters)
%   [phases_rad, velocity_out, fit_info] = calibrate_babelbrain_phase_profile(tr, parameters, opts)
%
% Input:
%   tr          - parameters.transducer(idx) struct with .annular geometry
%                 (elem_n, elem_id_mm, elem_od_mm, curv_radius_mm,
%                 dist_geom_ep_mm, elem_phase_rad, elem_p_amp) and
%                 .focal_distance_ep, .freq_hz.
%   parameters  - PRESTUS config (uses medium_properties.water; optional
%                 calibration.opt_* knobs below are read via a local copy).
%   opts        - optional struct:
%     .eval_max_mm   - axial evaluation extent from bowl [mm]
%                      (default max(130, 1.6*(focal_ep + bowl_off))).
%     .step_mm       - axial sample step [mm] (default 1.0).
%     .opt_limits    - [lo hi] bowl-distance window for the fit [mm]
%                      (default: full evaluated grid = whole profile).
%     .lambda        - L2 regularisation toward the library phases
%                      (default 0.5; keeps the fit a correction, not a
%                      free re-solve).
%     .max_fun_evals - optimiser budget (default 3e4; the Rayleigh forward
%                      model is ~ms/eval, so this is the main cost knob).
%     .seed          - RNG seed for reproducibility (default 1).
%
% Output:
%   phases_rad   - [1 x N_elem] profile-optimised drive phases [rad].
%   velocity_out - re-fitted scalar particle velocity [m/s] (informational;
%                  the absolute drive amplitude for BabelBrain still comes
%                  from source_amp_pa / elem_p_amp).
%   fit_info     - struct: target_peak_ep_mm, fit_peak_ep_mm, fit_flhm_ep_mm
%                  (-3 dB FLHM centre), library_peak_ep_mm, min_err, n_pts,
%                  opt_limits_mm.
%
% See also: PERFORM_GLOBAL_SEARCH, PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE,
%           RAYLEIGH_AXIAL_INTENSITY, EXPORT_BABELBRAIN_WEIGHTS,
%           CALIBRATE_FOCAL_STEERING, COMPUTE_FLHM_CENTER

arguments
    tr          (1,1) struct
    parameters  (1,1) struct
    opts        (1,1) struct = struct()
end

    an       = tr.annular;
    n_elem   = an.elem_n;
    c        = parameters.medium_properties.water.sound_speed;
    rho      = parameters.medium_properties.water.density;

    % Exit-plane <-> bowl offset: bowl distance = ep distance + (curv_radius - dist_geom_ep)
    bowl_off_mm = an.curv_radius_mm - an.dist_geom_ep_mm;
    focal_ep_mm = tr.focal_distance_ep;

    if ~isfield(opts, 'eval_max_mm') || isempty(opts.eval_max_mm)
        opts.eval_max_mm = max(130, 1.6 * (focal_ep_mm + bowl_off_mm));
    end
    if ~isfield(opts, 'step_mm')       || isempty(opts.step_mm);       opts.step_mm = 1.0;   end
    if ~isfield(opts, 'lambda')        || isempty(opts.lambda);        opts.lambda = 0.5;    end
    if ~isfield(opts, 'max_fun_evals') || isempty(opts.max_fun_evals); opts.max_fun_evals = 3e4; end
    if ~isfield(opts, 'seed')          || isempty(opts.seed);          opts.seed = 1;        end

    % Axial grid measured from the bowl surface (rayleigh_axial_intensity and
    % focusedAnnulusONeil both take bowl-referenced positions).
    dist_bowl_mm = (opts.step_mm : opts.step_mm : opts.eval_max_mm)';
    n_pts        = numel(dist_bowl_mm);

    % Library drive phases (the shared, per-solver-agnostic calibration result).
    lib_phases = an.elem_phase_rad(:)';

    % Initial particle velocity from the calibrated reference-ring pressure.
    if isfield(an, 'elem_p_amp') && ~isempty(an.elem_p_amp)
        velocity0 = an.elem_p_amp(1) / (rho * c);
    else
        velocity0 = 0.05;   % m/s, benign starting scale (shape fit, not amplitude)
    end

    % --- Target: O'Neil continuous-bowl ideal at the library drive ------------
    % Mirror the O'Neil branch of phase_optimization_annulus_full_curve exactly
    % (same (dist-0.5) sampling convention) so target and forward model live on
    % the same axis grid.
    p_axial_target = focusedAnnulusONeil( ...
        an.curv_radius_mm / 1e3, ...
        [an.elem_id_mm; an.elem_od_mm] / 1e3, ...
        repmat(velocity0, 1, n_elem), ...
        lib_phases, ...
        tr.freq_hz, c, rho, ...
        (dist_bowl_mm - 0.5) * 1e-3);
    target_I = p_axial_target(:).^2 / (2 * c * rho) * 1e-4;   % [W/cm^2]

    profile_target = struct('axial_distance_bowl', dist_bowl_mm, ...
                            'axial_intensity',     target_I);

    % --- Fit: Rayleigh forward model, whole-curve, phases re-optimised ---------
    % Build a LOCAL parameters copy configured for the reusable optimiser stack.
    p = parameters;
    p.transducer          = tr;                      % perform_global_search reads .transducer.annular
    p.transducer.annular  = an;
    p.calibration         = struct();
    p.calibration.forward_model            = 'rayleigh';     % discrete-ring forward model
    p.calibration.opt_weights              = 0;              % uniform => whole profile
    if isfield(opts, 'opt_limits') && ~isempty(opts.opt_limits)
        p.calibration.opt_limits           = opts.opt_limits;
    else
        p.calibration.opt_limits           = [dist_bowl_mm(1), dist_bowl_mm(end)];
    end
    p.calibration.opt_use_initial_phases   = true;           % start from library phases
    p.calibration.opt_regularization_lambda= opts.lambda;    % stay a correction, not a re-solve
    p.calibration.opt_seed                 = opts.seed;
    p.calibration.opt_max_fun_evals        = opts.max_fun_evals;
    p.calibration.opt_upper_velocity       = max(0.2, 4 * velocity0);

    % perform_global_search unconditionally saves a debug figure to
    % parameters.io.outputs_folder/img_calibration. The handoff `parameters` does
    % not carry that field, so route it to a scratch dir (the figure is a
    % throwaway here; the fit result is what we use).
    if ~isfield(p, 'io') || ~isstruct(p.io)
        p.io = struct();
    end
    if ~isfield(p.io, 'outputs_folder') || isempty(p.io.outputs_folder)
        if isfield(opts, 'outputs_folder') && ~isempty(opts.outputs_folder)
            p.io.outputs_folder = opts.outputs_folder;
        else
            p.io.outputs_folder = tempname;   % created (recursively) by perform_global_search
        end
    end

    [opt_phases, velocity_out, min_err] = perform_global_search(p, profile_target, velocity0);

    phases_rad = mod(opt_phases(:)', 2*pi);

    % --- Diagnostics -----------------------------------------------------------
    ep = @(bowl_mm) bowl_mm - bowl_off_mm;   % bowl distance -> exit-plane distance
    I_fit = rayleigh_axial_intensity(phases_rad, velocity_out, p, dist_bowl_mm);
    I_lib = rayleigh_axial_intensity(lib_phases, velocity0, p, dist_bowl_mm);
    [~, i_target] = max(target_I);
    [~, i_fit]    = max(I_fit);
    [~, i_lib]    = max(I_lib);
    fit_flhm_mm   = compute_flhm_center(I_fit(:), dist_bowl_mm);

    fit_info = struct( ...
        'target_peak_ep_mm',  ep(dist_bowl_mm(i_target)), ...
        'fit_peak_ep_mm',     ep(dist_bowl_mm(i_fit)), ...
        'fit_flhm_ep_mm',     ep(fit_flhm_mm), ...
        'library_peak_ep_mm', ep(dist_bowl_mm(i_lib)), ...
        'min_err',            min_err, ...
        'n_pts',              n_pts, ...
        'opt_limits_mm',      p.calibration.opt_limits);

end
