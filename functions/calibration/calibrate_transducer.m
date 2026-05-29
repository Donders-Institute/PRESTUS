function result = calibrate_transducer(parameters, charac_data, ...
    available_depths_ep, dist_bowl, desired_intensity, dist_bowl_offset, options)
% CALIBRATE_TRANSDUCER  Unified analytical entry point for transducer calibration
%
% Routes to the appropriate phase optimisation backend and forward model
% based on caller-specified options. Does NOT run k-Wave simulations —
% returns analytically calibrated phases and velocity that can be passed
% to CALIBRATION_TRANSDUCER for full simulation-based validation.
%
% Calibration mode:
%
%   'single_ref'  (default)
%     Runs PERFORM_GLOBAL_SEARCH at one reference depth. A hardware
%     correction offset (delta_phase_hw) is derived by subtracting the
%     geometric phases from the optimised phases at that depth. All other
%     depths are covered analytically via SET_REAL_PHASES + delta_phase_hw
%     without further optimisation.
%
% Two forward models:
%
%   'oneil'    (default)
%     O'Neil closed-form model via FOCUSEDANNULUSONEIL. Fast (~µs per
%     evaluation). Assumes uniform piston on a continuous spherical bowl.
%
%   'rayleigh'
%     Rayleigh–Sommerfeld numerical integral via RAYLEIGH_AXIAL_INTENSITY.
%     Slower but more accurate for large steering offsets and ring gaps.
%     Use for cross-validation against BabelBrain.
%
% Use as:
%   result = calibrate_transducer(parameters, charac_data, ...
%       available_depths_ep, dist_bowl, desired_intensity, dist_bowl_offset)
%   result = calibrate_transducer(..., options)
%
% Input:
%   parameters          - PRESTUS config with transducer.annular geometry,
%                         medium_properties.water, calibration.opt_* fields
%   charac_data         - [N_pts x N_depths] intensity matrix [W/cm²],
%                         columns indexed by available_depths_ep
%   available_depths_ep - [1 x N_depths] focal depths from exit plane [mm]
%   dist_bowl           - [N_pts x 1] axial positions from bowl [mm]
%   desired_intensity   - target peak Isppa [W/cm²]
%   dist_bowl_offset    - curv_radius_mm - dist_geom_ep_mm [mm]
%   options             - struct or name-value pairs with fields:
%     .Mode             - 'single_ref' (default; only supported mode)
%     .ForwardModel     - 'oneil' (default) | 'rayleigh'
%     .RefDepth         - reference depth [mm] for Mode='single_ref'
%                         (default: median of available_depths_ep)
%     .PhaseTable       - phase table struct for SET_REAL_PHASES
%                         (required for Mode='single_ref')
%     .Tran             - transducer struct from load_equipment_config
%                         (required for Mode='single_ref')
%     .ExportBabelBrain - path string to write OptimizedWeightsFile (.h5),
%                         or '' to skip (default: '')
%
% Output:
%   result  - struct with fields:
%     .opt_phases_rad       - [1 x N_elem] optimised phases [rad]
%     .opt_phases_deg       - [1 x N_elem] optimised phases [deg]
%     .opt_velocity         - calibrated scalar particle velocity [m/s]
%     .min_err              - final objective value
%     .mode                 - mode used
%     .forward_model        - forward model used
%     .delta_phase_hw_deg   - [1 x N_elem] hardware correction [deg]
%     .ref_depth_ep_mm      - reference depth used
%
% Example:
%   result = calibrate_transducer(parameters, charac_data, ...
%       available_depths_ep, dist_bowl, 32.5, dist_bowl_offset, ...
%       struct('RefDepth', 80, 'PhaseTable', phase_table, 'Tran', tran));
%
% See also: PERFORM_GLOBAL_SEARCH, RAYLEIGH_AXIAL_INTENSITY, SET_REAL_PHASES,
%           FIT_VELOCITY_TO_INTENSITY, EXPORT_BABELBRAIN_WEIGHTS,
%           CALIBRATION_TRANSDUCER

arguments
    parameters          (1,1) struct
    charac_data         (:,:) {mustBeNumeric}
    available_depths_ep (1,:) {mustBeNumeric}
    dist_bowl           (:,1) {mustBeNumeric}
    desired_intensity   (1,1) {mustBeNumeric}
    dist_bowl_offset    (1,1) {mustBeNumeric}
    options             (1,1) struct = struct()
end

%% Parse options ---------------------------------------------------------
mode          = get_opt(options, 'Mode',             'single_ref');
forward_model = get_opt(options, 'ForwardModel',     'oneil');
ref_depth     = get_opt(options, 'RefDepth',         median(available_depths_ep));
phase_table   = get_opt(options, 'PhaseTable',       []);
tran          = get_opt(options, 'Tran',             []);
export_path   = get_opt(options, 'ExportBabelBrain', '');

assert(strcmp(mode, 'single_ref'), ...
    'calibrate_transducer: only Mode=''single_ref'' is supported');
assert(ismember(forward_model, {'oneil', 'rayleigh'}), ...
    'calibrate_transducer: ForwardModel must be ''oneil'' or ''rayleigh''');
assert(~isempty(phase_table), ...
    'calibrate_transducer: PhaseTable is required');
assert(~isempty(tran), ...
    'calibrate_transducer: Tran is required');

fprintf('\n=== calibrate_transducer | ForwardModel: %s ===\n', forward_model);

%% Inject forward-model override into parameters -------------------------
parameters.calibration.forward_model = forward_model;

%% Calibrate -------------------------------------------------------------
ref_col = find(abs(available_depths_ep - ref_depth) < 0.5, 1);
        if isempty(ref_col)
            error(['calibrate_transducer: RefDepth %.1f mm not found. ' ...
                'Available: %s'], ref_depth, mat2str(available_depths_ep));
        end

        % Empirical profile at reference depth, scaled to desired intensity
        emp_ref = charac_data(:, ref_col);
        emp_ref(emp_ref < 0) = 0;
        emp_ref = emp_ref * (desired_intensity / max(emp_ref));

        profile_ref.axial_intensity     = emp_ref;
        profile_ref.axial_distance_bowl = dist_bowl;

        % Parameters for reference depth
        params_ref = parameters;
        params_ref.transducer.focal_distance_ep          = ref_depth;
        params_ref.transducer.focal_distance_bowl        = ref_depth + dist_bowl_offset;
        params_ref.calibration.desired_focal_distance_ep = ref_depth;
        params_ref.calibration.desired_intensity         = desired_intensity;
        params_ref = load_parameters(params_ref);

        % Geometric phases at reference depth
        geo_phases_ref_deg = set_real_phases(phase_table, tran, ref_depth, params_ref);
        params_ref.transducer.annular.elem_phase_deg = geo_phases_ref_deg;
        params_ref.transducer.annular.elem_phase_rad = geo_phases_ref_deg * pi / 180;

        % Global search at reference depth — seed from geometric phases so
        % the recovered delta is a small hardware-specific correction, not
        % an arbitrary alternative solution far from the geometric prediction.
        params_ref.calibration.opt_use_initial_phases = true;
        fprintf('Running global search at reference depth %.1f mm...\n', ref_depth);
        [opt_phases_rad, opt_velocity, min_err] = perform_global_search(...
            params_ref, profile_ref, ...
            get_initial_velocity(parameters));

        % Velocity correction to match desired intensity
        [opt_velocity, ~, ~] = fit_velocity_to_intensity(...
            params_ref, profile_ref, opt_phases_rad, opt_velocity, desired_intensity);

        % Hardware correction offset
        geo_unwrapped    = unwrap(geo_phases_ref_deg * pi/180) * 180/pi;
        opt_unwrapped    = unwrap(opt_phases_rad / pi * 180);
        delta_hw_deg     = mod(opt_unwrapped - geo_unwrapped, 360);

        fprintf('Hardware correction [deg]: %s\n', mat2str(round(delta_hw_deg)));

result.delta_phase_hw_deg = delta_hw_deg;
result.ref_depth_ep_mm    = ref_depth;

%% Package result --------------------------------------------------------
result.opt_phases_rad = opt_phases_rad;
result.opt_phases_deg = opt_phases_rad / pi * 180;
result.opt_velocity   = opt_velocity;
result.min_err        = min_err;
result.mode           = mode;
result.forward_model  = forward_model;

fprintf('\nCalibration complete.\n');
fprintf('  Phases [deg]: %s\n', mat2str(round(result.opt_phases_deg)));
fprintf('  Velocity:     %.4f m/s\n', opt_velocity);
fprintf('  Error:        %.6f\n', min_err);

%% Optional BabelBrain export -------------------------------------------
if ~isempty(export_path)
    export_babelbrain_weights(opt_phases_rad, opt_velocity, parameters, export_path);
end

end

%% Local helpers ---------------------------------------------------------
function v = get_opt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end

function v0 = get_initial_velocity(parameters)
    if isfield(parameters.calibration, 'initial_velocity') && ...
            ~isempty(parameters.calibration.initial_velocity)
        v0 = parameters.calibration.initial_velocity;
    else
        v0 = 0.05;
    end
end
