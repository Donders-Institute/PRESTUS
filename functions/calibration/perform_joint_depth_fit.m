function [opt_delta_phases, opt_velocity, min_err] = perform_joint_depth_fit(...
    parameters, charac_data, available_depths_ep, dist_bowl, dist_bowl_offset, ...
    desired_intensity, phase_table, tran, options)
% PERFORM_JOINT_DEPTH_FIT  Multi-depth joint optimisation of hardware phase correction
%
% Optimises a single per-element delta phase (hardware correction) and
% particle velocity simultaneously across all calibration depths, by
% minimising the sum of profile errors over depths. This mirrors
% BabelBrain's TxCalibration approach: the transfer matrix is built with
% geometric (Rayleigh-steered) phases at each depth and a single correction
% vector b is fitted across all depths.
%
% The drive signal at depth j is:
%   phases_j = geo_phases(j) + delta_phases
%
% where geo_phases(j) are the geometric phases for focal depth j computed
% by SET_REAL_PHASES, and delta_phases is the depth-independent hardware
% correction being optimised. This is the same parameterisation used by
% CALIBRATE_TRANSDUCER Mode='single_ref' for cross-depth extrapolation,
% but here all depths constrain the fit jointly instead of extrapolating
% from one reference depth.
%
% Use as:
%   [opt_delta_phases, opt_velocity, min_err] = perform_joint_depth_fit(...
%       parameters, charac_data, available_depths_ep, dist_bowl, dist_bowl_offset, ...
%       desired_intensity, phase_table, tran)
%   [opt_delta_phases, opt_velocity, min_err] = perform_joint_depth_fit(...
%       ..., options)
%
% Input:
%   parameters          - PRESTUS config with transducer.annular geometry,
%                         medium_properties.water, calibration.opt_* fields
%   charac_data         - [N_pts x N_depths] intensity matrix [W/cm²]
%   available_depths_ep - [1 x N_depths] focal depths from exit plane [mm]
%   dist_bowl           - [N_pts x 1] axial positions from bowl [mm]
%   dist_bowl_offset    - curv_radius_mm - dist_geom_ep_mm [mm]
%   desired_intensity   - target peak Isppa [W/cm²] used for per-depth scaling
%   phase_table         - phase table struct for SET_REAL_PHASES
%   tran                - transducer struct from load_equipment_config
%   options             - struct with optional fields:
%     .ForwardModel     - 'rayleigh' (default) | 'oneil'
%     .Lambda           - L2 regularisation weight on delta phases (default: 0)
%
% Output:
%   opt_delta_phases - [1 x N_elem] depth-independent hardware phase
%                      correction [rad], relative to geometric phases
%   opt_velocity     - calibrated scalar particle velocity [m/s]
%   min_err          - final normalised objective value
%
% See also: CALIBRATE_TRANSDUCER, PERFORM_GLOBAL_SEARCH, SET_REAL_PHASES,
%           RAYLEIGH_AXIAL_INTENSITY

arguments
    parameters          (1,1) struct
    charac_data         (:,:) {mustBeNumeric}
    available_depths_ep (1,:) {mustBeNumeric}
    dist_bowl           (:,1) {mustBeNumeric}
    dist_bowl_offset    (1,1) {mustBeNumeric}
    desired_intensity   (1,1) {mustBeNumeric}
    phase_table         (1,1) struct
    tran                (1,1) struct
    options             (1,1) struct = struct()
end

%% Parse options
forward_model = get_opt(options, 'ForwardModel', 'rayleigh');
lambda        = get_opt(options, 'Lambda',        0);

assert(ismember(forward_model, {'oneil', 'rayleigh'}), ...
    'perform_joint_depth_fit: ForwardModel must be ''oneil'' or ''rayleigh''');

n_elem   = parameters.transducer.annular.elem_n;
n_depths = numel(available_depths_ep);
water    = parameters.medium_properties.water;

fprintf('\n=== perform_joint_depth_fit | ForwardModel: %s | Lambda: %.4g | %d depths ===\n', ...
    forward_model, lambda, n_depths);

%% Pre-compute geometric phases and normalised empirical profiles per depth
geo_phases_all = zeros(n_depths, n_elem);  % [n_depths x n_elem], rad
emp_all        = zeros(size(charac_data)); % [N_pts x n_depths], W/cm²
opt_limits_all = zeros(n_depths, 2);       % non-NaN range per depth [mm from bowl]

for j = 1:n_depths
    depth_j = available_depths_ep(j);

    % Geometric phases at this depth
    params_j = parameters;
    params_j.transducer.focal_distance_ep          = depth_j;
    params_j.transducer.focal_distance_bowl        = depth_j + dist_bowl_offset;
    params_j.calibration.desired_focal_distance_ep = depth_j;
    params_j.calibration.desired_intensity         = desired_intensity;
    params_j = load_parameters(params_j);

    geo_deg_j = set_real_phases(phase_table, tran, depth_j, params_j);
    geo_phases_all(j, :) = geo_deg_j * pi / 180;

    % Normalise empirical profile to desired_intensity
    emp_j = charac_data(:, j);
    emp_j(emp_j < 0) = 0;
    if max(emp_j) > 0
        emp_j = emp_j * (desired_intensity / max(emp_j));
    end
    emp_all(:, j) = emp_j;

    % Optimisation range: non-NaN positions
    valid = ~isnan(emp_j) & (emp_j >= 0);
    if any(valid)
        opt_limits_all(j, :) = [min(dist_bowl(valid)), max(dist_bowl(valid))];
    else
        opt_limits_all(j, :) = [dist_bowl(1), dist_bowl(end)];
    end
end

% Override opt_limits from config if provided
if isfield(parameters.calibration, 'opt_limits') && ~isempty(parameters.calibration.opt_limits)
    opt_limits_all = repmat(parameters.calibration.opt_limits(:)', n_depths, 1);
end

if isfield(parameters.calibration, 'opt_weights') && parameters.calibration.opt_weights ~= 0
    profile_weights = parameters.calibration.opt_weights;
else
    profile_weights = 0;
end

%% Objective: sum of per-depth profile errors + L2 regularisation
    function total_err = joint_objective(p)
        delta  = p(1:n_elem);
        vel    = p(end);
        total_err = 0;
        for jj = 1:n_depths
            phases_jj = geo_phases_all(jj, :) + delta;   % absolute phases for depth jj
            I_jj = eval_fwd_local(parameters, phases_jj, vel, dist_bowl, forward_model, water);
            lim  = opt_limits_all(jj, :);
            mask = dist_bowl >= lim(1) & dist_bowl <= lim(2);
            emp_jj = emp_all(:, jj);

            % Weight vector
            if profile_weights == 0
                w = ones(size(dist_bowl));
            else
                [~, ci] = get_flhm_center_position(dist_bowl, emp_jj);
                sigma = dist_bowl(ci) / profile_weights;
                w = normpdf(dist_bowl, dist_bowl(ci), sigma);
            end
            w = w / sum(w(mask));

            diff2 = (I_jj - emp_jj).^2 .* w;
            total_err = total_err + mean(diff2(mask));
        end
        total_err = total_err / n_depths;
        % L2 regularisation on magnitude of delta phases
        if lambda > 0
            total_err = total_err + lambda * mean(delta.^2);
        end
    end

%% Initial guess and bounds
if isfield(parameters.calibration, 'initial_velocity') && ...
        ~isempty(parameters.calibration.initial_velocity)
    v0 = parameters.calibration.initial_velocity;
else
    v0 = 0.05;
end

if ~isfield(parameters.calibration, 'opt_upper_velocity') || ...
        isempty(parameters.calibration.opt_upper_velocity)
    parameters.calibration.opt_upper_velocity = 0.2;
end

% Start from zero delta (= geometric steering), like BabelBrain
initial_guess = [zeros(1, n_elem), v0];
lower_bounds  = [-pi * ones(1, n_elem), 0.001];
upper_bounds  = [ pi * ones(1, n_elem), parameters.calibration.opt_upper_velocity];

if isfield(parameters.calibration, 'opt_seed')
    rng(parameters.calibration.opt_seed, 'twister');
end

%% Optimise
fprintf('Optimising %d delta phases + velocity across %d depths...\n', n_elem, n_depths);

if ~isfield(parameters.calibration, 'opt_method') || ...
        strcmp(parameters.calibration.opt_method, 'FEXminimize')
    opt_options = setoptimoptions('popsize', 2000, 'FinDiffType', 'central', 'TolCon', 1e-8);
    [p_opt, min_err] = minimize(@joint_objective, initial_guess, [], [], [], [], ...
        lower_bounds, upper_bounds, [], opt_options);
elseif strcmp(parameters.calibration.opt_method, 'GlobalSearch')
    gs = GlobalSearch;
    prob = createOptimProblem('fmincon', ...
        'x0', initial_guess, 'objective', @joint_objective, ...
        'lb', lower_bounds, 'ub', upper_bounds, ...
        'options', optimoptions('fmincon', 'OptimalityTolerance', 1e-8));
    [p_opt, min_err] = run(gs, prob);
else
    error('perform_joint_depth_fit: unknown opt_method ''%s''', ...
        parameters.calibration.opt_method);
end

opt_delta_phases = p_opt(1:n_elem);
opt_velocity     = p_opt(end);

fprintf('Joint fit complete. delta_phases [deg]: %s\n', ...
    mat2str(round(opt_delta_phases / pi * 180)));
fprintf('Velocity: %.4f m/s  Error: %.6f\n', opt_velocity, min_err);

%% Save diagnostic figure
img_folder = fullfile(parameters.io.outputs_folder, 'img_calibration');
if ~exist(img_folder, 'dir'); mkdir(img_folder); end

figure('Visible', 'off');
hold on;
colors = lines(n_depths);
for j = 1:n_depths
    phases_j = geo_phases_all(j, :) + opt_delta_phases;
    I_j = eval_fwd_local(parameters, phases_j, opt_velocity, dist_bowl, forward_model, water);
    plot(dist_bowl - dist_bowl_offset, emp_all(:, j), '-', 'Color', colors(j,:));
    plot(dist_bowl - dist_bowl_offset, I_j, '--', 'Color', colors(j,:));
end
xlabel('Depth from exit plane [mm]');
ylabel('Intensity [W/cm²]');
title('Joint depth fit: empirical (solid) vs model (dashed)');
legend(arrayfun(@(d) sprintf('%d mm', d), available_depths_ep, 'UniformOutput', false), ...
    'Location', 'best');
saveas(gcf, fullfile(img_folder, 'JointDepthFit.png'));
close(gcf);

end

%% Local helpers -----------------------------------------------------------
function I = eval_fwd_local(params, phases_rad, velocity, dist_bowl_mm, fwd_model, water)
    if strcmp(fwd_model, 'rayleigh')
        I = rayleigh_axial_intensity(phases_rad, velocity, params, dist_bowl_mm);
    else
        p = focusedAnnulusONeil( ...
            params.transducer.annular.curv_radius_mm / 1e3, ...
            [params.transducer.annular.elem_id_mm; params.transducer.annular.elem_od_mm] / 1e3, ...
            repmat(velocity, 1, params.transducer.annular.elem_n), ...
            phases_rad, ...
            params.transducer.freq_hz, ...
            water.sound_speed, ...
            water.density, ...
            (dist_bowl_mm - 0.5) * 1e-3);
        I = p.^2 / (2 * water.sound_speed * water.density) * 1e-4;
    end
end

function v = get_opt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
