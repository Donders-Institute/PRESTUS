%% benchmark_phase_precession  Compare unconstrained vs. linear vs. monotonic phase optimisation
%
% Runs perform_global_search under each of the three opt_phase_precession
% modes on a shared synthetic target profile and reports optimisation error,
% wall-clock time, and the recovered phase pattern for each mode.
%
% The synthetic target is generated from the O'Neil model with a known
% ground-truth phase vector so that the true minimum error is approximately
% zero and recovery accuracy is assessable.
%
% Results are written to the console and saved as a figure and a MAT file
% in the current working directory.
%
% Usage:
%   Run from the PRESTUS root directory (so addpath calls resolve).
%   Adjust the CONFIGURATION section below to match your transducer geometry.

close all; clear;

%% Paths
func_path = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(func_path, 'functions')));
addpath(genpath(fullfile(func_path, 'external')));
addpath(genpath(fullfile(func_path, 'config')));

%% ---- CONFIGURATION -------------------------------------------------------
% Transducer geometry — edit to match your annular transducer
n_elem          = 4;            % number of annular elements
curv_radius_mm  = 64;           % bowl radius of curvature [mm]
freq_hz         = 500e3;        % driving frequency [Hz]
elem_id_mm      = [0, 14.5, 24.0, 30.5];    % inner diameters [mm]
elem_od_mm      = [14.0, 23.5, 30.0, 36.0]; % outer diameters [mm]

% Ground-truth phases used to synthesise the target profile
gt_phases_deg   = [0, 45, 90, 135];  % one per element [deg]
gt_velocity     = 0.05;              % particle velocity [m/s]

% Axial range for the synthetic profile
axial_mm        = (1:0.5:120)';  % [mm] from transducer surface

% Optimisation settings shared across all modes
opt_seed        = 42;
opt_upper_vel   = 0.2;
opt_method      = 'FEXminimize';  % 'FEXminimize' | 'GlobalSearch'
opt_weights     = 0;              % 0 = uniform weighting
opt_limits      = [5, 100];       % distance limits for error computation [mm]

% Number of repetitions per mode (for timing stability)
n_reps          = 3;
% ---------------------------------------------------------------------------

%% Build minimal parameters struct
parameters = struct();

parameters.transducer.annular.elem_n           = n_elem;
parameters.transducer.annular.curv_radius_mm   = curv_radius_mm;
parameters.transducer.annular.elem_id_mm       = elem_id_mm;
parameters.transducer.annular.elem_od_mm       = elem_od_mm;
parameters.transducer.freq_hz                  = freq_hz;

parameters.medium_properties.water.sound_speed = 1500; % m/s
parameters.medium_properties.water.density     = 1000; % kg/m³

parameters.calibration.opt_method         = opt_method;
parameters.calibration.opt_seed           = opt_seed;
parameters.calibration.opt_upper_velocity = opt_upper_vel;
parameters.calibration.opt_weights        = opt_weights;
parameters.calibration.opt_limits         = opt_limits;

% Dummy output folder (perform_global_search saves a figure there)
tmp_out = fullfile(tempdir, 'prestus_benchmark');
if ~exist(tmp_out, 'dir'); mkdir(tmp_out); end
parameters.io.outputs_folder = tmp_out;

%% Synthesise ground-truth target profile via O'Neil model
gt_phases_rad = gt_phases_deg / 180 * pi;

p_gt = focusedAnnulusONeil(...
    curv_radius_mm / 1e3, ...
    [elem_id_mm; elem_od_mm] / 1e3, ...
    repmat(gt_velocity, 1, n_elem), ...
    gt_phases_rad, ...
    freq_hz, ...
    parameters.medium_properties.water.sound_speed, ...
    parameters.medium_properties.water.density, ...
    (axial_mm - 0.5) * 1e-3);

profile_target.axial_distance_bowl = axial_mm;
profile_target.axial_intensity = p_gt.^2 / ...
    (2 * parameters.medium_properties.water.sound_speed * ...
     parameters.medium_properties.water.density) * 1e-4;

%% Define modes to benchmark
modes  = {false, 'linear', 'monotonic'};
labels = {'unconstrained', 'linear', 'monotonic'};
n_modes = numel(modes);

results(n_modes) = struct('mode', [], 'errors', [], 'times', [], ...
    'opt_phases_deg', [], 'opt_velocity', []);

%% Run benchmark
fprintf('\n=== Phase precession benchmark ===\n');
fprintf('Transducer: %d elements, %.0f kHz, Rc=%.0f mm\n', n_elem, freq_hz/1e3, curv_radius_mm);
fprintf('Ground truth phases: %s deg,  velocity: %.4f m/s\n\n', ...
    mat2str(gt_phases_deg), gt_velocity);

for m = 1:n_modes
    parameters.calibration.opt_phase_precession = modes{m};
    errs  = nan(1, n_reps);
    times = nan(1, n_reps);
    last_phases   = [];
    last_velocity = [];

    fprintf('--- Mode: %s ---\n', labels{m});
    for r = 1:n_reps
        t0 = tic;
        [opt_phases, opt_velocity, min_err] = perform_global_search(...
            parameters, profile_target, gt_velocity * 0.8); % start 20% below truth
        times(r) = toc(t0);
        errs(r)  = min_err;
        last_phases   = opt_phases;
        last_velocity = opt_velocity;
        fprintf('  rep %d/%d: error=%.6f, time=%.1f s\n', r, n_reps, min_err, times(r));
    end

    results(m).mode           = labels{m};
    results(m).errors         = errs;
    results(m).times          = times;
    results(m).opt_phases_deg = last_phases / pi * 180;
    results(m).opt_velocity   = last_velocity;
end

%% Print summary table
fprintf('\n%-16s  %10s  %10s  %10s  %10s\n', ...
    'Mode', 'mean_err', 'min_err', 'mean_t(s)', 'min_t(s)');
fprintf('%s\n', repmat('-', 1, 62));
for m = 1:n_modes
    fprintf('%-16s  %10.6f  %10.6f  %10.1f  %10.1f\n', ...
        results(m).mode, ...
        mean(results(m).errors), ...
        min(results(m).errors), ...
        mean(results(m).times), ...
        min(results(m).times));
end
fprintf('\nGround truth phases  : %s deg\n', mat2str(round(gt_phases_deg)));
for m = 1:n_modes
    fprintf('Recovered (%s)%s: %s deg\n', ...
        results(m).mode, ...
        repmat(' ', 1, max(0, 14-length(results(m).mode))), ...
        mat2str(round(results(m).opt_phases_deg)));
end

%% Plot
fig = figure('Name', 'Phase precession benchmark', 'Position', [100 100 900 600]);
colors = lines(n_modes);

% Error comparison
subplot(1,3,1);
hold on;
for m = 1:n_modes
    bh = bar(m, mean(results(m).errors));
    bh.FaceColor = colors(m,:);
    errorbar(m, mean(results(m).errors), std(results(m).errors), 'k', 'LineWidth', 1.2);
end
set(gca, 'XTick', 1:n_modes, 'XTickLabel', labels, 'XTickLabelRotation', 20);
ylabel('Mean optimisation error');
title('Optimisation error');
box off;

% Time comparison
subplot(1,3,2);
hold on;
for m = 1:n_modes
    bh = bar(m, mean(results(m).times));
    bh.FaceColor = colors(m,:);
    errorbar(m, mean(results(m).times), std(results(m).times), 'k', 'LineWidth', 1.2);
end
set(gca, 'XTick', 1:n_modes, 'XTickLabel', labels, 'XTickLabelRotation', 20);
ylabel('Wall-clock time (s)');
title('Computation time');
box off;

% Recovered phases vs. ground truth
subplot(1,3,3);
hold on;
plot(1:n_elem, mod(gt_phases_deg, 360), 'k--o', 'LineWidth', 1.5, 'DisplayName', 'ground truth');
for m = 1:n_modes
    plot(1:n_elem, mod(results(m).opt_phases_deg, 360), '-o', ...
        'Color', colors(m,:), 'LineWidth', 1.2, 'DisplayName', labels{m});
end
xlabel('Element index');
ylabel('Phase (deg)');
title('Recovered vs. ground-truth phases');
legend('Location', 'best');
ylim([0, 360]);
box off;

sgtitle(sprintf('Phase precession benchmark  |  %d elem, %.0f kHz, %d reps', ...
    n_elem, freq_hz/1e3, n_reps));

fig_path = fullfile(pwd, 'benchmark_phase_precession.png');
saveas(fig, fig_path);
fprintf('\nFigure saved: %s\n', fig_path);

%% Save results
mat_path = fullfile(pwd, 'benchmark_phase_precession.mat');
save(mat_path, 'results', 'gt_phases_deg', 'gt_velocity', 'profile_target', 'parameters');
fprintf('Results saved: %s\n', mat_path);
