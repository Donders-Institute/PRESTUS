function thermal_plot_protocol(params_thermal, parameters, varargin)
% THERMAL_PLOT_PROTOCOL  Plot the full ultrasound stimulation protocol timeline
%
% Generates a two-panel figure: the top panel shows a single pulse-train
% period with sine-wave ON intervals and shaded OFF regions; the bottom
% panel shows rectangular patches for each PT repetition across the full
% PTRD duration (up to 20 equidistant patches are drawn when many
% repetitions exist). The figure is saved as a PNG to
% parameters.io.dir_output and closed.
%
% Use as:
%   thermal_plot_protocol(params_thermal, parameters)
%   thermal_plot_protocol(params_thermal, parameters, Name, Value, ...)
%
% Input:
%   params_thermal - struct from THERMAL_PARAMETERS; must contain
%                    pt_dt [s], pt_on_steps_n, pt_on_steps_dur [s],
%                    pt_off_steps_n, pt_off_steps_dur [s],
%                    n_pulses_per_pt, n_ptri_reps, ptri_off_step_dur [s],
%                    ptri_off_steps_n, post_ptri_step_dur [s],
%                    post_ptri_steps_n, pd [s], pri [s], ptd [s],
%                    ptri [s], ptrd [s], dc, prf [Hz]
%   parameters     - PRESTUS config; must contain io.dir_output, io.output_affix,
%                    subject_id, simulation.medium
%
% Optional name-value pairs:
%   PulseFreqColor - (1,:) char, color for pulse ON regions (default: 'r')
%   OffColor       - (1,3) double, RGB color for OFF shading (default: [0.8 0.8 0.8])
%   FigureSize     - (1,2) double, figure size in inches (default: [14 10])
%   UpsampleFactor - (1,1) double, sine points per timestep (default: 50)
%
% See also: THERMAL_PARAMETERS, THERMAL_SIMULATION

p = inputParser;
addParameter(p, 'PulseFreqColor', 'r', @ischar);
addParameter(p, 'OffColor', [0.8 0.8 0.8], @(x) isnumeric(x) && length(x)==3);
addParameter(p, 'FigureSize', [14 10], @(x) length(x)==2);
addParameter(p, 'UpsampleFactor', 50, @isscalar);
parse(p, varargin{:});

pulse_color = p.Results.PulseFreqColor;
off_color = p.Results.OffColor;
fig_size = p.Results.FigureSize;
upsample = p.Results.UpsampleFactor;

if ~isfield(params_thermal, 'post_ptri_dur')
    params_thermal.post_ptri_dur = params_thermal.post_ptri_step_dur * params_thermal.post_ptri_steps_n;
end

CYCLES_PER_PD = 10;

function [t_high, signal_high] = generate_segment(t_start, total_duration, dt_coarse, is_on, upsample)
    n_steps = round(total_duration / dt_coarse);
    dt_high = total_duration / (n_steps * upsample);
    
    t_high = [];
    signal_high = [];
    t_local = 0;
    
    for step_i = 1:n_steps
        t_rel_step = (0:dt_high:dt_coarse-dt_high);
        t_abs_step = t_start + t_local + t_rel_step;
        if is_on
            phase = 2*pi * CYCLES_PER_PD * (t_local + t_rel_step) / total_duration;
            signal_step = sin(phase);
        else
            signal_step = zeros(size(t_rel_step));
        end
        t_high = [t_high; t_abs_step'];
        signal_high = [signal_high; signal_step'];
        t_local = t_local + dt_coarse;
    end
    
    t_high(end) = t_start + total_duration;
    if is_on
        signal_high(end) = sin(2*pi * CYCLES_PER_PD);
    end
end

%% Prepare Plot 1: Pulse train

total_on_dur  = params_thermal.pt_on_steps_n  * params_thermal.pt_on_steps_dur;
total_off_dur = params_thermal.pt_off_steps_n * params_thermal.pt_off_steps_dur;

t_first_pulse = []; signal_first_pulse = []; t_current = 0;

[t_on, s_on] = generate_segment(t_current, total_on_dur, params_thermal.pt_dt, true, upsample);
t_first_pulse = [t_first_pulse; t_on]; signal_first_pulse = [signal_first_pulse; s_on];
t_current = t_current + total_on_dur;

if params_thermal.pt_off_steps_n > 0 && params_thermal.pt_off_steps_dur>0
    [t_off, s_off] = generate_segment(t_current, total_off_dur, params_thermal.pt_dt, false, upsample);
    t_first_pulse = [t_first_pulse; t_off]; signal_first_pulse = [signal_first_pulse; s_off];
end

%% Prepare Plot 2: Patches for ON period in pulse repetitions

% max. 20 equi-distant ON patches across full PTRD
MAX_PATCHES = 20;
total_ptri_reps = params_thermal.n_ptri_reps;

if total_ptri_reps <= MAX_PATCHES
    % All patches
    step_reps = 1;
    n_ptri_to_plot = total_ptri_reps;
else
    % Equi-distant: select MAX_PATCHES evenly spaced
    step_reps = round(total_ptri_reps / MAX_PATCHES);
    n_ptri_to_plot = MAX_PATCHES;
end

t_on_patches = [];
t_current = 0;

for rep_i = 1:step_reps:total_ptri_reps
    if rep_i > total_ptri_reps, break; end
    t_on_start = (rep_i-1) * params_thermal.ptri;
    t_on_end = t_on_start + params_thermal.ptd;
    t_on_patches = [t_on_patches; t_on_start t_on_end];
end

% Full x-axis span
t_total_end = total_ptri_reps * params_thermal.ptri + params_thermal.post_ptri_dur;

%% Plot the thermal protocol

h = figure('Position', [100 100 fig_size(1)*100 fig_size(2)*100], 'Color', 'w');
% PLOT 1: First pulse (sine ON + flat OFF)
subplot(2,1,1);
    h1 = plot(t_first_pulse, signal_first_pulse, pulse_color, 'LineWidth', 2.5); hold on;
    
    % Shade pulse OFF only
    off_mask = abs(signal_first_pulse) < 1e-10;
    if any(off_mask)
        off_starts = find(diff([false; off_mask]) > 0);
        off_ends = find(diff([off_mask; false]) < 0);
        for i = 1:min(length(off_starts), length(off_ends))
            ts = t_first_pulse(off_starts(i):off_ends(i));
            fill(ts([1 end end 1]), [-1.2 -1.2 1.2 1.2], off_color, ...
                 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
    end
    
    title(sprintf('Pulse Detail: PRI=%.1fms'), 'FontSize', 13);
    ylabel('Signal'); grid on; ylim([-1.2 1.2]);
    legend(h1, 'US Pulse', 'Location', 'southeast');

% PLOT 2: Red ON patches (PT periods only) - NO SINE
subplot(2,1,2);
    hold on;
    
    % Plot patches
    hold on;
    for i = 1:size(t_on_patches, 1)
        t_patch = t_on_patches(i, :);
        fill([t_patch(1) t_patch(2) t_patch(2) t_patch(1)], [0 0 1 1], pulse_color, ...
             'FaceAlpha', 0.7, 'EdgeColor', pulse_color, 'LineWidth', 1);
    end

    % Title & axis
    title(sprintf('%d/%d PTs shown every %dth (%.1fs total + %.1fs post)', ...
        size(t_on_patches,1), total_ptri_reps, step_reps, params_thermal.ptrd, params_thermal.post_ptri_dur));
    xlim([0 t_total_end]);
    xlim([0 t_total_end]);
    sgtitle(sprintf('Thermal Protocol: PD=%.1fms/PRI=%.0fms/%d pulses-PT (%.1f%% DC)', ...
        params_thermal.pd*1e3, params_thermal.pri*1e3, params_thermal.n_pulses_per_pt, params_thermal.dc*100), ...
        'FontSize', 15, 'FontWeight', 'bold');

% save plot
output_plot_filename = fullfile(parameters.io.dir_img,...
    sprintf('sub-%03d_%s_thermal_protocol%s.png',...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
saveas(h, output_plot_filename, 'png')
close(h);

%% Print protocol summary

fprintf('Pulse: %.1fms ON + %.1fms OFF = %.1fms PRI ✓\n', ...
    total_on_dur*1e3, total_off_dur*1e3, (total_on_dur+total_off_dur)*1e3);
fprintf('PTRD: %d PT × %.2fs = %.1fs + %.1fs post ✓\n', params_thermal.n_ptri_reps, ...
    params_thermal.ptd, params_thermal.ptrd, params_thermal.post_ptri_dur);

end