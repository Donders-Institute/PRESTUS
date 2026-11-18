function thermal_plot_protocol(params_thermal, parameters, varargin)
% THERMAL_PLOT_PROTOCOL  Plot the sonication protocol of a thermal simulation
%
% Draws the protocol at three time scales, top to bottom: the whole session
% (every pulse train, then the steady-state period), one pulse-train interval
% (its pulses) and one pulse repetition interval. Each row shows the on/off
% envelope (not the carrier wave). The bottom row also marks the thermal
% simulation steps that resolve the pulse, and the middle row the coarse
% steps between trains, so the discretisation can be checked at a glance.
% With a THERMAL_SCHEDULE of several sources (async coupling), a fourth row
% shows one row of on-blocks per transducer over one cycle. The figure is
% saved as a PNG to parameters.io.dir_img and closed.
%
% Use as:
%   thermal_plot_protocol(params_thermal, parameters)
%   thermal_plot_protocol(params_thermal, parameters, Name, Value, ...)
%
% Input:
%   params_thermal - struct from THERMAL_PARAMETERS; must contain
%                    pd, pri, ptd, ptri, ptrd [s], dc, prf [Hz],
%                    n_pulses_per_pt, n_ptri_reps,
%                    pt_on_steps_n, pt_on_steps_dur [s],
%                    pt_off_steps_n, pt_off_steps_dur [s],
%                    ptri_off_steps_n, ptri_off_step_dur [s],
%                    post_ptri_steps_n, post_ptri_step_dur [s];
%                    optional ramp_dur [s], ramp_shape, ramp_env, ramp_eta
%   parameters     - PRESTUS config; must contain io.dir_img, io.output_affix,
%                    subject_id, simulation.medium
%
% Optional name-value pairs:
%   PulseFreqColor - colour of the on periods (default: [0.80 0.22 0.17])
%   OffColor       - (1,3) RGB of the off shading (default: [0.92 0.92 0.92])
%   FigureSize     - (1,2) figure size in inches (default: [11 8])
%   Schedule       - struct from THERMAL_SCHEDULE (default: []); with more
%                    than one source, adds the per-transducer row
%   UpsampleFactor - ignored; kept so existing calls still run
%
% See also: THERMAL_PARAMETERS, THERMAL_SIMULATION, THERMAL_SCHEDULE

p = inputParser;
addParameter(p, 'PulseFreqColor', [0.80 0.22 0.17], @(x) ischar(x) || (isnumeric(x) && numel(x) == 3));
addParameter(p, 'OffColor', [0.92 0.92 0.92], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'FigureSize', [11 8], @(x) numel(x) == 2);
addParameter(p, 'UpsampleFactor', 50, @isscalar);
addParameter(p, 'Schedule', [], @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
sched = p.Results.Schedule;
on_col = p.Results.PulseFreqColor;
off_col = p.Results.OffColor;
cool_col = [0.85 0.91 0.97];
fig_size = p.Results.FigureSize;

pt = params_thermal;
if ~isfield(pt, 'post_ptri_dur')
    pt.post_ptri_dur = pt.post_ptri_step_dur * pt.post_ptri_steps_n;
end
ptri_off = max(pt.ptri - pt.ptd, 0);
multi_source = ~isempty(sched) && isfield(sched, 'n_sources') && sched.n_sources > 1;
n_rows = 3 + multi_source;

% Envelopes are drawn as one patch per row. Blocks or gaps narrower than
% MIN_PX pixels would alias into stripes: such a row draws the span as one
% band in the duty-cycle-weighted colour and says so.
MIN_PX = 2;
row_px = fig_size(1) * 150 * 0.9;          % approximate plot width at the saved resolution

h = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
           'Position', [1 1 fig_size(1) fig_size(2)]);
tl = tiledlayout(h, n_rows, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Row 1: the whole session
ax = nexttile(tl);
t_end = pt.ptrd + pt.post_ptri_dur;
[sc, unit] = time_unit(t_end);
starts = (0:pt.n_ptri_reps-1)' * pt.ptri;
if pt.post_ptri_dur > 0
    shade(ax, pt.ptrd * sc, t_end * sc, cool_col);
end
if ptri_off > 0
    shade(ax, 0, pt.ptrd * sc, off_col);
end
what = plural(pt.n_ptri_reps, 'pulse train');
if resolvable(pt.ptd, ptri_off, t_end, row_px, MIN_PX)
    blocks(ax, starts * sc, (starts + pt.ptd) * sc, on_col);
else
    blocks(ax, 0, pt.ptrd * sc, band_colour(on_col, off_col, pt.ptd / pt.ptri));
    what = [what ' (drawn as a band)'];
end
style_row(ax, [0 t_end * sc], unit);
if pt.post_ptri_dur > 0
    text(ax, (pt.ptrd + pt.post_ptri_dur / 2) * sc, 0.5, sprintf('steady-state\n%s', fmt_t(pt.post_ptri_dur)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.25 0.35 0.5]);
end
if pt.post_ptri_dur > 0
    cool = sprintf(', then %s steady-state', fmt_t(pt.post_ptri_dur));
else
    cool = ', no steady-state period';
end
title(ax, sprintf('Session: %s over %s (PTRD)%s', what, fmt_t(pt.ptrd), cool), 'FontWeight', 'normal');

%% Row 2: one pulse-train interval
ax = nexttile(tl);
[sc, unit] = time_unit(pt.ptri);
p_starts = (0:pt.n_pulses_per_pt-1)' * pt.pri;
shade(ax, 0, pt.ptri * sc, off_col);
what = plural(pt.n_pulses_per_pt, 'pulse');
if resolvable(pt.pd, pt.pri - pt.pd, pt.ptri, row_px, MIN_PX)
    blocks(ax, p_starts * sc, (p_starts + pt.pd) * sc, on_col);
else
    blocks(ax, 0, pt.ptd * sc, band_colour(on_col, off_col, pt.dc));
    what = [what ' (drawn as a band)'];
end
style_row(ax, [0 pt.ptri * sc], unit);
if ptri_off > 0 && pt.ptri_off_steps_n > 0 && pt.ptri_off_steps_n <= 200
    step_ticks(ax, (pt.ptd + (1:pt.ptri_off_steps_n-1) * pt.ptri_off_step_dur) * sc);
end
if ptri_off > 0
    off_note = sprintf(', then %s off in %s of %s', fmt_t(ptri_off), ...
        plural(pt.ptri_off_steps_n, 'step'), fmt_t(pt.ptri_off_step_dur));
else
    off_note = ', no pause between trains';
end
title(ax, sprintf('One pulse train: %s in %s (PTD)%s; repeats every %s (PTRI)', ...
      what, fmt_t(pt.ptd), off_note, fmt_t(pt.ptri)), 'FontWeight', 'normal');

%% Row 3: one pulse repetition interval, with the simulation steps
ax = nexttile(tl);
[sc, unit] = time_unit(pt.pri);
shade(ax, 0, pt.pri * sc, off_col);
ramped = isfield(pt, 'ramp_dur') && pt.ramp_dur > 0;
if ramped
    % ramped envelope (THERMAL_PARAMETERS ramp_env); the heat follows its square
    tt = linspace(0, pt.pd, 400);
    patch(ax, [tt fliplr(tt)] * sc, [pt.ramp_env(tt) zeros(size(tt))], on_col, 'EdgeColor', 'none');
else
    blocks(ax, 0, pt.pd * sc, on_col);
end
style_row(ax, [0 pt.pri * sc], unit);
edges = [(0:pt.pt_on_steps_n) * pt.pt_on_steps_dur, ...
         pt.pd + (1:pt.pt_off_steps_n) * pt.pt_off_steps_dur];
if numel(edges) <= 400
    step_ticks(ax, edges(2:end-1) * sc);
end
if pt.pt_off_steps_n > 0
    steps = sprintf('%d on + %d off simulation steps', pt.pt_on_steps_n, pt.pt_off_steps_n);
else
    steps = sprintf('%d simulation %s', pt.pt_on_steps_n, ternary(pt.pt_on_steps_n == 1, 'step', 'steps'));
end
if ramped
    steps = sprintf('%s; %s ramp %s per edge (heat x %.3f)', steps, pt.ramp_shape, fmt_t(pt.ramp_dur), pt.ramp_eta);
end
title(ax, sprintf('One pulse: %s on (PD) every %s (PRI); %s', fmt_t(pt.pd), fmt_t(pt.pri), steps), ...
      'FontWeight', 'normal');

%% Row 4 (async): one row of on-blocks per transducer over one cycle
if multi_source
    ax = nexttile(tl);
    hold(ax, 'on');
    src = sched.source;
    switch sched.mode
        case 'pulses', t_win = 2 * src(1).pri;
        case 'trains', t_win = src(1).ptri;
        otherwise,     t_win = max([src.offset_s] + [src.ptri]);
    end
    t_win = min(t_win, sched.t_sonication);
    [sc, unit] = time_unit(t_win);
    cols = lines(sched.n_sources);
    for i = 1:sched.n_sources
        w = src(i).pulse_windows;
        w = w(w(:,1) < t_win, :);
        if size(w, 1) > 100            % too many pulses to draw: their trains
            w = src(i).train_windows;
            w = w(w(:,1) < t_win, :);
        end
        y0 = sched.n_sources - i;
        if ~isempty(w)
            blocks(ax, w(:,1) * sc, min(w(:,2), t_win) * sc, cols(i,:), y0 + [0.15 0.85]);
        end
    end
    yticks(ax, (0:sched.n_sources-1) + 0.5);
    yticklabels(ax, arrayfun(@(i) sprintf('T%d (+%s)', i, fmt_t(src(i).offset_s)), ...
        sched.n_sources:-1:1, 'UniformOutput', false));
    xlim(ax, [0 t_win * sc]); ylim(ax, [0 sched.n_sources]);
    xlabel(ax, sprintf('Time [%s]', unit));
    box(ax, 'off');
    title(ax, sprintf('Async schedule ''%s'': %d transducers, overall duty cycles %s %%', ...
        sched.mode, sched.n_sources, mat2str(round(100 * [src.duty_cycle], 3, 'significant'))), ...
        'FontWeight', 'normal');
end

title(tl, sprintf(['Sonication protocol (thermal simulation)\nPD %s · PRI %s · duty cycle %.3g%% · PRF %.4g Hz\n' ...
      '{\\fontsize{10}\\rm dotted lines: thermal simulation steps}'], ...
      fmt_t(pt.pd), fmt_t(pt.pri), pt.dc * 100, pt.prf), 'FontWeight', 'bold', 'FontSize', 13);

% save plot
output_plot_filename = fullfile(parameters.io.dir_img, ...
    sprintf('sub-%03d_%s_thermal_protocol%s.png', ...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
if exist('exportgraphics', 'file') == 2
    exportgraphics(h, output_plot_filename, 'Resolution', 150);
else
    saveas(h, output_plot_filename, 'png');
end
close(h);

end

%% ── helpers ─────────────────────────────────────────────────────────────
function blocks(ax, t0, t1, col, y)
% Rectangles [t0(k) t1(k)] x y as one patch.
if nargin < 5, y = [0 1]; end
t0 = t0(:); t1 = t1(:); n = numel(t0);
X = [t0 t1 t1 t0]'; Y = repmat([y(1) y(1) y(2) y(2)]', 1, n);
patch(ax, X, Y, col, 'EdgeColor', 'none');
hold(ax, 'on');
end

function shade(ax, t0, t1, col)
patch(ax, [t0 t1 t1 t0], [0 0 1 1], col, 'EdgeColor', 'none');
hold(ax, 'on');
end

function step_ticks(ax, t)
% Dotted lines across the envelope at simulation step boundaries.
if isempty(t), return; end
t = t(:)';
plot(ax, [t; t], repmat([0; 1], 1, numel(t)), ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.75);
end

function ok = resolvable(t_on, t_off, span, row_px, min_px)
% True when on blocks and the gaps between them are each at least min_px wide.
px = row_px / span;
ok = t_on * px >= min_px && (t_off <= 0 || t_off * px >= min_px);
end

function c = band_colour(on_col, off_col, frac)
% The mean colour of a row that is on for a fraction frac of the time.
if ischar(on_col), on_col = validatecolor(on_col); end
c = frac * on_col + (1 - frac) * off_col;
end

function style_row(ax, xl, unit)
xlim(ax, xl); ylim(ax, [0 1.05]);
yticks(ax, []); box(ax, 'off');
ax.XGrid = 'off';
xlabel(ax, sprintf('Time [%s]', unit));
end

function [sc, unit] = time_unit(t_max)
% Axis scale and unit for a span of t_max seconds.
if t_max < 1
    sc = 1e3; unit = 'ms';
elseif t_max <= 300
    sc = 1; unit = 's';
else
    sc = 1/60; unit = 'min';
end
end

function s = fmt_t(t)
% A duration in readable units.
if t <= 0
    s = '0 s';
elseif t < 1e-3
    s = sprintf('%.3g µs', t * 1e6);
elseif t < 1
    s = sprintf('%.3g ms', t * 1e3);
elseif t <= 300
    s = sprintf('%.3g s', t);
else
    s = sprintf('%.3g min', t / 60);
end
end

function s = plural(n, word)
s = sprintf('%d %s', n, ternary(n == 1, word, [word 's']));
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
