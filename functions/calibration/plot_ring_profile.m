function fig = plot_ring_profile(tran, rp, info, out_file)
% PLOT_RING_PROFILE  Show a fitted within-ring velocity profile and what it changes
%
% Left: the bowl in cross-section, each sub-ring drawn with its weight
% (line width ~ |u|, colour = phase). Middle: |u| and phase against the
% radius, ring boundaries marked. Right (with info from FIT_RING_PROFILE):
% measured profiles against the uniform-ring and the profile model at up to
% six depths, with R² for each, and the R² of every depth for both models.
%
% Use as:
%   plot_ring_profile(tran, rp)
%   plot_ring_profile(tran, rp, info)
%   fig = plot_ring_profile(tran, rp, info, out_file)   % also saves a PNG
%
% Input:
%   tran     - equipment transducer struct (annular, without profile)
%   rp       - ring profile (FIT_RING_PROFILE, or lib.global_model.ring_profile)
%   info     - optional second output of FIT_RING_PROFILE ([] = profile only)
%   out_file - optional PNG path
%
% See also: FIT_RING_PROFILE, RING_PROFILE_TRANSDUCER

arguments
    tran (1,1) struct
    rp (1,1) struct
    info = []
    out_file {mustBeTextScalar} = ''
end

    if nargin < 3; info = []; end
    if nargin < 4; out_file = ''; end
    trp = ring_profile_transducer(tran, rp);
    a   = trp.transducer.annular; a0 = tran.transducer.annular;
    u   = a.elem_weight; R = a.curv_radius_mm;
    r_in = a.elem_id_mm / 2; r_out = a.elem_od_mm / 2; rc = (r_in + r_out) / 2;
    cm = phase_colormap();
    ph = mod(angle(u) * 180 / pi + 180, 360) - 180;
    col = @(p) cm(max(1, min(size(cm, 1), round((p + 180) / 360 * (size(cm, 1) - 1)) + 1)), :);

    has_info = ~isempty(info);
    fig = figure('Color', 'w', 'Position', [60 60 1400 ternary(has_info, 820, 420)]);
    tl = tiledlayout(fig, ternary(has_info, 2, 1), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Within-ring velocity profile: %d sub-rings per ring, smoothness %g', ...
        rp.n_sub, rp.smooth), 'FontWeight', 'bold');

    % Bowl cross-section
    ax = nexttile(tl); hold(ax, 'on'); axis(ax, 'equal');
    for e = 1:numel(u)
        s = linspace(r_in(e), r_out(e), 12);
        zb = R - sqrt(R^2 - s.^2);
        for sg = [-1 1]
            plot(ax, zb, sg * s, '-', 'Color', col(ph(e)), 'LineWidth', 1 + 5 * abs(u(e)) / max(abs(u)));
        end
    end
    for b = unique([a0.elem_id_mm, a0.elem_od_mm]) / 2
        zb = R - sqrt(R^2 - b^2);
        plot(ax, zb + [-1 1], [b b], 'k-', 'LineWidth', 0.5); plot(ax, zb + [-1 1], -[b b], 'k-', 'LineWidth', 0.5);
    end
    colormap(ax, cm); clim(ax, [-180 180]); cb = colorbar(ax); cb.Label.String = 'phase (°)';
    xlabel(ax, 'depth from apex (mm)'); ylabel(ax, 'radius (mm)');
    title(ax, 'Bowl: width = |u|, colour = phase');

    % |u| and phase against radius
    ax = nexttile(tl); hold(ax, 'on');
    yyaxis(ax, 'left');
    stairs(ax, [r_in, r_out(end)], [abs(u), abs(u(end))], 'LineWidth', 2);
    yline(ax, 1, ':'); ylabel(ax, '|u| (ring mean = 1)');
    yyaxis(ax, 'right');
    plot(ax, rc, ph, 'o-', 'LineWidth', 1.2); ylabel(ax, 'phase (°)');
    for b = a0.elem_od_mm / 2; xline(ax, b, '-', 'Color', [0.7 0.7 0.7]); end
    xlabel(ax, 'radius (mm)'); title(ax, 'Weights across each ring (grey: ring edges)');
    grid(ax, 'on');

    if has_info
        % R² per depth
        ax = nexttile(tl); hold(ax, 'on');
        d = info.depths_ep_mm;
        plot(ax, d, info.r2_uniform, 'o-', 'LineWidth', 1.5, 'DisplayName', 'uniform rings');
        plot(ax, d, info.r2_profile, 's-', 'LineWidth', 1.5, 'DisplayName', 'ring profile');
        if any(isfinite(info.r2_heldout))
            plot(ax, d, info.r2_heldout, 'd--', 'LineWidth', 1.2, 'DisplayName', 'ring profile, depth held out');
        end
        xlabel(ax, 'focal depth setting (mm)'); ylabel(ax, sprintf('R² (z ≥ %g mm)', rp.near_field_skip_mm));
        legend(ax, 'Location', 'best', 'Box', 'off'); grid(ax, 'on'); title(ax, 'Fit per depth');
        % Example depths
        [~, order] = sort(info.r2_profile - info.r2_uniform, 'descend');
        pick = sort(order(1:min(3, numel(order))));
        for k = pick
            ax = nexttile(tl); hold(ax, 'on');
            plot(ax, info.z_ep_mm, info.I_meas(:, k), 'k', 'LineWidth', 1.6, 'DisplayName', 'hydrophone');
            plot(ax, info.z_ep_mm, info.I_uniform(:, k), '--', 'LineWidth', 1.4, ...
                'DisplayName', sprintf('uniform rings (R² %.3f)', info.r2_uniform(k)));
            plot(ax, info.z_ep_mm, info.I_profile(:, k), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('ring profile (R² %.3f)', info.r2_profile(k)));
            xlabel(ax, 'z from exit plane (mm)'); ylabel(ax, 'I / I_{max}');
            title(ax, sprintf('%.1f mm setting', d(k))); legend(ax, 'Location', 'best', 'Box', 'off');
            grid(ax, 'on');
        end
    end
    if ~isempty(out_file)
        exportgraphics(fig, out_file, 'Resolution', 150);
    end
end

function cm = phase_colormap()
% Cyclic blue–white–orange–white–blue map for phase
    t = linspace(0, 1, 256)';
    b = [0.165 0.471 0.839]; o = [0.922 0.408 0.204]; w = [0.96 0.96 0.96];
    cm = zeros(256, 3);
    for k = 1:256
        x = t(k);
        if x < 0.25;      cm(k, :) = b + (w - b) * (x / 0.25);
        elseif x < 0.5;   cm(k, :) = w + (o - w) * ((x - 0.25) / 0.25);
        elseif x < 0.75;  cm(k, :) = o + (w - o) * ((x - 0.5) / 0.25);
        else;             cm(k, :) = w + (b - w) * ((x - 0.75) / 0.25);
        end
    end
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end
