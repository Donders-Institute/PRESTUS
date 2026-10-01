function plot_calibration_report(prof, report, diag, title_str, png_path)
% PLOT_CALIBRATION_REPORT  Measured vs library vs calibrated axial profiles, one panel per depth
%
% Figure for the result of CALIBRATE_LIBRARY_ENTRY, saved as a PNG.
%
% Use as:
%   plot_calibration_report(prof, report, diag, title_str, png_path)
%
% Input:
%   prof      - profiles the entry was calibrated on (CALIBRATE_LIBRARY_ENTRY output)
%   report    - report table (r2_library, r2_refined per depth)
%   diag      - diag struct (.before from EVALUATE_LIBRARY_ENTRY, .z_after, .I_after)
%   title_str - figure title
%   png_path  - output PNG
%
% See also: CALIBRATE_LIBRARY_ENTRY, REFINE_TRANSDUCER_LIBRARY, CALIBRATION_TRANSDUCER

    before = diag.before;
    n  = numel(prof.depths_ep_mm);
    fh = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 min(1800, 330 * n) 320]);
    tl = tiledlayout(fh, 1, n, 'TileSpacing', 'compact');
    for k = 1:n
        ax = nexttile(tl); hold(ax, 'on');
        plot(ax, before.z_ep_mm, before.I_meas(:, k), 'k-', 'LineWidth', 1.6, 'DisplayName', 'Hydrophone');
        plot(ax, before.z_ep_mm, before.I_model(:, k), '-', 'Color', [0 0.45 0.70], ...
            'DisplayName', sprintf('Library (R² %.2f)', report.r2_library(k)));
        plot(ax, diag.z_after, diag.I_after(:, k), '-', 'Color', [0.85 0.37 0.01], 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Refined (R² %.2f)', report.r2_refined(k)));
        xline(ax, prof.depths_ep_mm(k), ':', 'HandleVisibility', 'off');
        ttl = sprintf('%.1f mm', prof.depths_ep_mm(k));
        if ~strcmp(prof.profile_type{k}, 'measured'); ttl = sprintf('%s (%s profile)', ttl, prof.profile_type{k}); end
        title(ax, ttl, 'FontWeight', 'normal');
        xlabel(ax, 'mm from exit plane'); grid(ax, 'on'); box(ax, 'on');
        if k == 1; ylabel(ax, 'Normalised intensity'); legend(ax, 'Location', 'northeast', 'FontSize', 7); end
    end
    title(tl, title_str);
    exportgraphics(fh, png_path, 'Resolution', 150); close(fh);
    fprintf('Report figure: %s\n', png_path);
end
