function pct_skullmapping(simnibs_path, path_pct)
% PCT_SKULLMAPPING  Compute pseudo-CT mapping for cortical and trabecular bone from UTE histograms
%
% Fits the subject-specific UTE->HU line from the cortical and trabecular
% UTE histograms (PCT_SKULLMAPPING_FIT, default settings) and saves the
% coefficients (pCT_skull_mapping.txt, rounded) and a figure to path_pct.
%
% Use as:
%   pct_skullmapping(simnibs_path, path_pct)
%
% Input:
%   simnibs_path - path to the subject's m2m folder (e.g. .../m2m_sub-009)
%   path_pct     - path to the pseudoCT working directory
%
% See also: PCT_SKULLMAPPING_FIT, PCT_SOFT_TISSUE_PEAK, PCT_SKULLEXPAND

arguments
    simnibs_path (1,:) char
    path_pct     (1,:) char
end

    %% Load UTE image and tissue mask
    ute_corr  = niftiread(fullfile(path_pct, 'UTE_STnorm.nii.gz'));
    mask_corr = niftiread(fullfile(simnibs_path, 'final_tissues.nii.gz'));

    %% Fit (pct_skullmapping_fit, default settings)
    [m, b, fi] = pct_skullmapping_fit(ute_corr, mask_corr);
    x = fi.x; y1 = fi.y_cortical; y2 = fi.y_trabecular;
    xLeftY1 = fi.x_points(1); xCross = fi.x_points(2); xRightY2 = fi.x_points(3);
    y3 = hist(ute_corr(mask_corr == 0), x); y3 = y3 ./ max(y3);    % air, for the plot only
    fprintf('Points: cortical lower width %.3f, crossover %.3f, trabecular upper width %.3f\n', fi.x_points);
    disp(['Least Squares Linear Fit: y = ', num2str(m), ' * x + ', num2str(b)]);

    %% Visualize the histograms with key points

    h = figure('Position', [100, 100, 600, 500]);
    subplot(3,1,[1,2])
        xx = linspace(0,0.1,2000);
        yy = 0.*xx-1000;
        xx1 = linspace(0.1,1,8000);
        pCT_kosciessa = m.*xx1+b;
        % Calculate alternative mappings
        pCT_miscouridou = -2085 * xx1 + 2329;
        pCT_carpino = -2194 * xx1 + 2236;
        pCT_wiesinger = -2000 * (xx1 - 1) + 42; % Adjusted for (UTE-1)
        pCT_treeby = -2929.6 * xx1 + 3274.9;
        xx2 = linspace(1,1.5,5000);
        yy2 = 0.*xx2+42;
        plot(xx,yy, 'LineWidth', 2, 'Color', [.5 .5 .5]);
        hold on
        p1 = plot(xx1,pCT_kosciessa, 'LineWidth', 2, 'Color', 'r');
        p2 = plot(xx1,pCT_miscouridou, 'LineWidth', 1);
        p3 = plot(xx1,pCT_carpino, 'LineWidth', 1);
        p4 = plot(xx1,pCT_wiesinger, 'LineWidth', 1);
        p5 = plot(xx1,pCT_treeby, 'LineWidth', 1);
        hold on
        plot(xx2,yy2, 'LineWidth', 2, 'Color', [.5 .5 .5]);
        % Add vertical lines with specific y-limits
        line([0.1, 0.1], [-1000, pCT_kosciessa(1)], 'Color', [.5 .5 .5], 'LineStyle', '-', 'LineWidth', 1.5);
        line([1, 1], [42, pCT_kosciessa(end)], 'Color', [.5 .5 .5], 'LineStyle', '-', 'LineWidth', 1.5);
        xlabel('UTE [a.u.]')
        ylabel('pCT [HU]')
        legend([p1, p2, p3, p4, p5], ...
            {'kosciessa'; 'miscouridou'; 'carpino'; 'wiesinger'; 'treeby'}, ...
            'location', 'NorthEast', 'orientation', 'vertical')
        legend('boxoff')
        title(['pCT = ', num2str(round(m)), ' UTE + ', num2str(round(b))])
    subplot(3,1,3)
        hold on;
        p1 = plot(x,y1, 'LineWidth', 2, 'Color', 'k');
        p2 = plot(x,y2, 'LineWidth', 2, 'Color', 'r');
        p3 = plot(x,y3, 'LineWidth', 2, 'Color', [.5 .5 .5]);
        % Plot vertical lines for crossover and FWQM points
        plot([xLeftY1 xLeftY1], [0 1], 'k:', 'LineWidth', 1);
        plot([xCross xCross], [0 1], 'k:', 'LineWidth', 1);
        plot([xRightY2 xRightY2], [0 1], 'k:', 'LineWidth', 1);
        xlabel('UTE'); ylabel({'Counts';'(max-norm.)'})
        legend([p1, p2, p3], {'cortical'; 'trabecular'; 'air'}, 'location', 'NorthEast')
        legend('boxoff')

    figureName = ['pCT_skull_mapping'];
    saveas(h, fullfile(path_pct, figureName), 'png');

    %% Export values to text
    fprintf(['Mapping formula: y = ', num2str(m), ' * x + ', num2str(b), '\n']);

    txt_file = fopen(fullfile(path_pct, 'pCT_skull_mapping.txt'), 'w');
    fprintf(txt_file, '%.2f\n %.2f\n', round(m), round(b));
    fclose(txt_file);

end