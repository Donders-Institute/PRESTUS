function pct_skullmapping(subject_id, base_path)
    
    subject_folder = fullfile(base_path, sprintf('m2m_sub-%03s', subject_id));
    ute_file = fullfile(subject_folder, 'UTE_reg_thr0_corr_norm.nii.gz');
    ute_corr = niftiread(ute_file);

    % Note: hard-coded to simnibs labels
    idx_cortical = 7;
    idx_trabecular = 8;
    idx_air = 0;

    % constrain to cortical and trabecular bone
    mask_file = fullfile(subject_folder, 'final_tissues.nii.gz');
    mask_corr = niftiread(mask_file);

    mask_cortical = logical(mask_corr==idx_cortical);
    ute_cortical = ute_corr(mask_cortical);

    mask_trabecular = logical(mask_corr==idx_trabecular);
    ute_trabecular = ute_corr(mask_trabecular);

    mask_air = logical(mask_corr==idx_air);
    ute_air = ute_corr(mask_air);

    x = []; y1 = []; y2 = [];
    x = [0:0.01:1.5];
    y1 = hist(ute_cortical, x);
    y2 = hist(ute_trabecular, x);
    y3 = hist(ute_air, x);

    % For skull distributions, interpolate values between 0.95 and 1.05 
    % (exclude potential soft tissue partial volumes)

    interpPoints = [0.90:0.01:1.1];
    idx_interp = find(x>=interpPoints(1) & x<=interpPoints(end));
    y1(idx_interp) = NaN;
    y2(idx_interp) = NaN;
    y1(idx_interp) = interp1(x, y1, interpPoints, 'spline');
    y2(idx_interp) = interp1(x, y2, interpPoints, 'spline');

    % normalize to maximum

    y1 = y1./max(y1);
    y2 = y2./max(y2);
    y3 = y3./max(y3);

    %% Step 1: Find maxima of both histograms
    [pks1, locs1] = findpeaks(y1, x);
    [pks2, locs2] = findpeaks(y2, x);
    
    % Identify global maxima for each histogram
    [maxPeak1, idxMax1] = max(pks1);
    [maxPeak2, idxMax2] = max(pks2);
    
    xMax1 = locs1(idxMax1); % x-coordinate of max for y1
    xMax2 = locs2(idxMax2); % x-coordinate of max for y2
    
    disp(['Maxima of y1 at x = ', num2str(xMax1)]);
    disp(['Maxima of y2 at x = ', num2str(xMax2)]);
    
    %% Step 2: Restrict range to [xMax1, xMax2]
    rangeIdx = (x >= xMax1 & x <= xMax2); % Logical index for restricted range
    xRestricted = x(rangeIdx);
    y1Restricted = y1(rangeIdx);
    y2Restricted = y2(rangeIdx);
    
    %% Step 3: Find crossover point (x1) within restricted range
    diffHistRestricted = y1Restricted - y2Restricted; % Difference between histograms
    crossIdx = find(diffHistRestricted(1:end-1) .* diffHistRestricted(2:end) < 0); % Detect sign change
    
    % Check if crossIdx is empty
    if isempty(crossIdx)
        warning('No crossover point found within the restricted range. Choosing the midway point between tissue peaks. Visual inspection recommended.');
        xCross = round(mean([xMax1, xMax2]));
    else
        xCross = interp1(diffHistRestricted(crossIdx:crossIdx+1), ...
                     xRestricted(crossIdx:crossIdx+1), 0);
    end
    
    disp(['Crossover point (x1): ', num2str(xCross)]);
    
    %% Step 3: Calculate FWHM for y2 (x2)
    % Find peak value and quarter-maximum value for y2
    [peakValueY2, peakIdxY2] = max(y2);
    halfMaxY2 = peakValueY2 / 2;
    
    % Find indices where y2 crosses half-maximum on both sides of the peak
    leftIdxY2 = find(y2(1:peakIdxY2) < halfMaxY2, 1, 'last');
    rightIdxY2 = find(y2(peakIdxY2:end) < halfMaxY2, 1, 'first') + peakIdxY2 - 1;
    
    % Interpolate to find exact crossing points
    xLeftY2 = interp1(y2(leftIdxY2:leftIdxY2+1), x(leftIdxY2:leftIdxY2+1), halfMaxY2);
    xRightY2 = interp1(y2(rightIdxY2-1:rightIdxY2), x(rightIdxY2-1:rightIdxY2), halfMaxY2);
    
    fwhmWidthY2 = xRightY2 - xLeftY2; % Calculate FWHM width
    
    disp(['FWHM width for y2: ', num2str(fwhmWidthY2)]);
    disp(['Left crossing point (xLeft): ', num2str(xLeftY2)]);
    disp(['Right crossing point (xRight): ', num2str(xRightY2)]);
    
    %% Assemble points for fitting
    
    X_CROSSOVER = xCross;
    X_FWQM = xRightY2;
    
    % assumptions about HU values at key points
    Y_CROSSOVER = 700;
    Y_FWQM = 300;
    
    % Calculate slope (m)
    m = (Y_FWQM - Y_CROSSOVER) / (X_FWQM - X_CROSSOVER);
    
    % Calculate intercept (c)
    c = Y_CROSSOVER - m * X_CROSSOVER;
    
    % Display the linear equation
    disp(['Linear equation: y = ', num2str(m), ' * x + ', num2str(c)]);
    
    %% Visualize the histograms with key points

    h = figure('Position', [100, 100, 600, 500]);
    subplot(3,1,[1,2])
        xx = linspace(0,0.1,2000);
        yy = 0.*xx-1000;
        xx1 = linspace(0.1,1,8000);
        pCT_kosciessa = m.*xx1+c;
        % Calculate alternative mappings
        pCT_miscouridou = -2085 * xx1 + 2329;
        pCT_carpino = -2194 * xx1 + 2236;
        pCT_wiesinger = -2000 * (xx1 - 1) + 42; % Adjusted for (UTE-1)
        pCT_treeby = -2929.6 * xx1 + 3247.9;
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
        title(['pCT = ', num2str(round(m)), ' UTE + ', num2str(round(c))])
    subplot(3,1,3)
        hold on;
        p1 = plot(x,y1, 'LineWidth', 2, 'Color', 'k');
        p2 = plot(x,y2, 'LineWidth', 2, 'Color', 'r');
        p3 = plot(x,y3, 'LineWidth', 2, 'Color', [.5 .5 .5]);
        % Plot vertical lines for crossover and FWQM points
        plot([X_CROSSOVER X_CROSSOVER], [0 1], 'k:', 'LineWidth', 1);
        plot([X_FWQM X_FWQM], [0 1], 'k:', 'LineWidth', 1);
        xlabel('UTE'); ylabel({'Counts';'(max-norm.)'})
        legend([p1, p2, p3], {'cortical'; 'trabecular'; 'air'}, 'location', 'NorthEast')
        legend('boxoff')

    figureName = ['pCT_mapping_skull'];
    saveas(h, fullfile(subject_folder, 'pseudoCT_backup', figureName), 'png');
    
    %% Export values to text
    
    % saves coefficients in txt file in m2m folder
    fprintf(['Mapping formula: y = ', num2str(m), ' * x + ', num2str(c)]);

    txt_file = fopen(fullfile(subject_folder, 'pct_skull_mapping.txt'), 'w');
    fprintf(txt_file, '%.2f\n %.2f\n', round(m), round(c));
    fclose(txt_file);

    exit;
end