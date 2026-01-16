function pct_skullmapping(subject_id, base_path)

% PCT_SKULLMAPPING Computes pseudo-CT mapping for cortical and trabecular bone using UTE histograms.
%
% This function processes UTE (ultrashort echo time) images and tissue masks 
% to compute pseudo-CT mappings for cortical and trabecular bone. It identifies 
% key points such as maxima, FWHM widths, and crossover points between cortical 
% and trabecular histograms. The function also generates visualizations and 
% exports the mapping formula to a text file.
%
% Input:
%   subject_id - String specifying the subject ID.
%   base_path  - String specifying the base path to the subject's m2m folder.
%
% Output:
%   None. The function saves visualizations and mapping coefficients in the subject's folder.

    %% Load UTE image and tissue mask
    subject_folder = fullfile(base_path, sprintf('m2m_sub-%03s', subject_id));
    ute_file = fullfile(subject_folder, 'UTE_reg_thr0_corr_norm.nii.gz');
    ute_corr = niftiread(ute_file);

    % Hard-coded SimNIBS labels for tissue types
    idx_cortical = 7; % Cortical bone
    idx_trabecular = 8; % Trabecular bone
    idx_air = 0; % Air

    % Load tissue mask
    mask_file = fullfile(subject_folder, 'final_tissues.nii.gz');
    mask_corr = niftiread(mask_file);

    %% Extract UTE values for each tissue type
    mask_cortical = logical(mask_corr == idx_cortical);
    ute_cortical = ute_corr(mask_cortical);

    mask_trabecular = logical(mask_corr == idx_trabecular);
    ute_trabecular = ute_corr(mask_trabecular);

    mask_air = logical(mask_corr == idx_air);
    ute_air = ute_corr(mask_air);

    %% Compute histograms for UTE values
    x = [0:0.01:1.5]; % Bin edges for histograms
    y1 = hist(ute_cortical, x); % Cortical bone histogram
    y2 = hist(ute_trabecular, x); % Trabecular bone histogram
    y3 = hist(ute_air, x); % Air histogram

    %% Interpolate values between 0.90 and 1.10 to exclude partial volumes
    interpPoints = [0.90:0.01:1.1];
    idx_interp = find(x >= interpPoints(1) & x <= interpPoints(end));
    y1(idx_interp) = NaN;
    y2(idx_interp) = NaN;
    y1(idx_interp) = interp1(x, y1, interpPoints, 'spline');
    y2(idx_interp) = interp1(x, y2, interpPoints, 'spline');

    %% Normalize histograms to their maximum values
    y1 = y1 ./ max(y1);
    y2 = y2 ./ max(y2);
    y3 = y3 ./ max(y3);

    %% Step 1: Find maxima of cortical and trabecular histograms
    [pks1, locs1] = findpeaks(y1, x); % Peaks for cortical histogram
    [pks2, locs2] = findpeaks(x .* y2, x); % Peaks for trabecular histogram (scaled by x)
    
    [~, idxMax1] = max(pks1);
    [~, idxMax2] = max(pks2);
    
    xMax1 = locs1(idxMax1); % Cortical peak position
    xMax2 = locs2(idxMax2); % Trabecular peak position
    
    disp(['Maxima of cortical histogram at x = ', num2str(xMax1)]);
    disp(['Maxima of trabecular histogram at x = ', num2str(xMax2)]);

    %% Step 2: Calculate FWQM for y1 (x1, cortical bone)

    % Find peak value and quarter-maximum value for y1
    [~, peakIdxY1] = max(y1);
    peakValueY1 = y1(peakIdxY1);
    halfMaxY1 = peakValueY1 / 4;
    
    % Find indices where y1 crosses half-maximum on both sides of the peak
    leftIdxY1 = find(y1(1:peakIdxY1) < halfMaxY1, 1, 'last');
    rightIdxY1 = find(y1(peakIdxY1:end) < halfMaxY1, 1, 'first') + peakIdxY1 - 1;
    
    % Interpolate to find exact crossing points
    xLeftY1 = interp1(y1(leftIdxY1:leftIdxY1+1), x(leftIdxY1:leftIdxY1+1), halfMaxY1);
    xRightY1 = interp1(y1(rightIdxY1-1:rightIdxY1), x(rightIdxY1-1:rightIdxY1), halfMaxY1);
    
    fwhmWidthY1 = xRightY1 - xLeftY1; % Calculate FWHM width
    
    disp(['FWHM width for y1: ', num2str(fwhmWidthY1)]);
    disp(['Left crossing point (xLeft): ', num2str(xLeftY1)]);
    disp(['Right crossing point (xRight): ', num2str(xRightY1)]);

    %% Step 3: Find crossover point (x1) within restricted range

    % Restrict range to [xMax1, xMax2]
    rangeIdx = (x >= xMax1 & x <= xMax2); % Logical index for restricted range
    xRestricted = x(rangeIdx);
    y1Restricted = y1(rangeIdx);
    y2Restricted = y2(rangeIdx);

    diffHistRestricted = y1Restricted - y2Restricted; % Difference between histograms
    % Define the window size for consecutive samples
    windowSize = 5;
    % Preallocate array for moving averages
    movingAvgDiff = zeros(length(diffHistRestricted) - windowSize + 1, 1);
    % Compute moving average of differences over the specified window size
    for i = 1:length(movingAvgDiff)
        movingAvgDiff(i) = mean(diffHistRestricted(i:i+windowSize-1));
    end
    % Find the index of the minimum absolute moving average difference
    [minDiff, minIdx] = min(abs(movingAvgDiff));
    % Adjust index to account for window size
    crossIdx = minIdx + floor(windowSize / 2);

    % Display results
    fprintf('Min. difference: %.2f\n', minDiff);
    fprintf('Crossover index: %d\n', crossIdx);

    % Check if crossIdx is empty
    if isempty(crossIdx)
        warning('No crossover point found within the restricted range. Choosing the midway point between tissue peaks. Visual inspection recommended.');
        xCross = round(mean([xMax1, xMax2]));
    else
        xCross = x(crossIdx + find(rangeIdx, 1, 'first')-1);
    end
    
    disp(['Crossover point (x1): ', num2str(xCross)]);
    
    %% Step 3: Calculate FWHM for y2 (x2, trabecular bone)
    % Find peak value and half-maximum value for y2
    [~, peakIdxY2] = max(y2);
    peakValueY2 = y2(peakIdxY2);
    halfMaxY2 = peakValueY2 / 4;
    
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

    % Define the three points
    X = [xLeftY1, xCross, xRightY2]; % X-coordinates (cortical lower FWHM; crossover; trabecular higher FWQM)
    Y = [1900, 700, 300];            % Y-coordinates
    
    % Solve for least squares fit using polyfit (degree 1 for linear fit)
    coefficients = polyfit(X, Y, 1);
    
    % Extract slope (m) and intercept (b)
    m = coefficients(1); % Slope
    b = coefficients(2); % Intercept
    
    % Display the linear equation
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
    saveas(h, fullfile(subject_folder, 'pseudoCT', figureName), 'png');
    
    %% Export values to text
    
    % saves coefficients in txt file in m2m folder
    fprintf(['Mapping formula: y = ', num2str(m), ' * x + ', num2str(b)]);

    txt_file = fopen(fullfile(subject_folder, 'pseudoCT', 'pCT_skull_mapping.txt'), 'w');
    fprintf(txt_file, '%.2f\n %.2f\n', round(m), round(b));
    fclose(txt_file);

    exit;
end