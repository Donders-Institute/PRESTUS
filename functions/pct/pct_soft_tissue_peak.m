function pct_soft_tissue_peak(simnibs_folder, path_pct)
% PCT_SOFT_TISSUE_PEAK  Identify the soft tissue peak from a UTE image intensity distribution
%
% Constrains a bias-corrected UTE image to grey and white matter voxels,
% finds the peak of the log-intensity distribution, and saves the peak
% value and histogram figure to path_pct.
%
% Use as:
%   pct_soft_tissue_peak(simnibs_folder, path_pct)
%
% Input:
%   simnibs_folder - path to the subject's m2m folder (contains final_tissues.nii.gz)
%   path_pct       - pCT processing folder (contains UTE_reg_thr0_corr.nii.gz)
%
% See also: PCT_SKULLMAPPING, FIT_PAIRWISELINEAR

arguments
    simnibs_folder (1,:) char
    path_pct       (1,:) char
end

    % Define paths to subject-specific files
    ute_file = fullfile(path_pct, 'UTE_reg_thr0_corr.nii.gz');
    ute_corr = niftiread(ute_file);

    % Constrain analysis to grey matter (label 1) and white matter (label 2)
    mask_file = fullfile(simnibs_folder, 'final_tissues.nii.gz');
    mask_corr = niftiread(mask_file);
    mask_corr = logical(mask_corr == 1 | mask_corr == 2); % Hard-coded SimNIBS labels
    ute_corr = ute_corr(mask_corr); % Apply mask to UTE image

    % Visualize histograms of UTE intensity distribution
    h = figure;
    tiledlayout(2, 1);

    % Histogram of raw UTE intensities
    nexttile;
    histogram(ute_corr, 'BinLimits', [0, 15000]);
    title('Intensity distribution UTE (GM/WM)');
    xlabel('UTE');

    % Histogram of logarithmic UTE intensities
    nexttile;
    histogram(-log(ute_corr), 'BinLimits', [-15, 0]);
    xlabel('-log(UTE)');

    % Save histogram visualization as an image
    figureName = 'pCT_histogram';
    saveas(h, fullfile(path_pct, figureName), 'png');

    % Extract soft tissue peak from logarithmic intensity distribution
    log_ute_corr = log(ute_corr);
    [counts, bins] = histcounts(log_ute_corr); % Compute histogram counts and bin edges
    [~, max_idx] = max(counts); % Find bin with maximum count
    peak_value = bins(max_idx); % Peak bin in logarithmic space

    % Convert peak value back to raw intensity scale
    peak = exp(peak_value);
    
    % Display soft tissue peak value
    fprintf('The soft tissue peak is at the intensity %.2f \n', peak);

    % Save soft tissue peak value to a text file in the pseudoCT folder
    txt_file = fopen(fullfile(path_pct, 'pCT_soft_tissue_value.txt'), 'w');
    fprintf(txt_file, '%.2f\n', peak);
    fclose(txt_file);

end