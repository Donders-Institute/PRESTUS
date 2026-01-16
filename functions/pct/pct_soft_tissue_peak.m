function pct_soft_tissue_peak(subject_id, base_path)

% PCT_SOFT_TISSUE_PEAK Identifies the soft tissue peak from a UTE image's intensity distribution.
%
% This function processes a bias field-corrected UTE image constrained to grey and white matter 
% to analyze its intensity distribution. It identifies the soft tissue peak by analyzing the 
% logarithmic distribution of intensity values and saves the result in a text file.
%
% Input:
%   subject_id - String specifying the subject ID (e.g., '001').
%   base_path  - String specifying the base path to the subject's m2m folder.
%
% Output:
%   None. The function saves histograms and the soft tissue peak value in the subject's folder.

    % Define paths to subject-specific files
    subject_folder = fullfile(base_path, sprintf('m2m_sub-%03s', subject_id));
    ute_file = fullfile(subject_folder, 'UTE_reg_thr0_corr.nii.gz');
    ute_corr = niftiread(ute_file);

    % Constrain analysis to grey matter (label 1) and white matter (label 2)
    mask_file = fullfile(subject_folder, 'final_tissues.nii.gz');
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
    saveas(h, fullfile(subject_folder, 'pseudoCT', figureName), 'png');

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
    txt_file = fopen(fullfile(subject_folder, 'pseudoCT', 'pCT_soft_tissue_value.txt'), 'w');
    fprintf(txt_file, '%.2f\n', peak);
    fclose(txt_file);

end