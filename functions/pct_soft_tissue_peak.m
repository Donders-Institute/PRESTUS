
% Function to create an histogram of the intensities distribution from a bias
% field corrected UTE image and find the soft tissue peak

function pct_soft_tissue_peak(subject_id, base_path)
    
    subject_folder = fullfile(base_path, sprintf('m2m_sub-%03s', subject_id));
    ute_file = fullfile(subject_folder, 'UTE_reg_thr0_corr.nii.gz');
    ute_corr = niftiread(ute_file);

    % constrain to grey and white matter
    mask_file = fullfile(subject_folder, 'final_tissues.nii.gz');
    mask_corr = niftiread(mask_file);
    mask_corr = logical(mask_corr==1 | mask_corr==2); % Note: hard-coded to simnibs labels
    ute_corr = ute_corr(mask_corr);

    % Visually check the histogram
    h = figure;
    tiledlayout(2, 1)
    nexttile
    histogram(ute_corr, 'BinLimits', [0, 15000])
    title('Intensity distribution UTE (GM/WM)');
    xlabel('UTE');
    nexttile
    histogram(-log(ute_corr), 'BinLimits', [-15, 0])
    xlabel('-log(UTE)');
    figureName = ['pCT_histogram'];
    saveas(h, fullfile(subject_folder, 'pseudoCT', figureName), 'png');
    
    % Extract the value of the soft tissue peak from the histogram
    % distribution of the logarithm of the soft tissue intensity values
    log_ute_corr = log(ute_corr);
    [counts, bins] = histcounts(log_ute_corr);
    [~, max_idx] = max(counts);
    peak_value = bins(max_idx); % peak bin in log-ute
    
    % Saved soft tissue peak value in previously created txt file in m2m folder
    peak = exp(peak_value);
    fprintf('The soft tissue peak is at the intensity %.2f \n', peak);
    txt_file = fopen(fullfile(subject_folder, 'pseudoCT', 'pCT_soft_tissue_value.txt'), 'w');
    fprintf(txt_file, '%.2f\n', peak);
    fclose(txt_file);

    exit;
end