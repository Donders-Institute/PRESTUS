% developed by Eleonora Carpino
% eleonora.carpino@icloud.com
% 14 June 2023

% Function to create an histogram of the intensities distribution from a bias
% field corrected UTE image and find the soft tissue peak

function find_soft_tissue_peak(subject_id, base_path)
    
    subject_folder = fullfile(base_path, sprintf('m2m_sub-%03s', subject_id));
    ute_file = fullfile(subject_folder, 'UTE_reg_thr0_corr.nii.gz');
    ute_corr = niftiread(ute_file);

    % Visually check the histogram
    figure;
    tiledlayout(2, 1)
    nexttile
    histogram(ute_corr, 'BinLimits', [0, 9000])
    title('Intensity distribution UTE');
    xlabel('UTE');
    nexttile
    histogram(-log(ute_corr), 'BinLimits', [-9, 0])
    xlabel('-log(UTE)');
    
    % Extract the value of the soft tissue peak from the histogram
    % distribution of the logarithm of the soft tissue intensity values
    log_ute_corr = log(ute_corr);
    [counts, bins] = histcounts(log_ute_corr);
    [~, max_idx] = max(counts);
    peak_value = bins(max_idx);
    
    % Saved soft tissue peak value in previously created txt file in m2m folder
    peak = exp(peak_value);
    fprintf('The soft tissue peak is at the intensity %.2f \n', peak);
    txt_file = fopen(fullfile(subject_folder, 'soft_tissue_value.txt'), 'w');
    fprintf(txt_file, '%.2f\n', peak);
    fclose(txt_file);
end