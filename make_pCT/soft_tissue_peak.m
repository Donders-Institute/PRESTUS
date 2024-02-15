% developed by Eleonora Carpino
% eleonora.carpino@icloud.com
% 14 June 2023

% Code to make an histogram of the intensities distribution from a bias
% field corrected UTE image and find the soft tissue peak

% Insert the subject numbers in the subject list
subject_list = [43,45,51];

for subject_id = subject_list
    
    % Modify the subject_folder path so that it is the segmentation output folder (m2m
    % folder in SimNIBS)
    subject_folder = fullfile('/home/action/elecar/MyProject/new_segmentations/',sprintf('m2m_sub-%1$03d/',subject_id));
    ute_corr = niftiread(fullfile(subject_folder,'UTE_reg_thr0_corr.nii.gz'));

    % To visually check the histogram
    tiledlayout(2,1)
    nexttile
    histogram(ute_corr,'BinLimits',[0,9000])
    title('Intensity distribution UTE');
    xlabel('UTE');
    nexttile
    histogram(-log(ute_corr),'BinLimits',[-9,0])
    xlabel('-log(UTE)');
    
    % Extract the value of the soft tissue peak from the histogram
    % distribution of the logarithm of the soft tissue intensity values
    log_ute_corr = log(ute_corr);
    %hist = histogram(log_ute_corr);
    [counts, bins] = histcounts(log_ute_corr);
    [~, max_idx] = max(counts);
    peak_value = bins(max_idx);
    
    % The value of the soft tissue peak is saved in a previously created
    % txt file called 'soft_tissue_value.txt'. The file is in the same m2m
    % folder
    peak = exp(peak_value);
    fprintf('The soft tissue peak is at the intensity %.2f \n', peak);
    txt_file = fopen(fullfile(subject_folder,'soft_tissue_value.txt'),'w');
    fprintf(txt_file,'%.2f\n',peak);
    fclose(txt_file);


end