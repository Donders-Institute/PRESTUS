function combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters, options)

% COMBINE_PLOTS_BY_SUFFIX Combines subject-specific plots into a single montage image.
%
% This function reads individual subject plots (with a specified suffix), adds 
% text labels to each plot, and combines them into a single montage image. The 
% montage is created using the `montage` command-line tool. 
%
% Input:
%   suffix        - String specifying the suffix of the plot filenames to combine.
%   outputs_path  - String specifying the directory where output files are stored.
%   subject_list  - List of subject IDs for which plots will be combined.
%   parameters    - Struct containing additional configuration options (e.g., `subject_subfolder`).
%
% Options:
%   font_size     - Font size for text labels added to each plot (default: 64).
%
% Output:
%   A combined montage image is saved in the `outputs_path` directory with the 
%   filename `all_<suffix>.png`.

    arguments
        suffix string 
        outputs_path string 
        subject_list 
        parameters struct
        options.font_size = 64 % Default font size for text labels
    end

    % Ensure temporary directory exists for intermediate files
    if ~exist(sprintf('%s/tmp/', outputs_path), 'dir')
        mkdir(sprintf('%s/tmp/', outputs_path));
    end

    % Loop through each subject in the subject list
    for subject_i = 1:length(subject_list)
        % Determine the path to the subject's output directory
        if parameters.subject_subfolder
            sub_path = fullfile(outputs_path, sprintf('sub-%03i', subject_list(subject_i)));
        else
            sub_path = fullfile(outputs_path);
        end

        % Construct the filename of the subject's plot
        old_name = fullfile(sub_path, sprintf('sub-%03i_%s.png', subject_list(subject_i), suffix));

        % Read the image and ensure it matches the original image size by padding if necessary
        orig_imsize = size(imread(old_name)); % Get original image size
        I = imread(old_name); % Read image file
        I = padarray(I, max([0 0 0; orig_imsize - size(I)], [], 1), I(1), 'pre'); % Pad image if needed

        % Add a text label to the top-left corner of the image with subject ID
        I = insertText(I, [0 5], sprintf('P%02i', subject_list(subject_i)), ...
            'FontSize', options.font_size, 'BoxOpacity', 0, 'TextColor', 'white');

        % Save the modified image to the temporary directory
        new_name = sprintf('%s/tmp/%i.png', outputs_path, subject_list(subject_i));
        imwrite(I, new_name);
    end

    % Combine all modified images into a single montage using the `montage` command-line tool
    system(sprintf('cd "%s/tmp/"; montage -background black -gravity south -mode concatenate $(ls -1 *.png | sort -g) ../all_%s.png', ...
        convertCharsToStrings(outputs_path), suffix));

    % Clean up temporary files and remove the temporary directory
    system(sprintf('cd "%s/tmp/"; rm *.png', outputs_path));
    rmdir(sprintf('%s/tmp/', outputs_path));
end

