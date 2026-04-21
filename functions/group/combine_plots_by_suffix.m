function combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters, options)

% COMBINE_PLOTS_BY_SUFFIX  Combine per-subject plots into a labelled montage image
%
% Reads individual subject plot files matching the given suffix, adds
% subject ID text labels, and assembles them into a single montage PNG
% saved in outputs_path as all_<suffix>.png.
%
% Use as:
%   combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters)
%   combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters, options)
%
% Input:
%   suffix       - filename suffix identifying the plots to combine
%   outputs_path - directory containing subject output files
%   subject_list - array of subject identifiers
%   parameters   - (1,1) simulation parameters struct
%   options      - name-value options: font_size (default: 64)
%
% See also: CREATE_GROUP_MNI_PLOTS

    arguments
        suffix       (1,:) char
        outputs_path (1,:) char
        subject_list
        parameters   (1,1) struct
        options.font_size (1,1) double = 64
    end

    % Ensure temporary directory exists for intermediate files
    if ~exist(sprintf('%s/tmp/', outputs_path), 'dir')
        mkdir(sprintf('%s/tmp/', outputs_path));
    end

    % Loop through each subject in the subject list
    for subject_i = 1:length(subject_list)
        % Determine the path to the subject's output directory
        if parameters.path.subject_subfolder
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

