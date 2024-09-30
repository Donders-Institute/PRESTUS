function combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters, options)
    arguments
        suffix string 
        outputs_path string 
        subject_list 
        parameters struct

        options.font_size = 64
    end
    if ~exist(sprintf('%s/tmp/', outputs_path),'dir')
        mkdir(sprintf('%s/tmp/', outputs_path));
    end
    
    for subject_i = 1:length(subject_list)
        if parameters.subject_subfolder
            sub_path = fullfile(outputs_path, sprintf('sub-%03i', subject_list(subject_i)));
        else
            sub_path = fullfile(outputs_path);
        end
        old_name = fullfile(sub_path, sprintf('sub-%03i_%s.png', subject_list(subject_i), suffix));
        orig_imsize = size(imread(old_name));
        I = imread(old_name);
        I = padarray(I, max([0 0 0; orig_imsize-size(I)],[], 1), I(1), 'pre');
        I = insertText(I,[0 5],sprintf('P%02i',subject_list(subject_i)),'FontSize', options.font_size, 'BoxOpacity',0,'TextColor','white');
        new_name = sprintf('%s/tmp/%i.png', outputs_path, subject_list(subject_i));
        imwrite(I, new_name);
    end

    system(sprintf('cd "%s/tmp/"; montage -background black -gravity south  -mode concatenate $(ls -1 *.png | sort -g)  ../all_%s.png', convertCharsToStrings(outputs_path), suffix));

    % remove temporary directory
    system(sprintf('cd "%s/tmp/"; rm *.png', outputs_path));
    rmdir(sprintf('%s/tmp/', outputs_path));
end
