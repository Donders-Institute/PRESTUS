function combine_plots_by_suffix(suffix, outputs_path, subject_list, parameters, options)
    arguments
        suffix string 
        outputs_path string 
        subject_list 
        parameters struct

        options.font_size = 64
    end
    if ~exist(sprintf('%stmp/', outputs_path),'dir')
        mkdir(sprintf('%stmp/', outputs_path));
    end
    system(sprintf('cd "%stmp/"; rm *.png', outputs_path));
    if parameters.subject_subfolder
        orig_imsize = size(imread(fullfile(outputs_path, sprintf('sub-%03i/sub-%03i_%s.png', subject_list(1),subject_list(1), suffix))));
%        imarray = uint8(zeros(size(imread(fullfile(outputs_path, sprintf('sub-%03i/sub-%03i_%s.png', subject_list(end), subject_list(end), suffix))))));
    else
        orig_imsize = size(imread(fullfile(outputs_path, sprintf('sub-%03i_%s.png', subject_list(1),  suffix))));
%        imarray = uint8(zeros(size(imread(fullfile(outputs_path, sprintf('sub-%03i_%s.png', subject_list(end), suffix))))));
    end

%    imarray = repmat(imarray, 1, 1, 1, length(subject_list));
    for subject_i = 1:length(subject_list)
        if parameters.subject_subfolder
            old_name = fullfile(outputs_path, sprintf('sub-%03i/sub-%03i_%s.png', subject_list(subject_i), subject_list(subject_i), suffix));
        else
            old_name = fullfile(outputs_path, sprintf('sub-%03i_%s.png', subject_list(subject_i), suffix));
        end
        I = imread(old_name);
        I = padarray(I, max([0 0 0; orig_imsize-size(I)],[], 1), I(1), 'pre');

        I = insertText(I,[0 5],sprintf('P%02i',subject_list(subject_i)),'FontSize', options.font_size, 'BoxOpacity',0,'TextColor','white');
        %imarray(1:size(I, 1), 1:size(I, 2), 1:size(I, 3), subject_i) = I;
        new_name = sprintf('%stmp/%i.png', outputs_path, subject_list(subject_i));
        imwrite(I, new_name);
    end

    system(sprintf('cd "%stmp/"; montage -background black -gravity south  -mode concatenate $(ls -1 *.png | sort -g)  ../all_%s.png', convertCharsToStrings(outputs_path), suffix));
end
