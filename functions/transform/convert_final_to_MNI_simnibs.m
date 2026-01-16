function convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters, options)

% CONVERT_FINAL_TO_MNI_SIMNIBS Converts an image to MNI space using SimNIBS.
%
% This function uses the SimNIBS `subject2mni` command to transform an image 
% from subject-specific space to MNI space. It also removes unnecessary affixes 
% added by SimNIBS to the output filename.
%
% Input:
%   path_to_input_img  - String specifying the path to the input image in subject space.
%   m2m_folder         - String specifying the path to the `m2m` folder containing transformation data.
%   path_to_output_img - String specifying the desired path for the output image in MNI space.
%   parameters         - Struct containing pipeline configuration parameters:
%                        * simnibs_bin_path: Path to SimNIBS binaries.
%                        * ld_library_path: Optional path for library linking (if required).
%
% Options:
%   interpolation_order - Integer specifying the interpolation order for resampling (default: 1).
%
% Output:
%   The transformed image is saved at `path_to_output_img`.

    arguments
        path_to_input_img string
        m2m_folder string
        path_to_output_img string
        parameters struct
        options.interpolation_order = 1 % Default interpolation order
    end

    % Check if LD_LIBRARY_PATH is specified and construct the export command if needed
    if isfield(parameters, 'ld_library_path')
        ld_command = sprintf('export LD_LIBRARY_PATH="%s"; ', parameters.ld_library_path);
    else
        ld_command = ''; % No library linking required
    end

    % Run SimNIBS `subject2mni` command to transform the image to MNI space
    system(sprintf('%s%s/subject2mni --in %s --out %s --m2mpath %s --interpolation_order %d;', ...
        ld_command, parameters.simnibs_bin_path, path_to_input_img, path_to_output_img, m2m_folder, options.interpolation_order));

    % Handle unnecessary affix added by SimNIBS to the output filename
    if ~matches(path_to_output_img, '_MNI.nii.gz') 
        simnibs_name = strrep(path_to_output_img, '.nii.gz', '_MNI.nii.gz'); % Generate SimNIBS default filename
    end

    % Rename the file to match the desired output filename
    system(sprintf('mv %s %s', simnibs_name, path_to_output_img));
end
