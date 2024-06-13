function convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters, options)
    % runs subject2mni from SimNIBS
    arguments
        path_to_input_img string
        m2m_folder string
        path_to_output_img string
        parameters struct
        options.interpolation_order = 1
    end
    
    if isfield(parameters,'ld_library_path')
        ld_command = sprintf('export LD_LIBRARY_PATH="%s"; ', parameters.ld_library_path);
    else
        ld_command = '';
    end
    system(sprintf('%s%s/subject2mni --in %s --out %s --m2mpath %s --interpolation_order %d;', ld_command, parameters.simnibs_bin_path, path_to_input_img, path_to_output_img, m2m_folder, options.interpolation_order))
    if ~matches(path_to_output_img,'_MNI.nii.gz') %SimNIBS adds an unneeded affix to the file name, removing it
        simnibs_name = strrep(path_to_output_img, '.nii.gz', '_MNI.nii.gz');
    end
    system(sprintf('mv %s %s', simnibs_name, path_to_output_img));

end