function convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters, options)
% CONVERT_FINAL_TO_MNI_SIMNIBS  Transform an image to MNI space using the SimNIBS subject2mni command
%
% Calls the SimNIBS subject2mni CLI tool to warp a NIfTI image from
% subject-specific (T1 conform) space to MNI space. An optional
% LD_LIBRARY_PATH export is prepended when parameters.hpc.ld_library_path
% is set. Because SimNIBS appends '_MNI' to the output filename by default,
% the result is renamed to path_to_output_img after the system call.
%
% Use as:
%   convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
%   convert_final_to_MNI_simnibs(..., 'interpolation_order', 0)
%
% Input:
%   path_to_input_img  - path to input NIfTI in subject space
%   m2m_folder         - path to SimNIBS m2m folder with registration data (toMNI subdirectory)
%   path_to_output_img - desired output path in MNI space (.nii.gz)
%   parameters         - PRESTUS config; uses startup.simnibs_bin_path and hpc.ld_library_path
%   options.interpolation_order - resampling order passed to subject2mni
%                                 (0=nearest, 1=linear; default: 1)
%
% See also: SIMULATION_NIFTI, CONVERT_FINAL_TO_MNI_MATLAB, SUBJECT2MNI_COORDS_LDFIX

    arguments
        path_to_input_img string
        m2m_folder string
        path_to_output_img string
        parameters struct
        options.interpolation_order = 1 % Default interpolation order
    end

    % Check if LD_LIBRARY_PATH is specified and construct the export command if needed
    % (Unix only; Windows does not use LD_LIBRARY_PATH and has no `export` builtin)
    if isfield(parameters.hpc, 'ld_library_path') && ~isempty(parameters.hpc.ld_library_path) && ~ispc
        ld_command = sprintf('export LD_LIBRARY_PATH="%s"; ', parameters.hpc.ld_library_path);
    else
        ld_command = ''; % No library linking required
    end

    % Run SimNIBS `subject2mni` command to transform the image to MNI space.
    % Build binary path with fullfile (cross-platform separator) and quote all
    % path arguments so spaces in paths don't break the command.
    subject2mni_bin = fullfile(parameters.startup.simnibs_bin_path, 'subject2mni');
    system(sprintf('%s"%s" --in "%s" --out "%s" --m2mpath "%s" --interpolation_order %d', ...
        ld_command, subject2mni_bin, path_to_input_img, path_to_output_img, m2m_folder, options.interpolation_order));

    % Handle unnecessary affix added by SimNIBS to the output filename
    if ~matches(path_to_output_img, '_MNI.nii.gz')
        simnibs_name = strrep(path_to_output_img, '.nii.gz', '_MNI.nii.gz'); % Generate SimNIBS default filename
    end

    % Rename the file to match the desired output filename.
    % Use MATLAB's movefile instead of `mv` so this works on Windows as well.
    movefile(simnibs_name, path_to_output_img);
end
