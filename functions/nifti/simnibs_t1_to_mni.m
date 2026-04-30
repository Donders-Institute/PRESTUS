function simnibs_t1_to_mni(parameters, is_layered, m2m_folder)
% SIMNIBS_T1_TO_MNI  Register subject T1 to MNI space if not already done
%
% charm does not produce a T1-in-MNI file, so this step is run post-hoc once
% and existence-guarded so it is never repeated.
%
% MNI output is controlled by parameters.io.save_MNI (default: true via
% should_save_output; global fallback: parameters.io.save_matrices).
%
% Use as:
%   simnibs_t1_to_mni(parameters, is_layered, m2m_folder)
%
% Input:
%   parameters - PRESTUS config struct
%   is_layered - logical, true when medium is 'layered'
%   m2m_folder - path to SimNIBS m2m folder
%
% See also: NIFTI_MEDIUM, CONVERT_FINAL_TO_MNI_SIMNIBS, SHOULD_SAVE_OUTPUT

if ~is_layered || ~should_save_output(parameters.io, 'save_MNI')
    return
end

path_to_input_img  = fullfile(m2m_folder, 'T1.nii.gz');
path_to_output_img = fullfile(m2m_folder, 'toMNI', 'T1_to_MNI_post-hoc.nii.gz');

if ~exist(path_to_output_img, 'file')
    convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters);
end
end
