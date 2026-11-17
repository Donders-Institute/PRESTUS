function simnibs_t1_to_mni(parameters, is_layered, m2m_folder)
% SIMNIBS_T1_TO_MNI  Register subject T1 to MNI space if not already done
%
% charm does not produce a T1-in-MNI file, so this step is run post-hoc once.
% The image is written next to the segmentation (GET_T1_MNI_FILE), not into
% the read-only m2m folder. Skipping is controlled by confirm_overwriting so
% overwrite_files = 'always' forces re-registration and 'never' skips an
% existing file.
%
% MNI output is controlled by parameters.io.save_MNI (default: true via
% should_save_output; not affected by parameters.io.save_matrices).
%
% Use as:
%   simnibs_t1_to_mni(parameters, is_layered, m2m_folder)
%
% Input:
%   parameters - PRESTUS config struct (path.seg, subject_id)
%   is_layered - logical, true when medium is 'layered'
%   m2m_folder - path to SimNIBS m2m folder (read only)
%
% See also: NIFTI_MEDIUM, CONVERT_FINAL_TO_MNI_SIMNIBS, SHOULD_SAVE_OUTPUT, GET_T1_MNI_FILE

if ~is_layered || ~should_save_output(parameters.io, 'save_MNI')
    return
end

path_to_input_img  = fullfile(m2m_folder, 'T1.nii.gz');
path_to_output_img = get_t1_mni_file(parameters.path.seg, parameters.subject_id);

if confirm_overwriting(path_to_output_img, parameters)
    out_dir = fileparts(path_to_output_img);
    if ~isfolder(out_dir); mkdir(out_dir); end
    mni_warp_method = 'simnibs';
    if isfield(parameters, 'analysis') && isfield(parameters.analysis, 'mni_warp_method')
        mni_warp_method = parameters.analysis.mni_warp_method;
    end
    convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters, 'method', mni_warp_method);
end
end
