function nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, m2m_folder)
% NIFTI_TO_MNI  Convert a T1-space NIfTI to MNI space via SimNIBS
%
% Guards are applied before calling convert_final_to_MNI_simnibs: skips when
% io.save_MNI is false, the file already exists (overwrite guard), or the
% medium is not layered.
%
% MNI output is controlled by parameters.io.save_MNI (default: true via
% should_save_output; global fallback: parameters.io.save_matrices).
%
% Use as:
%   nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, m2m_folder)
%
% Input:
%   orig_file_gz - path to source NIfTI in T1 space (with '.nii.gz' extension)
%   mni_file     - output path in MNI space (with '.nii.gz' extension)
%   parameters   - PRESTUS config struct
%   is_layered   - logical, true when medium is 'layered'
%   m2m_folder   - path to SimNIBS m2m folder
%
% See also: NIFTI_WRITE_VOLUME, CONVERT_FINAL_TO_MNI_SIMNIBS, SHOULD_SAVE_OUTPUT

if ~is_layered || ~should_save_output(parameters.io, 'save_MNI') || ...
        ~confirm_overwriting(mni_file, parameters)
    return
end
convert_final_to_MNI_simnibs(orig_file_gz, m2m_folder, mni_file, parameters, 'interpolation_order', 0);
end
