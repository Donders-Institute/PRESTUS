function nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, output_mni, m2m_folder)
% NIFTI_TO_MNI  Convert a T1-space NIfTI to MNI space via SimNIBS
%
% Guards are applied before calling convert_final_to_MNI_simnibs: skips when
% MNI output is disabled, the file already exists (overwrite guard), or the
% medium is not layered.
%
% Use as:
%   nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, output_mni, m2m_folder)
%
% Input:
%   orig_file_gz - path to source NIfTI in T1 space (with '.nii.gz' extension)
%   mni_file     - output path in MNI space (with '.nii.gz' extension)
%   parameters   - PRESTUS config struct
%   is_layered   - logical, true when medium is 'layered'
%   output_mni   - logical, whether MNI output is requested
%   m2m_folder   - path to SimNIBS m2m folder
%
% See also: NIFTI_WRITE_VOLUME, CONVERT_FINAL_TO_MNI_SIMNIBS

if ~output_mni || ~is_layered || ~confirm_overwriting(mni_file, parameters)
    return
end
convert_final_to_MNI_simnibs(orig_file_gz, m2m_folder, mni_file, parameters, 'interpolation_order', 0);
end
