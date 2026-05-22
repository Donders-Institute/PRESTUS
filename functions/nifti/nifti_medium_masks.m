function nifti_medium_masks(parameters, planimg, medium_masks, is_layered, m2m_folder)
% NIFTI_MEDIUM_MASKS  Write medium_masks pipeline output to NIfTI
%
% Back-transforms tissue label map to T1 space and writes to dir_nii_T1w.
% Optionally converts to MNI space (controlled by parameters.io.save_MNI).
%
% Use as:
%   nifti_medium_masks(parameters, planimg, medium_masks, is_layered, m2m_folder)
%
% Input:
%   parameters   - PRESTUS config struct
%   planimg      - struct with fields: t1_image_orig, inv_transf, t1_header
%   medium_masks - (:,:,:) uint8, simulation-grid tissue labels
%   is_layered   - logical, true when medium is 'layered'
%   m2m_folder   - path to SimNIBS m2m folder
%
% See also: NIFTI_MEDIUM, NIFTI_PCT, NIFTI_TO_T1W, NIFTI_TO_MNI

masks_file    = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_medium_masks', ...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
masks_file_gz = strcat(masks_file, '.nii.gz');
mni_masks_file = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_medium_masks.nii.gz', ...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

nifti_to_t1w(uint8(medium_masks), masks_file, parameters, planimg, ...
    'IsLayered', is_layered, 'Datatype', 'uint8', 'BitsPerPixel', 8, 'Resampler', 'nearest');
nifti_to_mni(masks_file_gz, mni_masks_file, parameters, is_layered, m2m_folder);
end
