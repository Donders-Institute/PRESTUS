function nifti_pct(parameters, planimg, pseudoCT, m2m_folder)
% NIFTI_PCT  Write pseudo-CT pipeline output to NIfTI
%
% Back-transforms Hounsfield-unit skull map to T1 space and writes to
% dir_nii_T1w. Optionally converts to MNI space (controlled by parameters.io.save_MNI).
% Only called when parameters.pct.enabled == 1 (caller's responsibility).
%
% Use as:
%   nifti_pct(parameters, planimg, pseudoCT, m2m_folder)
%
% Input:
%   parameters - PRESTUS config struct
%   planimg    - struct with fields: t1_image_orig, inv_transf, t1_header
%   pseudoCT   - (:,:,:) single, simulation-grid Hounsfield values
%   m2m_folder - path to SimNIBS m2m folder
%
% See also: NIFTI_MEDIUM, NIFTI_MEDIUM_MASKS, NIFTI_TO_T1W, NIFTI_TO_MNI

pct_file     = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_pseudoCT', ...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
pct_file_gz  = strcat(pct_file, '.nii.gz');
pct_mni_file = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_pseudoCT.nii.gz', ...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

nifti_to_t1w(single(pseudoCT), pct_file, parameters, planimg, 'Resampler', 'nearest');
nifti_to_mni(pct_file_gz, pct_mni_file, parameters, true, m2m_folder);
end
