function nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, m2m_folder, options)
% NIFTI_TO_MNI  Convert a T1-space NIfTI to MNI space via SimNIBS
%
% Guards are applied before calling convert_final_to_MNI_simnibs: skips when
% io.save_MNI is false, the file already exists (overwrite guard), or the
% medium is not layered.
%
% MNI output is controlled by parameters.io.save_MNI (default: true via
% should_save_output; global fallback: parameters.io.save_matrices).
%
% When FillValue is provided, voxels in the MNI output with value 0 (the
% SimNIBS background default for out-of-FOV regions) are replaced with that
% value after warping.
%
% Use as:
%   nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, m2m_folder)
%   nifti_to_mni(..., 'FillValue', water_value)
%
% Input:
%   orig_file_gz - path to source NIfTI in T1 space (with '.nii.gz' extension)
%   mni_file     - output path in MNI space (with '.nii.gz' extension)
%   parameters   - PRESTUS config struct
%   is_layered   - logical, true when medium is 'layered'
%   m2m_folder   - path to SimNIBS m2m folder
%
% Options:
%   FillValue    - scalar; replaces zero-valued voxels (SimNIBS background)
%                  in the MNI output with this value (default: [], no fill)
%
% See also: NIFTI_WRITE_VOLUME, CONVERT_FINAL_TO_MNI_SIMNIBS, SHOULD_SAVE_OUTPUT

arguments
    orig_file_gz  string
    mni_file      string
    parameters    struct
    is_layered    logical
    m2m_folder    string
    options.FillValue (1,1) double = NaN  % NaN sentinel = "no fill"
end

if ~is_layered || ~should_save_output(parameters.io, 'save_MNI') || ...
        ~confirm_overwriting(mni_file, parameters)
    return
end
convert_final_to_MNI_simnibs(orig_file_gz, m2m_folder, mni_file, parameters, 'interpolation_order', 0);

% Replace zero-fill voxels (SimNIBS background for out-of-T1-FOV regions)
% with the requested water property value so the MNI map is consistent with
% the T1w map, where nifti_to_t1w already applies the same fill.
if ~isnan(options.FillValue) && isfile(mni_file)
    hdr  = niftiinfo(mni_file);
    vol  = niftiread(hdr);
    vol(vol == 0) = cast(options.FillValue, class(vol));
    niftiwrite(vol, mni_file, hdr, 'Compressed', true);
end
end
