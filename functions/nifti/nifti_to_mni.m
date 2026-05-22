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
% Two optional fill modes address out-of-T1w-FOV voxels (those SimNIBS
% zero-fills because they fall outside the T1w image boundary):
%
%   'constant'  Replace out-of-FOV voxels with FillValue. When FovMask is
%               provided it identifies background precisely; otherwise voxels
%               equal to 0 are replaced. FillValue=NaN (default) means no
%               fill when FovMask is absent; when FovMask is provided, NaN
%               is written to background voxels.
%               Safe when 0 cannot be a valid voxel value, or when FovMask
%               is available (e.g. pCT, acoustic properties).
%
%   'nearest'   Replace out-of-FOV voxels (identified via FovMask) with the
%               value of the nearest in-FOV voxel. Requires FovMask.
%
% Use as:
%   nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, m2m_folder)
%   nifti_to_mni(..., 'FillMethod', 'constant', 'FillValue', val)
%   nifti_to_mni(..., 'FillMethod', 'constant', 'FillValue', NaN, 'FovMask', mask)
%   nifti_to_mni(..., 'FillMethod', 'nearest',  'FovMask',  mask)
%
% Input:
%   orig_file_gz - path to source NIfTI in T1 space (with '.nii.gz' extension)
%   mni_file     - output path in MNI space (with '.nii.gz' extension)
%   parameters   - PRESTUS config struct
%   is_layered   - logical, true when medium is 'layered'
%   m2m_folder   - path to SimNIBS m2m folder
%
% Options:
%   FillMethod   - 'constant' or 'nearest' (default: 'constant')
%   FillValue    - scalar used when FillMethod='constant' (default: NaN)
%   FovMask      - logical MNI-space array; 1=valid T1w source, 0=background
%                  Required for FillMethod='nearest'; optional for 'constant'.
%                  See MNI_FOV_MASK.
%
% See also: NIFTI_WRITE_VOLUME, CONVERT_FINAL_TO_MNI_SIMNIBS, MNI_FOV_MASK

arguments
    orig_file_gz   string
    mni_file       string
    parameters     struct
    is_layered     logical
    m2m_folder     string
    options.FillMethod  string        = 'constant'
    options.FillValue   (1,1) double  = NaN    % NaN = no fill
    options.FovMask                   = []
end

if ~is_layered || ~should_save_output(parameters.io, 'save_MNI') || ...
        ~confirm_overwriting(mni_file, parameters)
    return
end
mni_warp_method = 'simnibs';
if isfield(parameters, 'analysis') && isfield(parameters.analysis, 'mni_warp_method')
    mni_warp_method = parameters.analysis.mni_warp_method;
end
convert_final_to_MNI_simnibs(orig_file_gz, m2m_folder, mni_file, parameters, 'interpolation_order', 0, 'method', mni_warp_method);

if ~isfile(mni_file)
    return
end

switch options.FillMethod
    case 'constant'
        if isnan(options.FillValue) && isempty(options.FovMask)
            return
        end
        hdr = niftiinfo(mni_file);
        vol = niftiread(hdr);
        if ~isempty(options.FovMask)
            bg = ~options.FovMask;
        else
            bg = (vol == 0);
        end
        if isnan(options.FillValue)
            % NaN requires a float type; cast if necessary
            if ~isfloat(vol); vol = single(vol); hdr.Datatype = 'single'; hdr.BitsPerPixel = 32; end
            vol(bg) = single(NaN);
        else
            vol(bg) = cast(options.FillValue, class(vol));
        end
        niftiwrite(vol, strrep(mni_file, '.nii.gz', ''), hdr, 'Compressed', true);

    case 'nearest'
        if isempty(options.FovMask)
            return
        end
        hdr = niftiinfo(mni_file);
        vol = niftiread(hdr);
        bg  = ~options.FovMask;
        if any(bg, 'all')
            vol_tmp = vol;
            vol_tmp(bg) = NaN;
            vol_tmp = nearest_fill_nan(vol_tmp);
            vol(bg) = vol_tmp(bg);
            niftiwrite(vol, strrep(mni_file, '.nii.gz', ''), hdr, 'Compressed', true);
        end
end
end
