function nifti_write_volume(data, orig_file, parameters, planimg, opts)
% NIFTI_WRITE_VOLUME  Back-transform one volume to T1 space and write NIfTI
%
% Confirms overwriting, applies tformarray to bring data from simulation-grid
% space to T1 space, and writes a compressed NIfTI. MNI conversion is the
% caller's responsibility (see NIFTI_TO_MNI).
%
% Use as:
%   nifti_write_volume(data, orig_file, parameters, planimg)
%   nifti_write_volume(..., 'Datatype', 'uint8', 'BitsPerPixel', 8, 'Resampler', 'nearest')
%
% Input:
%   data       - numeric array in simulation-grid space
%   orig_file  - output path without '.nii.gz' extension (T1 space)
%   parameters - PRESTUS config struct
%   planimg    - struct with fields: t1_image_orig, inv_transf, t1_header
%   opts       - name-value options (see arguments block)
%
% See also: NIFTI_TO_MNI, NIFTI_ACOUSTIC, NIFTI_MEDIUM, NIFTI_THERMAL

arguments
    data
    orig_file         (1,:) char
    parameters        (1,1) struct
    planimg           (1,1) struct
    opts.IsLayered    (1,1) logical = true
    opts.Datatype     (1,:) char    = 'single'
    opts.BitsPerPixel               = []
    opts.Resampler    (1,:) char    = 'cubic'
end

orig_file_gz = strcat(orig_file, '.nii.gz');

if ~confirm_overwriting(orig_file_gz, parameters)
    return
end

if opts.IsLayered
    orig_hdr = planimg.t1_header;
    orig_hdr.Datatype = opts.Datatype;
    if ~isempty(opts.BitsPerPixel)
        orig_hdr.BitsPerPixel = opts.BitsPerPixel;
    end
    data_backtransf = tformarray(data, planimg.inv_transf, ...
        makeresampler(opts.Resampler, 'fill'), [1 2 3], [1 2 3], ...
        size(planimg.t1_image_orig), [], 0);
    niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
else
    niftiwrite(data, orig_file, 'Compressed', true);
end
end
