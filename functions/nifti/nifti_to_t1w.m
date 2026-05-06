function nifti_to_t1w(data, orig_file, parameters, planimg, opts)
% NIFTI_TO_T1W  Back-transform one volume to T1 space and write NIfTI
%
% Confirms overwriting, applies tformarray to bring data from simulation-grid
% space to T1 space, and writes a compressed NIfTI. MNI conversion is the
% caller's responsibility (see NIFTI_TO_MNI).
%
% Use as:
%   nifti_to_t1w(data, orig_file, parameters, planimg)
%   nifti_to_t1w(..., 'Datatype', 'uint8', 'BitsPerPixel', 8, 'Resampler', 'nearest')
%   nifti_to_t1w(..., 'FillMethod', 'constant', 'FillValue', water_val)
%   nifti_to_t1w(..., 'FillMethod', 'nearest')
%
% Input:
%   data       - numeric array in simulation-grid space
%   orig_file  - output path without '.nii.gz' extension (T1 space)
%   parameters - PRESTUS config struct
%   planimg    - struct with fields: t1_image_orig, inv_transf, t1_header
%   opts       - name-value options (see arguments block)
%
% Options:
%   FillMethod   - 'constant' (default) or 'nearest'.
%                  'constant': out-of-grid voxels get FillValue.
%                  'nearest':  out-of-grid voxels are filled with the value
%                              of the nearest in-grid voxel (bwdist). Use when
%                              FillValue is ambiguous (e.g. pCT where 0 HU is
%                              a valid tissue value).
%   FillValue    - scalar fill value for FillMethod='constant' (default: 0)
%
% See also: NIFTI_TO_MNI, NIFTI_ACOUSTIC, NIFTI_MEDIUM, NIFTI_THERMAL

arguments
    data
    orig_file           (1,:) char
    parameters          (1,1) struct
    planimg             (1,1) struct
    opts.IsLayered      (1,1) logical = true
    opts.Datatype       (1,:) char    = 'single'
    opts.BitsPerPixel               = []
    opts.Resampler      (1,:) char    = 'cubic'
    opts.FillMethod     (1,:) char    = 'constant'
    opts.FillValue      (1,1) double  = 0
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

    if strcmp(opts.FillMethod, 'nearest')
        % Fill with NaN as sentinel so out-of-grid voxels can be detected
        % after the transform, then replace each with its nearest in-grid value.
        data_backtransf = single(tformarray(single(data), planimg.inv_transf, ...
            makeresampler(opts.Resampler, 'fill'), [1 2 3], [1 2 3], ...
            size(planimg.t1_image_orig), [], NaN));
        bg = isnan(data_backtransf);
        if any(bg, 'all')
            [~, idx]             = bwdist(~bg);
            data_backtransf(bg)  = data_backtransf(idx(bg));
        end
        data_backtransf = cast(data_backtransf, opts.Datatype);
    else
        data_backtransf = tformarray(data, planimg.inv_transf, ...
            makeresampler(opts.Resampler, 'fill'), [1 2 3], [1 2 3], ...
            size(planimg.t1_image_orig), [], opts.FillValue);
    end

    niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
else
    if isfield(planimg, 'origin_ras_mm') && ~isempty(planimg.origin_ras_mm)
        % build a diagonal RAS affine from grid resolution and declared origin
        res = parameters.grid.resolution_mm;
        hdr = niftiinfo();
        hdr.ImageSize    = size(data);
        hdr.PixelDimensions = repmat(res, 1, ndims(data));
        hdr.Datatype     = opts.Datatype;
        T = diag([res res res 1]);
        T(1:3, 4) = planimg.origin_ras_mm(:);
        hdr.Transform    = affineTransform3d(T');
        niftiwrite(cast(data, opts.Datatype), orig_file, hdr, 'Compressed', true);
    else
        niftiwrite(data, orig_file, 'Compressed', true);
    end
end
end
