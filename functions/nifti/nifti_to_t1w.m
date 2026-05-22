function nifti_to_t1w(data, orig_file, parameters, planimg, opts)
% NIFTI_TO_T1W  Back-transform one volume to T1 space and write NIfTI
%
% Confirms overwriting, applies affine_resample_3d to bring data from simulation-grid
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
%                              of the nearest in-grid voxel (nearest_fill_nan). Use when
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
        data_backtransf = single(affine_resample_3d(single(data), planimg.inv_transf, ...
            size(planimg.t1_image_orig), opts.Resampler, NaN));
        data_backtransf = nearest_fill_nan(data_backtransf);
        data_backtransf = cast(data_backtransf, opts.Datatype);
    else
        data_backtransf = cast(affine_resample_3d(single(data), planimg.inv_transf, ...
            size(planimg.t1_image_orig), opts.Resampler, double(opts.FillValue)), opts.Datatype);
    end

    niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
else
    res = parameters.grid.resolution_mm;

    % For axisymmetric phantom runs the acoustic/thermal outputs are radially
    % expanded to 3D [2*Nr x 2*Nr x Nz]. Extract the central xz-slice so the
    % saved NIfTI matches the 2D phantom input dimensions [Nr_full x Nz].
    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1 && ndims(data) == 3
        data = squeeze(data(round(size(data,1)/2), :, :));
    end

    % Write without header first so niftiinfo can parse the file dimensions,
    % then rewrite with the full spatial header.
    niftiwrite(cast(data, opts.Datatype), orig_file, 'Compressed', true);
    hdr = niftiinfo([orig_file '.nii.gz']);
    hdr.PixelDimensions = repmat(res, 1, numel(hdr.PixelDimensions));
    hdr.Datatype        = opts.Datatype;
    if ~isempty(opts.BitsPerPixel)
        hdr.BitsPerPixel = opts.BitsPerPixel;
    end

    % Reuse the phantom input NIfTI transform when available so output NIfTIs
    % share the same world-space orientation (including any z-flip) and can be
    % overlaid directly on the phantom in a viewer.
    if isfield(planimg, 'phantom_header') && ~isempty(planimg.phantom_header)
        hdr.Transform = planimg.phantom_header.Transform;
    else
        if isfield(planimg, 'origin_ras_mm') && ~isempty(planimg.origin_ras_mm)
            origin = planimg.origin_ras_mm(:);
        else
            origin = zeros(3, 1);
        end
        T = diag([res res res 1]);
        T(1:3, 4) = origin;
        hdr.Transform = affine3d(T');
    end

    niftiwrite(cast(data, opts.Datatype), orig_file, hdr, 'Compressed', true);
end
end
