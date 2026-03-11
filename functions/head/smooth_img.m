function thresholded_img = smooth_img(unsmoothed_img, fwhm_mm, voxel_size_mm, threshold, method)
%SMOOTH_IMG - Smooth image with FWHM [mm]
%
% SYNOPSIS:
%   smoothed = smooth_img(binary_mask, 1.0, 0.5, 0.5, 'gaussian')
%
% INPUT:
%   unsmoothed_img  - Input image [Nx Ny Nz] or [Nx Ny]
%   fwhm_mm         - FWHM in mm (scalar → isotropic, or 1xN per dim)
%   voxel_size_mm   - Voxel spacing in mm (scalar or [dx dy (dz)])
%   threshold       - Post-smoothing binarization [0.1-0.9] | 0: continuous
%   method          - ['gaussian'|'box']
%
% OUTPUT:
%   thresholded_img - Smoothed (binary if threshold > 0)
%
% DEFAULT: fwhm_mm=1.0, voxel_size_mm=0.5
%
% HOW IT WORKS:
%   Convert FWHM(mm) → FWHM(voxels) = fwhm_mm ./ voxel_size_mm
%   GAUSSIAN: FWHM(voxels) → σ(voxels) = fwhm_voxels/2.35482 → imgaussfilt(3)
%   BOX:      Kernel size = round(FWHM(voxels)) per dim (no min size)

    arguments
        unsmoothed_img {mustBeNumericOrLogical}
        fwhm_mm (1,:) double {mustBeNonnegative} = 1.0
        voxel_size_mm (1,:) double {mustBePositive} = 0.5
        threshold (1,1) double {mustBeInRange(threshold, 0, 1)} = 0.5
        method string {mustBeMember(method, ["gaussian", "box"])} = "gaussian"
    end

    img = double(unsmoothed_img);
    ndims_img = ndims(img);

    % Normalize fwhm_mm and voxel_size_mm to length ndims_img
    if isscalar(fwhm_mm)
        fwhm_mm = repmat(fwhm_mm, 1, ndims_img);
    end
    if isscalar(voxel_size_mm)
        voxel_size_mm = repmat(voxel_size_mm, 1, ndims_img);
    end
    if numel(fwhm_mm) ~= ndims_img || numel(voxel_size_mm) ~= ndims_img
        error('fwhm_mm and voxel_size_mm must be scalar or match image dimensionality (%dD).', ndims_img);
    end

    % Convert FWHM(mm) → FWHM(voxels) per dimension
    fwhm_voxels = fwhm_mm ./ voxel_size_mm;

    % Convert FWHM(voxels) → σ(voxels) for Gaussian [web:9]
    sigma_voxels = fwhm_voxels / 2.35482;  % FWHM = 2.35482 × σ [web:9]

    if any(fwhm_voxels > 0)
        fprintf('Smoothing: FWHM=[%s] mm, voxel=[%s] mm → FWHM=[%s] vox → σ=[%s] vox\n', ...
            num2str(fwhm_mm), num2str(voxel_size_mm), num2str(fwhm_voxels), num2str(sigma_voxels));

        if ndims_img == 2
            sig = sigma_voxels(1:2);

            switch method
                case "gaussian"
                    % imgaussfilt expects sigma in pixels/voxels [web:1][web:2]
                    smoothed_img = imgaussfilt(img, sig);

                case "box"
                    % Kernel size = round(FWHM voxels), no minimum
                    ksz = round(fwhm_voxels(1:2));
                    if any(ksz == 0)
                        smoothed_img = img;  % No smoothing
                        return;
                    end
                    % Force odd kernel size
                    ksz = ksz + mod(ksz+1,2);
                    kernel = ones(ksz) / prod(ksz);
                    smoothed_img = imfilter(img, kernel, 'replicate');  % Works for ksz=1 [web:24]
            end

        elseif ndims_img == 3
            sig = sigma_voxels(1:3);

            switch method
                case "gaussian"
                    % imgaussfilt3 sigma also in voxels [web:4]
                    smoothed_img = imgaussfilt3(img, sig);

                case "box"
                    ksz = round(fwhm_voxels(1:3));
                    if any(ksz == 0)
                        smoothed_img = img;  % No smoothing
                        return;
                    end
                    ksz = ksz + mod(ksz+1,2);
                    kernel = ones(ksz) / prod(ksz);
                    smoothed_img = convn(img, kernel, 'same');  % Works for ksz=1 [web:21]
            end
        else
            error('Supports 2D/3D only (ndims=%d)', ndims_img);
        end
    else
        smoothed_img = img;  % No smoothing
    end

    % Post-process binarization (optional)
    if threshold > 0
        thresholded_img = smoothed_img > threshold;
    else
        thresholded_img = smoothed_img;
    end
end
