function thresholded_img = smooth_img(unsmoothed_img, fwhm_mm, voxel_size_mm, threshold, method)
% SMOOTH_IMG  Smooth an image with a Gaussian or box kernel defined by FWHM
%
% Converts FWHM in mm to voxel units and applies either a Gaussian
% or box-average filter. Optionally binarises the output
% with a user-specified threshold.
%
% Use as:
%   thresholded_img = smooth_img(unsmoothed_img, fwhm_mm, voxel_size_mm, threshold, method)
%
% Input:
%   unsmoothed_img - N-D numeric or logical image
%   fwhm_mm        - smoothing kernel FWHM [mm]; scalar or [1xNdim] (default: 1.0)
%   voxel_size_mm  - voxel spacing [mm]; scalar or [1xNdim] (default: 0.5)
%   threshold      - binarisation threshold 0–1; 0 returns continuous output (default: 0.5)
%   method         - 'gaussian', 'box', or 'off' to skip smoothing (default: 'gaussian')
%
% Output:
%   thresholded_img - smoothed image, binary when threshold > 0
%
% See also: SKULL_FILL_HOLES, HEAD_SMOOTH_AND_CROP

    arguments
        unsmoothed_img {mustBeNumericOrLogical}
        fwhm_mm (1,:) double {mustBeNonnegative} = 1.0
        voxel_size_mm (1,:) double {mustBePositive} = 0.5
        threshold (1,1) double {mustBeInRange(threshold, 0, 1)} = 0.5
        method string {mustBeMember(method, ["gaussian", "box", "off"])} = "gaussian"
    end

    img = double(unsmoothed_img);

    if method == "off"
        thresholded_img = img;
        return;
    end

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


        if ndims_img == 2
            sig = sigma_voxels(1:2);

            switch method
                case "gaussian"
                    hw = ceil(3 * max(sig));
                    [gx,gy] = ndgrid(-hw:hw, -hw:hw);
                    g = exp(-0.5*(gx.^2/sig(1)^2 + gy.^2/sig(2)^2));
                    g = g / sum(g(:));
                    pad = hw;
                    ip = padarray_replicate(img, [pad pad]);
                    smoothed_img = conv2(ip, g, 'valid');

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
                    pad = floor(ksz(1)/2);
                    ip = padarray_replicate(img, [pad pad]);
                    smoothed_img = conv2(ip, kernel, 'valid');
            end

        elseif ndims_img == 3
            sig = sigma_voxels(1:3);

            switch method
                case "gaussian"
                    hw = ceil(3 * max(sig));
                    x = (-hw:hw)';
                    g1 = exp(-0.5*(x.^2/sig(1)^2)); g1 = reshape(g1/sum(g1), [], 1, 1);
                    g2 = exp(-0.5*(x.^2/sig(2)^2)); g2 = reshape(g2/sum(g2), 1, [], 1);
                    g3 = exp(-0.5*(x.^2/sig(3)^2)); g3 = reshape(g3/sum(g3), 1, 1, []);
                    ip = padarray_replicate_3d(img, [hw hw hw]);
                    smoothed_img = convn(convn(convn(ip, g1, 'same'), g2, 'same'), g3, 'same');
                    smoothed_img = smoothed_img(hw+1:end-hw, hw+1:end-hw, hw+1:end-hw);

                case "box"
                    ksz = round(fwhm_voxels(1:3));
                    if any(ksz == 0)
                        smoothed_img = img;  % No smoothing
                        return;
                    end
                    ksz = ksz + mod(ksz+1,2);
                    kernel = ones(ksz) / prod(ksz);
                    pad = floor(ksz/2);
                    ip = padarray_replicate_3d(img, pad);
                    smoothed_img = convn(ip, kernel, 'same');
                    smoothed_img = smoothed_img(pad(1)+1:end-pad(1), pad(2)+1:end-pad(2), pad(3)+1:end-pad(3));
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

function out = padarray_replicate(img, pad)
% Replicate-pad a 2-D array by pad(1) rows and pad(2) cols
pr = pad(1); pc = pad(2);
out = [repmat(img(1,:),    pr, 1);
       img;
       repmat(img(end,:),  pr, 1)];
out = [repmat(out(:,1),   1, pc), out, repmat(out(:,end), 1, pc)];
end

function out = padarray_replicate_3d(img, pad)
% Replicate-pad a 3-D array by pad(1)/pad(2)/pad(3) voxels along dims 1/2/3
p1 = pad(1); p2 = pad(2); p3 = pad(3);
out = cat(1, repmat(img(1,:,:), p1, 1), img, repmat(img(end,:,:), p1, 1));
out = cat(2, repmat(out(:,1,:), 1, p2), out, repmat(out(:,end,:), 1, p2));
out = cat(3, repmat(out(:,:,1), 1, 1, p3), out, repmat(out(:,:,end), 1, 1, p3));
end
