function thresholded_img = smooth_img(unsmoothed_img, fwhm_voxels, threshold, method)
%SMOOTH_IMG - Smooth tissue masks with voxel-space FWHM
%
% SYNOPSIS:
%   smoothed = smooth_img(binary_mask, 2, 0.5, 'gaussian')
%
% INPUT:
%   unsmoothed_img  - Input image [Nx Ny Nz] or [Nx Ny]
%   fwhm_voxels     - FWHM in grid voxels
%   threshold       - Post-smoothing binarization [0.1-0.9] | 0: continuous
%   method          - ['gaussian'|'box']
%
% OUTPUT:
%   thresholded_img - Smoothed (binary if threshold > 0)
%
% DEFAULT: fwhm_voxels=2 (e.g., 1mm equivalent @ 0.5mm grid)
%
% HOW IT WORKS:
%   GAUSSIAN: FWHM → σ = fwhm/2.35482 → imgaussfilt3(img, σ)
%     Creates smooth tissue gradients (center-high, edge-low weights)
%   BOX:     Uniform average over round(fwhm_voxels)^3 cube
%     Equal weights across kernel.

    arguments
        unsmoothed_img {mustBeNumericOrLogical}
        fwhm_voxels (1,1) double {mustBePositive} = 2.0
        threshold (1,1) double {mustBeInRange(threshold, 0, 1)} = 0.5
        method string {mustBeMember(method, ["gaussian", "box"])} = "gaussian"
    end
    
    img = double(unsmoothed_img);
    ndims_img = ndims(img);
    
    % Convert FWHM(voxels) → σ(voxels) for Gaussian
    sigma_voxels = fwhm_voxels / 2.35482;  % FWHM = 2.35482 × σ
    
    if fwhm_voxels > 0
        fprintf('Smoothing: FWHM=%.1f voxels → σ=%.2f voxels\n', fwhm_voxels, sigma_voxels);
        
        if ndims_img == 2
            switch method
                case "gaussian"
                    % === GAUSSIAN: Weighted average, bell-shaped kernel ===
                    % Weights: exp(-r²/(2σ²)). Center=1.0, edges~0.01
                    smoothed_img = imgaussfilt(img, sigma_voxels);
                    
                case "box"  
                    % === BOX: Uniform average over cubic kernel ===
                    % All voxels in kernel get equal weight (1/N_total)
                    kernel_size = max(3, round(fwhm_voxels));  % e.g. fwhm=4 → 4×4
                    kernel = ones(kernel_size) / kernel_size^2;  % Sum=1.0
                    smoothed_img = imfilter(img, kernel, 'replicate');
            end
            
        elseif ndims_img == 3
            switch method
                case "gaussian"
                    % === 3D GAUSSIAN: Rotational symmetric smoothing ===
                    % Each voxel = weighted average of 6σ-radius neighborhood
                    smoothed_img = imgaussfilt3(img, sigma_voxels);
                    
                case "box"
                    % === 3D BOX: Equal average over cubic volume ===
                    % e.g. fwhm=4 → [4 4 4] cube = 64 voxels, each weight=1/64
                    kernel_size = max(3, round(fwhm_voxels));
                    kernel = ones([kernel_size kernel_size kernel_size]) / kernel_size^3;
                    smoothed_img = convn(img, kernel, 'same');
            end
        else
            error('Supports 2D/3D only (ndims=%d)', ndims_img);
        end
    else
        smoothed_img = img;  % No smoothing
    end
    
    % Post-process binarization (optional)
    if threshold > 0
        thresholded_img = smoothed_img > threshold;  % Hard threshold
    else
        thresholded_img = smoothed_img;  % Keep continuous probabilities
    end
end
