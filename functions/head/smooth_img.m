function thresholded_img = smooth_img(unsmoothed_img, smooth_window, threshold, method)
%SMOOTH_IMG - Smooth tissue masks
%
% SYNOPSIS:
%   smoothed = smooth_img(binary_mask, smooth_window, threshold, 'anisotropic')
%
% INPUT:
%   unsmoothed_img - [Nx Ny Nz] logical/double binary mask
%   smooth_window  - Controls smoothing strength (iterations for AD)
%   threshold      - Post-smoothing binarization [0.1-0.9]
%   method         - Filter type ['gaussian'|'box']
%
% OUTPUT:
%   thresholded_img   - [Nx Ny Nz] logical smoothed binary mask
%
% DEFAULT BEHAVIOR:
%   method='gaussian'
%   iterations=smooth_window*2, conduction=0.125 (skull-optimized)
    arguments
        unsmoothed_img {mustBeNumericOrLogical}
        smooth_window (1,1) double {mustBePositive} = 4
        threshold (1,1) double {mustBeInRange(threshold, 0, 1)} = 0.5
        method string {mustBeMember(method, ["gaussian", "box"])} = "gaussian"
    end
    
    img = double(unsmoothed_img);
    ndims_img = ndims(img);
    
    if ndims_img == 2
        switch method
            case "gaussian"
                sigma = smooth_window / 3;
                filter_size = max(2, 2 * round(smooth_window / 2) + 1);  % Force odd: 5,7,9... with minimal size 2
                smoothed_img = imgaussfilt(img, sigma, 'FilterSize', filter_size);
                disp(['Gaussian filtering: sigma ', num2str(sigma) ,' size ', num2str(filter_size), ' voxels']);
            case "box"
                kernel = ones(smooth_window, smooth_window) / smooth_window^2;
                smoothed_img = imfilter(img, kernel, 'replicate');
        end
    elseif ndims_img == 3
        switch method
            case "gaussian"
                filter_size = max(2, 2 * round(smooth_window / 2) + 1);  % Also fix for consistency
                smoothed_img = smooth3(img, 'gaussian', filter_size);
                disp(['Gaussian filtering: size ', num2str(filter_size), ' voxels']);
            case "box"
                kernel = ones(smooth_window, smooth_window, smooth_window) / smooth_window^3;
                smoothed_img = convn(img, kernel, 'same');
        end
    else
        error('Supports 2D/3D only (ndims=%d)', ndims_img);
    end
    
    thresholded_img = smoothed_img > threshold;
end
