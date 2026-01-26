function thresholded_img = smooth_img(unsmoothed_img, smooth_window, threshold, method)
%SMOOTH_IMG - Advanced edge-preserving smoothing for TUS tissue masks
% Improved: Anisotropic diffusion (default) + Gaussian fallback
%
% SYNOPSIS:
%   smoothed = smooth_img(binary_mask, smooth_window, threshold, 'anisotropic')
%
% INPUT:
%   unsmoothed_img - [Nx Ny Nz] logical/double binary mask
%   smooth_window  - Controls smoothing strength (iterations for AD)
%   threshold      - Post-smoothing binarization [0.1-0.9]
%   method         - Filter type ['anisotropic'|'gaussian'|'box']
%
% OUTPUT:
%   thresholded_img   - [Nx Ny Nz] logical smoothed binary mask
%
% DEFAULT BEHAVIOR:
%   method='anisotropic': imdiffusefilt, edge-preserving
%   iterations=smooth_window*2, conduction=0.125 (skull-optimized)

    arguments
        unsmoothed_img (:,:,:) {mustBeNumericOrLogical}
        smooth_window (1,1) double {mustBePositive} = 4
        threshold (1,1) double {mustBeInRange(threshold, 0, 1)} = 0.5
        method string {mustBeMember(method, ["anisotropic", "gaussian", "box"])} = "anisotropic"
    end
    
    img = double(unsmoothed_img);
    
    switch method
        case "anisotropic"
            iterations = round(smooth_window * 2.5);  % 10 for ws=4
            connectivity = 'minimal';
            smoothed_img = imdiffusefilt(img, ...
                'Connectivity', connectivity, 'NumberOfIterations', iterations);

        case "gaussian"
            size_vec = max(3, round(smooth_window / 2) * 2 + 1);  % Odd scalar
            smoothed_img = smooth3(img, 'gaussian', size_vec);
            
        case "box"
            % PRESTUS v0.4 implementation
            kernel = ones(smooth_window, smooth_window, smooth_window) / smooth_window^3;
            smoothed_img = convn(img, kernel, 'same');
    end
    
    thresholded_img = smoothed_img > threshold;
end
