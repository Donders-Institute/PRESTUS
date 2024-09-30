function smoothed_img = smooth_img(unsmoothed_img, windowSize, threshold)
% A function used to smooth and then binarize an image 
    arguments
       unsmoothed_img (:,:,:) double
       windowSize (1,1) double = 4 % size of the smoothing window
       threshold (1,1) double = 0.5 % binary threshold
    end
    kernel = ones(windowSize,windowSize,windowSize)/ windowSize ^ 3;
    blurryImage = convn(double(unsmoothed_img), kernel, 'same');
    smoothed_img = blurryImage > threshold; % Rethreshold
end
