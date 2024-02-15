function [min_dims, max_dims, grid_size] = get_crop_dims(image, margin)
    
    % This ensures that the CSF mask is not expanded outside the
    % bounds of the image dimensions

    I = find(image);
    [x,y,z]=ind2sub(size(image), I);
    min_dims = [min(x), min(y), min(z)]-margin;
    max_dims = [max(x), max(y), max(z)]+margin;
    grid_size = max_dims - min_dims + 1;
end