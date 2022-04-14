function [min_dims, max_dims, grid_size] = get_crop_dims(image, margin)
I = find(image);
[x,y,z]=ind2sub(size(image), I);
min_dims = [min(x), min(y), min(z)]-margin;
max_dims = [max(x), max(y), max(z)]+margin;
grid_size = max_dims - min_dims;
% make grid divisible by n
n = 4;
grid_size = grid_size + (n - rem(grid_size,n));
max_dims = min_dims + grid_size-1; 
end