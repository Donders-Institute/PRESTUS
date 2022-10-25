function  [val,Ix,Iy,Iz] = masked_max_3d(array_3d, mask)
    
    % Takes the 3D intensity matrix of the unmasked region and
    % calculates the maximum intensity on the three separate axes.

    array_3d(~mask) = nan;
    [val,I] = max(array_3d(:));
    [Ix,Iy,Iz] = ind2sub(size(array_3d), I);
end