function vol = nearest_fill_nan(vol)
% NEAREST_FILL_NAN  Replace NaN voxels with the value of the nearest non-NaN voxel
%
% Drop-in replacement for the bwdist-based nearest-fill pattern:
%   [~, idx] = bwdist(~isnan(vol));
%   vol(isnan(vol)) = vol(idx(isnan(vol)));
%
% Uses iterative 26-connected dilation so that each NaN voxel is
% replaced by the average of its non-NaN 26-neighbours. Repeats until no
% NaN voxels remain (or no progress is made). The averaging gives
% identical results to true nearest-neighbour when the fill front advances
% one voxel at a time, and is a good approximation otherwise.
%
% Use as:
%   vol = nearest_fill_nan(vol)
%
% Input/Output:
%   vol - numeric array (any dimensionality); modified in-place
%
% See also: NIFTI_TO_T1W, NIFTI_TO_MNI

    bg = isnan(vol);
    if ~any(bg, 'all')
        return
    end

    vol(bg) = 0;
    K = ones(3, 3, 3, 'single');

    while true
        fg        = single(~bg);
        nb_sum    = convn(fg .* single(vol), K, 'same');
        nb_cnt    = convn(fg, K, 'same');
        new_fg    = bg & (nb_cnt > 0);
        if ~any(new_fg, 'all')
            break
        end
        vol(new_fg) = nb_sum(new_fg) ./ nb_cnt(new_fg);
        bg(new_fg)  = false;
    end
end
