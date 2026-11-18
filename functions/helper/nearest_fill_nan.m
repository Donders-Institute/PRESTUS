function vol = nearest_fill_nan(vol)
% NEAREST_FILL_NAN  Replace NaN voxels with the value of the nearest non-NaN voxel
%
% Approximates the bwdist-based nearest fill
%   [~, idx] = bwdist(~isnan(vol));
%   vol(isnan(vol)) = vol(idx(isnan(vol)));
% by iterative 26-connected dilation: each pass sets the NaN voxels next to
% filled voxels to the mean of their filled 26-neighbours, until no NaN
% voxels remain (or no progress is made). Where those neighbours differ, the
% mean differs from the nearest value. An all-NaN input returns zeros.
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
