function mni_ras = neuronav_apply_deformation(ras_seg, warp_field, warp_affine, nearest_idx)
% NEURONAV_APPLY_DEFORMATION  Map segmentation-space RAS coordinates through a SimNIBS coordinate field
%
% Samples a SimNIBS Conform2MNI_nonl/MNI2Conform_nonl field directly
% (trilinear interpolation) instead of shelling out to
% subject2mni_coords/mni2subject_coords per point. The field stores
% absolute target-space RAS coordinates per voxel (not a displacement),
% and MATLAB's 1-based niftiinfo affine requires a +1 voxel offset to
% agree with SimNIBS's own (0-based) convention -- both were verified
% against subject2mni_coords -t nonl to sub-millimetre agreement.
%
% SimNIBS zero-pads the field outside its valid domain (the same
% all-zero convention the CLI tools use to report "no mapping"). Points
% landing in that region are redirected to the nearest voxel with a
% defined mapping via nearest_idx (see NEURONAV_BUILD_DEFORMATION_CACHE)
% and read directly (no interpolation, since neighbouring voxels may
% themselves be undefined there).
%
% Use as:
%   mni_ras = neuronav_apply_deformation(ras_seg, warp_field, warp_affine, nearest_idx)
%
% Input:
%   ras_seg     - [Nx3] RAS mm coordinates in the field's source space
%   warp_field  - [X Y Z 3] SimNIBS coordinate field (as read by niftiread)
%   warp_affine - [4x4] voxel->RAS affine for warp_field (niftiinfo Transform.T)
%   nearest_idx - [X Y Z] nearest-valid-voxel linear index, from
%                 NEURONAV_BUILD_DEFORMATION_CACHE(warp_field)
%
% Output:
%   mni_ras - [Nx3] RAS mm coordinates in the field's target space
%
% See also: NEURONAV_BUILD_DEFORMATION_CACHE, NEURONAV_CONVERT_NATIVE_TO_MNI

    N = size(ras_seg, 1);
    dims = size(warp_field, 1:3);
    mni_ras = nan(N, 3);

    Fx = griddedInterpolant({1:dims(1), 1:dims(2), 1:dims(3)}, warp_field(:,:,:,1), 'linear', 'nearest');
    Fy = griddedInterpolant({1:dims(1), 1:dims(2), 1:dims(3)}, warp_field(:,:,:,2), 'linear', 'nearest');
    Fz = griddedInterpolant({1:dims(1), 1:dims(2), 1:dims(3)}, warp_field(:,:,:,3), 'linear', 'nearest');

    for i = 1:N
        vox = (warp_affine' \ [ras_seg(i,:)'; 1])';
        vox = vox(1:3) + 1;                                % 0-based -> 1-based (MATLAB) offset
        vox_c = max(min(vox, dims), [1 1 1]);              % keep in-bounds for indexing
        vox_r = round(vox_c);

        is_valid = any(squeeze(warp_field(vox_r(1), vox_r(2), vox_r(3), :)) ~= 0);

        % Trilinear interpolation uses the 8 surrounding voxels. Next to the edge
        % of the field some have no mapping (all zero) and would pull the result
        % towards the origin: interpolate over the voxels with a mapping only.
        lo8 = min(max(floor(vox_c), 1), max(dims - 1, 1));
        fr  = vox_c - lo8;
        acc = zeros(1, 3); wsum = 0;
        for cx = 0:1
            for cy = 0:1
                for cz = 0:1
                    c  = [cx cy cz];
                    ix = min(lo8 + c, dims);
                    v  = double(squeeze(warp_field(ix(1), ix(2), ix(3), :))');
                    if any(v ~= 0)
                        w = prod(c .* fr + (1 - c) .* (1 - fr));
                        acc = acc + w * v; wsum = wsum + w;
                    end
                end
            end
        end

        if wsum >= 1 - 1e-9
            mni_ras(i,:) = [Fx(vox_c(1), vox_c(2), vox_c(3)), ...
                             Fy(vox_c(1), vox_c(2), vox_c(3)), ...
                             Fz(vox_c(1), vox_c(2), vox_c(3))];
        elseif is_valid && wsum > 0
            mni_ras(i,:) = acc / wsum;
        else
            lin = sub2ind(dims, vox_r(1), vox_r(2), vox_r(3));
            [vx, vy, vz] = ind2sub(dims, nearest_idx(lin));
            mni_ras(i,:) = squeeze(warp_field(vx, vy, vz, :))';
        end
    end
end
