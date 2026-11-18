function [out, info] = remap_acoustic_field(vol, M, out_dims, quantity)
% REMAP_ACOUSTIC_FIELD  Resample an acoustic volume onto another grid, by quantity
%
% Resamples vol (stored on a source grid or FOV box) onto a destination
% array of size out_dims under the forward voxel mapping M
% (GRID_REMAP_AFFINE; [vox_dst 1] = [vox_src 1] * M). How it is interpolated
% depends on what the volume holds:
%   'pressure'  peak or amplitude pressure >= 0: p^2 is interpolated
%               (trilinear) and the square root taken, so the intensity and
%               the heat source (both ~ p^2) keep their average, as the
%               separate thermal grid does (THERMAL_GRID_SETUP)
%   'negative'  peak negative pressure <= 0: as 'pressure' on |p|, sign kept
%   'complex'   complex amplitude: real and imaginary parts trilinear
%   'linear'    any other continuous quantity: trilinear
%   'mask'      logical or label volume: nearest neighbour
% Destination voxels outside the source volume are 0.
%
% Use as:
%   [out, info] = remap_acoustic_field(vol, M, out_dims, quantity)
%
% Input:
%   vol      - [Nx x Ny x Nz] source volume
%   M        - 4x4 forward affine, source voxel -> destination voxel
%   out_dims - [1x3] destination array size
%   quantity - 'pressure' | 'negative' | 'complex' | 'linear' | 'mask'
%
% Output:
%   out  - [out_dims] resampled volume (class of vol; logical for logical input)
%   info - .quantity, .coverage: fraction of destination voxels inside the
%          source volume
%
% See also: GRID_REMAP_AFFINE, GRID_WORLD_AFFINE, AFFINE_RESAMPLE_3D

    out_dims = double(out_dims(:)');
    cls = class(vol);
    vol = gather(vol);
    switch quantity
        case 'pressure'
            out = sqrt(max(affine_resample_3d(double(vol).^2, M, out_dims, 'linear', 0), 0));
        case 'negative'
            out = -sqrt(max(affine_resample_3d(double(vol).^2, M, out_dims, 'linear', 0), 0));
        case 'complex'
            out = complex(affine_resample_3d(double(real(vol)), M, out_dims, 'linear', 0), ...
                          affine_resample_3d(double(imag(vol)), M, out_dims, 'linear', 0));
        case 'linear'
            out = affine_resample_3d(double(vol), M, out_dims, 'linear', 0);
        case 'mask'
            out = affine_resample_3d(double(vol), M, out_dims, 'nearest', 0);
        otherwise
            error('remap_acoustic_field:quantity', ...
                'Unknown quantity ''%s'' (pressure, negative, complex, linear or mask).', quantity);
    end
    if strcmp(cls, 'logical')
        out = out > 0.5;
    else
        out = cast(out, cls);
    end
    inside = affine_resample_3d(ones(size(vol, 1), size(vol, 2), size(vol, 3)), M, out_dims, 'nearest', 0);
    info = struct('quantity', quantity, 'coverage', mean(inside(:) > 0));
end
