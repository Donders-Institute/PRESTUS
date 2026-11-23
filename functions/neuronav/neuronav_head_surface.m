function [F, V] = neuronav_head_surface(file)
% NEURONAV_HEAD_SURFACE  Outer head surface of an MNI-space NIfTI, for plots
%
% A label image (few integer values, e.g. final_tissues_MNI) is taken as
% label > 0; an intensity template (e.g. MNI152_T1) is thresholded at a
% fraction of its robust maximum. Only the largest component is kept, with
% its cavities filled, so the result is the outer shell. Sampled at ~2 mm
% and reduced to at most ~20000 faces. Used by the MNI views of
% NEURONAV_GUI and the PRESTUS Nexus.
%
% Use as:
%   [F, V] = neuronav_head_surface(file)
%
% Input:
%   file - NIfTI in MNI space
%
% Output:
%   F - [N_f x 3] faces (1-based vertex indices)
%   V - [N_v x 3] vertices, world mm
%
% See also: NEURONAV_FIND_MNI_TEMPLATE, NEURONAV_SUMMARY_TO_MNI

arguments
    file {mustBeTextScalar}
end

    info = niftiinfo(file);
    vol = double(niftiread(info));
    vol = vol(:, :, :, 1);
    step = max(1, round(2 / mean(info.PixelDimensions(1:3))));
    vol = vol(1:step:end, 1:step:end, 1:step:end);
    vals = unique(vol(1:53:end));
    if numel(vals) < 64 && all(vals == round(vals))
        mask = vol > 0;
    else
        vol = smooth3(vol, 'box', 3);
        mask = vol > 0.2 * prctile(vol(vol > 0), 98);
    end
    cc = bwconncomp(mask, 26);
    if cc.NumObjects == 0, error('no head voxels in %s', file); end
    [~, big] = max(cellfun(@numel, cc.PixelIdxList));
    mask = false(size(mask)); mask(cc.PixelIdxList{big}) = true;
    mask = imfill(mask, 'holes');
    for z = 1:size(mask, 3)   % cavities open at the cut-off neck
        mask(:, :, z) = imfill(mask(:, :, z), 'holes');
    end
    pad = zeros(size(mask) + 2);   % closes the surface at the volume edge
    pad(2:end-1, 2:end-1, 2:end-1) = mask;
    fv = isosurface(smooth3(pad, 'box', 3), 0.5);
    fv = reducepatch(fv, min(1, 20000 / size(fv.faces, 1)));
    % isosurface returns [column row slice] of the padded volume; the NIfTI
    % affine takes 0-based [i j k] (as in NEURONAV_APPLY_DEFORMATION).
    ijk = (fv.vertices(:, [2 1 3]) - 2) * step;
    V = [ijk, ones(size(ijk, 1), 1)] * info.Transform.T;
    V = V(:, 1:3);
    F = fv.faces;
end
