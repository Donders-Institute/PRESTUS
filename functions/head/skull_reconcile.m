function [medium_masks, bone, pseudoCT, info] = skull_reconcile(parameters, medium_masks, bone, pseudoCT, source)
% SKULL_RECONCILE  Make a bone mask the skull layer, reconciled with the segmentation
%
% With pct.bone_source 'hu' or 'mask', the skull is a bone mask that is not
% the segmentation's skull. The two disagree at the skull surface, and the
% same steps apply whatever the mask came from (SKULL_FROM_HU: an HU
% threshold; PCT_LOAD_BONE_MASK: a file):
%   1. median filter of the bone mask, pct.bone_median_mm wide (0 = off)
%   2. opening with a cube of pct.bone_open_mm edge (0 = off): removes bone
%      specks and thin spurs narrower than the cube
%   3. closing with a cube of pct.bone_close_mm edge (0 = off)
%   4. keep the largest 26-connected component
%   5. with a (pseudo-)CT: raise bone HU below pct.bone_hu_threshold to it, so
%      no bone voxel is softer than bone (closing, interpolation, or a mask
%      drawn on another grid put soft-tissue HU into bone)
%   6. bone outside the segmentation skull becomes skull; segmentation skull
%      outside the bone joins pct.bone_soft_layer: 'nearest' (default) gives
%      each voxel the medium of the nearest voxel that is neither segmentation
%      skull nor bone, so the outer side joins the skin and the inner side the
%      brain or CSF; a layer name assigns that one layer to all of them.
%
% Use as:
%   [medium_masks, bone, pseudoCT, info] = skull_reconcile(parameters, medium_masks, bone, pseudoCT, source)
%
% Input:
%   parameters   - (1,1) configuration struct (pct.bone_*, grid.resolution_mm,
%                  layers, medium_properties)
%   medium_masks - [Nx x Ny x Nz] medium label array (index into
%                  fieldnames(parameters.medium_properties)), on the simulation grid
%   bone         - [Nx x Ny x Nz] logical bone mask on the same grid
%   pseudoCT     - HU on the same grid, or [] (no HU floor)
%   source       - text for the log line, e.g. 'HU > 300' or 'mask'
%
% Output:
%   medium_masks - updated medium label array
%   bone         - the bone mask after steps 1-4
%   pseudoCT     - HU with the bone floor applied ([] if given empty)
%   info         - struct: n_bone, n_seg_skull, n_to_soft, n_added, n_floor
%                  (voxels), and to_soft (voxels per medium that segmentation
%                  skull joined)
%
% See also: SKULL_FROM_HU, PCT_LOAD_BONE_MASK, HEAD_SMOOTH_AND_CROP

    arguments
        parameters   (1,1) struct
        medium_masks (:,:,:) {mustBeNumeric}
        bone         (:,:,:) {mustBeNumericOrLogical}
        pseudoCT     = []
        source       {mustBeTextScalar} = 'mask'
    end
    if ~isequal(size(bone), size(medium_masks))
        error('skull_reconcile:size', 'The bone mask and the medium masks must share the simulation grid.');
    end
    if ~isempty(pseudoCT) && ~isequal(size(pseudoCT), size(medium_masks))
        error('skull_reconcile:ctSize', 'The (pseudo-)CT and the medium masks must share the simulation grid.');
    end

    res  = parameters.grid.resolution_mm;
    thr  = get_field(parameters, {'pct', 'bone_hu_threshold'}, 300);
    med  = get_field(parameters, {'pct', 'bone_median_mm'}, 0);
    opn  = get_field(parameters, {'pct', 'bone_open_mm'}, 0);
    clo  = get_field(parameters, {'pct', 'bone_close_mm'}, 0);
    soft = char(get_field(parameters, {'pct', 'bone_soft_layer'}, 'nearest'));

    labels = fieldnames(parameters.medium_properties);
    skull_i = find(strcmp(labels, 'skull'));
    if isempty(skull_i)
        error('skull_reconcile:noSkull', 'pct.bone_source ''%s'' needs a ''skull'' layer in medium_properties.', ...
            char(get_field(parameters, {'pct', 'bone_source'}, 'mask')));
    end
    soft_nearest = strcmp(soft, 'nearest');
    soft_i = find(strcmp(labels, soft));
    if ~soft_nearest && (isempty(soft_i) || ~isfield(parameters.layers, soft))
        error('skull_reconcile:badSoftLayer', ['pct.bone_soft_layer ''%s'' is neither ''nearest'' nor a ' ...
            'layer of this run.'], soft);
    end
    seg_skull = ismember(medium_masks, find(contains(labels, 'skull')));

    bone = logical(bone);
    if med > 0
        n = max(1, 2 * round((med / res - 1) / 2) + 1);          % odd width in voxels
        if n > 1, bone = medfilt3(bone, [n n n], 'zeros'); end
    end
    if opn > 0
        n = max(1, round(opn / res));
        if n > 1, bone = imopen(bone, true(n, n, n)); end
    end
    if clo > 0
        n = max(1, round(clo / res));
        if n > 1, bone = imclose(bone, true(n, n, n)); end
    end
    cc = bwconncomp(bone, 26);
    if cc.NumObjects > 1
        [~, k] = max(cellfun(@numel, cc.PixelIdxList));
        bone = false(size(bone));
        bone(cc.PixelIdxList{k}) = true;
    end

    n_floor = 0;
    if ~isempty(pseudoCT)
        low = bone & pseudoCT < thr;
        n_floor = nnz(low);
        pseudoCT(low) = thr;
    end

    to_soft = seg_skull & ~bone;
    info = struct('n_bone', nnz(bone), 'n_seg_skull', nnz(seg_skull), ...
        'n_to_soft', nnz(to_soft), 'n_added', nnz(bone & ~seg_skull), 'n_floor', n_floor);
    if soft_nearest
        % medium of the nearest voxel outside both skull definitions
        keep = ~(seg_skull | bone);
        if any(to_soft, 'all') && any(keep, 'all')
            [~, nearest] = bwdist(keep);
            medium_masks(to_soft) = medium_masks(nearest(to_soft));
        end
    else
        medium_masks(to_soft) = soft_i;
    end
    joined = medium_masks(to_soft);
    info.to_soft = struct();
    for k = unique(joined(joined > 0))'   % label 0: voxels outside every layer have no name
        info.to_soft.(labels{k}) = nnz(joined == k);
    end
    medium_masks(bone) = skull_i;

    floor_txt = 'no CT, no HU floor';
    if ~isempty(pseudoCT), floor_txt = sprintf('%d raised to %g HU', n_floor, thr); end
    fprintf(['[skull_reconcile] bone = %s (median %g mm, opening %g mm, closing %g mm): %d voxels (%s); ' ...
        'segmentation skull %d, of which %d -> %s; %d added\n'], char(source), med, opn, clo, ...
        info.n_bone, floor_txt, info.n_seg_skull, info.n_to_soft, soft_summary(info.to_soft, soft), info.n_added);
end

function txt = soft_summary(to_soft, soft)
% 'nearest (skin 120, brain 45)' or the layer name
    txt = soft;
    f = fieldnames(to_soft);
    if strcmp(soft, 'nearest') && ~isempty(f)
        parts = cellfun(@(n) sprintf('%s %d', n, to_soft.(n)), f, 'UniformOutput', false);
        txt = sprintf('nearest (%s)', strjoin(parts, ', '));
    end
end
