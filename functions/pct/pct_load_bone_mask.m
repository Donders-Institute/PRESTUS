function [bone, bone_file] = pct_load_bone_mask(parameters, tissues_mask_image, tissues_mask_header)
% PCT_LOAD_BONE_MASK  The bone mask the skull is built from, on the segmentation grid
%
% pct.bone_source selects the bone mask that PREPROC_HEAD and
% HEAD_MODEL_LOAD_OR_BUILD carry onto the simulation grid:
%   - 'segmentation' and 'hu': the segmentation's skull labels
%     (charm_seg_labels().bonemask); with 'hu', SKULL_FROM_HU replaces it on
%     the simulation grid
%   - 'mask': a bone mask from a file (path.bone_mask_pattern), for instance
%     the bone of another planning tool. Nonzero voxels are bone. It is
%     resampled (nearest neighbour) onto the segmentation grid through both
%     images' world (RAS mm) affines, so it must be co-registered to the
%     subject's T1w; voxels outside it are not bone. SKULL_RECONCILE makes it
%     the skull layer on the simulation grid.
%
% Use as:
%   [bone, bone_file] = pct_load_bone_mask(parameters, tissues_mask_image, tissues_mask_header)
%
% Input:
%   parameters          - PRESTUS parameters (pct.bone_source, path.bone_mask_pattern,
%                         path.anat, subject_id)
%   tissues_mask_image  - segmentation labels (as given to PREPROC_HEAD)
%   tissues_mask_header - its NIfTI header (after ENSURE_RAS_PLUS)
%
% Output:
%   bone      - logical bone mask on the segmentation grid
%   bone_file - the file read ('' for the segmentation's skull)
%
% See also: PCT_BONE_MASK_FILE, SKULL_RECONCILE, PCT_LOAD_CT

arguments
    parameters          (1,1) struct
    tissues_mask_image
    tissues_mask_header (1,1) struct
end

    bone_file = '';
    if ~strcmp(char(get_field(parameters, {'pct', 'bone_source'}, 'segmentation')), 'mask')
        bone = ismember(tissues_mask_image, charm_seg_labels().bonemask);
        return
    end

    bone_file = pct_bone_mask_file(parameters);
    hdr = niftiinfo(bone_file);
    [img, hdr] = ensure_ras_plus(niftiread(bone_file), hdr);
    % 1-based mask voxel -> world -> 1-based segmentation voxel, row-vector
    % convention ([vox0 1] * T = [ras 1]) as AFFINE_RESAMPLE_3D expects
    to0 = eye(4); to0(4, 1:3) = -1;
    to1 = eye(4); to1(4, 1:3) =  1;
    M = to0 * hdr.Transform.T / tissues_mask_header.Transform.T * to1;
    bone = affine_resample_3d(double(img ~= 0), M, size(tissues_mask_image, 1:3), 'nearest', 0) > 0.5;

    seg = ismember(tissues_mask_image, charm_seg_labels().bonemask);
    fprintf('[pct_load_bone_mask] %s: %d bone voxels on the segmentation grid; %d in the segmentation skull (%d)\n', ...
        bone_file, nnz(bone), nnz(bone & seg), nnz(seg));
    if ~any(bone(:))
        error('PRESTUS:pct:emptyBoneMask', 'The bone mask %s has no voxel on the segmentation grid.', bone_file);
    end
end
