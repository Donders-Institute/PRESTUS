function [medium_masks, bone_mask, pseudoCT, info] = skull_from_hu(parameters, medium_masks, pseudoCT)
% SKULL_FROM_HU  Define the skull layer by a Hounsfield threshold on the simulation grid
%
% With pct.bone_source 'hu', the skull is where the (pseudo-)CT exceeds
% pct.bone_hu_threshold on the simulation grid (HU already within
% pct.limits.hu_*), instead of the segmentation's skull labels. This follows
% BabelBrain's CT path (BabelDatasetPreps.py, v0.8.8). SKULL_RECONCILE then
% cleans the mask (median filter, closing, largest component), raises bone HU
% below the threshold to it, and reconciles the mask with the segmentation
% skull (pct.bone_soft_layer), as for a bone mask from a file.
%
% Use as:
%   [medium_masks, bone_mask, pseudoCT, info] = skull_from_hu(parameters, medium_masks, pseudoCT)
%
% Input:
%   parameters   - (1,1) configuration struct (pct.*, grid.resolution_mm, layers,
%                  medium_properties)
%   medium_masks - [Nx x Ny x Nz] medium label array (index into
%                  fieldnames(parameters.medium_properties)), on the simulation grid
%   pseudoCT     - [Nx x Ny x Nz] HU on the same grid
%
% Output:
%   medium_masks - updated medium label array
%   bone_mask    - [Nx x Ny x Nz] logical HU-defined bone
%   pseudoCT     - HU with the bone floor applied
%   info         - struct from SKULL_RECONCILE
%
% See also: SKULL_RECONCILE, HEAD_SMOOTH_AND_CROP, SKULL_FILL_HOLES

    arguments
        parameters   (1,1) struct
        medium_masks (:,:,:) {mustBeNumeric}
        pseudoCT     (:,:,:) {mustBeNumeric}
    end
    if isempty(pseudoCT) || ~isequal(size(pseudoCT), size(medium_masks))
        error('skull_from_hu:noCT', ['pct.bone_source ''hu'' needs a (pseudo-)CT on the simulation grid ' ...
            '(pct.enabled: 1).']);
    end
    thr = get_field(parameters, {'pct', 'bone_hu_threshold'}, 300);
    [medium_masks, bone_mask, pseudoCT, info] = skull_reconcile(parameters, medium_masks, ...
        pseudoCT > thr, pseudoCT, sprintf('HU > %g', thr));
end
