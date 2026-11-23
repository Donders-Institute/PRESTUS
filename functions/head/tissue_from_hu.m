function [medium_masks, info] = tissue_from_hu(parameters, medium_masks, pseudoCT)
% TISSUE_FROM_HU  Soft tissue against water by a CT threshold on the simulation grid
%
% With pct.soft_hu_threshold set, the (pseudo-)CT decides where tissue ends
% and water begins, instead of the segmentation alone. Outside the skull:
%   - tissue voxels below the threshold become water (for example air spaces
%     or a segmentation that reaches beyond the head)
%   - water voxels (label 0 or the water layer) at or above the threshold take
%     the medium of the nearest soft-tissue voxel (neither water nor skull)
% Skull voxels are left to the skull definition (pct.bone_source). Runs after
% the skull is set, on the HU already within pct.limits.hu_*. Tissue that the
% CSF-based crop (headmodel.csf_expansion) turned into water is tissue again
% where its HU reach the threshold.
%
% Use as:
%   [medium_masks, info] = tissue_from_hu(parameters, medium_masks, pseudoCT)
%
% Input:
%   parameters   - (1,1) configuration struct (pct.soft_hu_threshold, medium_properties)
%   medium_masks - [Nx x Ny x Nz] medium label array (index into
%                  fieldnames(parameters.medium_properties)), on the simulation grid
%   pseudoCT     - HU on the same grid
%
% Output:
%   medium_masks - updated medium label array (unchanged when the threshold is [])
%   info         - struct: threshold, n_to_water, n_to_tissue (voxels)
%
% See also: SKULL_FROM_HU, SKULL_RECONCILE, HEAD_SMOOTH_AND_CROP

    arguments
        parameters   (1,1) struct
        medium_masks (:,:,:) {mustBeNumeric}
        pseudoCT     = []
    end
    thr = get_field(parameters, {'pct', 'soft_hu_threshold'}, []);
    info = struct('threshold', thr, 'n_to_water', 0, 'n_to_tissue', 0);
    if isempty(thr)
        return
    end
    if isempty(pseudoCT) || ~isequal(size(pseudoCT), size(medium_masks))
        error('tissue_from_hu:noCT', ['pct.soft_hu_threshold needs a (pseudo-)CT on the simulation grid ' ...
            '(pct.enabled: 1).']);
    end

    labels  = fieldnames(parameters.medium_properties);
    water_i = find(strcmp(labels, 'water'));
    if isempty(water_i), water_i = 0; end
    skull    = ismember(medium_masks, find(contains(labels, 'skull')));
    is_water = medium_masks == 0 | medium_masks == water_i;
    soft     = ~skull & ~is_water;

    to_water  = soft & pseudoCT < thr;
    to_tissue = is_water & pseudoCT >= thr;
    if any(to_tissue, 'all') && any(soft, 'all')
        [~, nearest] = bwdist(soft);                 % nearest soft-tissue voxel, before any change
        medium_masks(to_tissue) = medium_masks(nearest(to_tissue));
    else
        to_tissue(:) = false;
    end
    medium_masks(to_water) = water_i;

    info.n_to_water  = nnz(to_water);
    info.n_to_tissue = nnz(to_tissue);
    fprintf('[tissue_from_hu] soft tissue = HU >= %g outside the skull: %d voxels -> water, %d water voxels -> tissue\n', ...
        thr, info.n_to_water, info.n_to_tissue);
end
