function hm = head_model_load_or_build(parameters, tissues_mask_image, tissues_mask_header, pseudoCT_image, pseudoCT_header)
% HEAD_MODEL_LOAD_OR_BUILD  Pose-independent head model of one subject, built once and cached
%
% The layered head model (medium masks smoothed and skull holes filled,
% segmentation, binary skull mask, pseudo-CT) does not depend on where the
% transducer is: only its rotation to the beam and the crop around head and
% transducer do (HEAD_MODEL_PLACE). With headmodel.reuse_across_placements,
% PREPROC_HEAD therefore builds the model once, unrotated (RAS+), rescaled to
% grid.resolution_mm, and every run of the subject (another pose, aim or
% affix) loads it instead of repeating the smoothing and hole filling, which
% dominate head preprocessing at fine resolution.
%
% The cache <dir_cache>/sub-XXX_<medium>_headmodel_<res>mm.mat stores the
% settings it depends on (CACHE_SIGNATURE kind 'headmodel': resolution,
% smoothing and hole-filling settings, layers, tissue labels, pseudo-CT
% settings, segmentation file date); a model built with other settings is
% rebuilt. It is written under a temporary name and renamed, so runs that
% build it at the same time do not leave a partial file.
%
% Use as:
%   hm = head_model_load_or_build(parameters, tissues_mask_image, tissues_mask_header, ...
%                                 pseudoCT_image, pseudoCT_header)
%
% Input:
%   parameters          - PRESTUS config (layered medium)
%   tissues_mask_image  - segmentation in the RAS+ T1 grid (PREPROC_HEAD)
%   tissues_mask_header - its NIfTI header
%   pseudoCT_image, pseudoCT_header - pseudo-CT in the same grid ([] without pCT)
%
% Output:
%   hm - struct: medium_masks, segmentation, bone_mask, pseudoCT ([] without
%        pCT), M (4x4: T1 voxel -> head-model voxel, row-vector convention as
%        the PREPROC_* transforms), resolution_mm, signature
%
% See also: HEAD_MODEL_PLACE, PREPROC_HEAD, PREPROC_RESCALE_ONLY, CACHE_SIGNATURE

arguments
    parameters          (1,1) struct
    tissues_mask_image
    tissues_mask_header (1,1) struct
    pseudoCT_image      = []
    pseudoCT_header     = [] %#ok<INUSA> header of the pseudo-CT (same grid as the segmentation)
end

    res = parameters.grid.resolution_mm;
    file = fullfile(parameters.io.dir_cache, sprintf('sub-%03d_%s_headmodel_%gmm.mat', ...
        parameters.subject_id, parameters.simulation.medium, res));
    sig = cache_signature(parameters, 'headmodel');

    if isfile(file)
        try
            stored = load(file, 'signature');
            if isequaln(stored.signature, sig)
                fprintf('[head_model] loading the head model: %s\n', file);
                hm = load(file);
                return
            end
            fprintf('[head_model] head model settings changed; rebuilding %s\n', file);
        catch ME
            fprintf('[head_model] could not read %s (%s); rebuilding\n', file, ME.message);
        end
    end

    log_timer('start', 'preproc_headmodel', parameters.io.dir_output);
    fprintf('[head_model] building the head model (unrotated, %.3g mm) ...\n', res);

    % Unrotated rescale to the grid resolution, built as PREPROC_RESCALE_ONLY
    % but from the image bounds alone, so the grid does not depend on any
    % transducer (labels nearest neighbour, pseudo-CT linear)
    sz = size(tissues_mask_image, 1:3);
    scale = double(tissues_mask_header.PixelDimensions(1)) / res;
    forward = makehgtform('translate', -(sz + 1) / 2);
    S = makehgtform('scale', scale);
    bounds = affine_bounds(forward' * S', [0 0 0; sz]);
    newdims = ceil(2 * max(abs(bounds), [], 1)) + 2;
    M = forward' * S' * makehgtform('translate', (newdims + 1) / 2)';
    segmentation = affine_resample_3d(tissues_mask_image, M, newdims, 'nearest', 0);
    bone_mask = affine_resample_3d(ismember(tissues_mask_image, charm_seg_labels().bonemask), M, newdims, 'nearest', 0);
    pseudoCT = [];
    if ~isempty(pseudoCT_image)
        if ~isequal(size(pseudoCT_image, 1:3), sz)
            error('PRESTUS:headModel:pctGrid', ...
                'The pseudo-CT and the segmentation must share a voxel grid for the shared head model.');
        end
        pseudoCT = affine_resample_3d(pseudoCT_image, M, newdims, 'linear', 0);
    end

    % Medium masks: smoothing and skull hole filling (as HEAD_SMOOTH_AND_CROP)
    log_timer('start', 'preproc_medium_mask', parameters.io.dir_output);
    medium_masks = preproc_medium_mask(segmentation, parameters);
    log_timer('stop', 'preproc_medium_mask');
    if any(contains(fieldnames(parameters.layers), 'skull'))
        medium_masks = skull_fill_holes(parameters, medium_masks, round(size(medium_masks) / 2), segmentation);
    end

    hm = struct('medium_masks', medium_masks, 'segmentation', segmentation, 'bone_mask', bone_mask, ...
        'pseudoCT', pseudoCT, 'M', M, 'resolution_mm', res, 'signature', sig);
    tmp = [file(1:end-4) '_' char(java.util.UUID.randomUUID) '.mat'];
    try
        save(tmp, '-struct', 'hm', '-v7.3');
        movefile(tmp, file, 'f');
        fprintf('[head_model] saved: %s\n', file);
    catch ME
        warning('PRESTUS:headModel:save', 'Could not save the head model (%s); it is rebuilt next time.', ME.message);
        if isfile(tmp), delete(tmp); end
    end
    log_timer('stop', 'preproc_headmodel');
end
