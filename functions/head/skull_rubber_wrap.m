function [SKULL_BALLON] = skull_rubber_wrap(parameters, BW, medium_masks, segmented_img)
% SKULL_RUBBER_WRAP  Fill skull holes where brain meets skin, as a rubber sheet around the skull would
%
% Closes the largest connected skull component with a sphere of radius
% headmodel.skull_wrap_radius_mm (default 5 mm): like a rubber sheet, the
% closing bridges holes narrower than about twice the radius. Closed voxels
% that touch both brain and skin (and are not brain) become skull, as do
% CSF voxels face-adjacent to skin within 2 voxels of them; single-voxel
% gaps along z are filled. With debug volumes on, the added voxels and the
% final mask are saved as NIfTI files.
%
% The closing is evaluated only where it can change the result: in blocks
% around the voxels next to both brain and skin, padded by twice the radius
% (the closing at a voxel depends on the skull within that distance). The
% cost therefore follows the size of the skull holes, not of the grid. The
% Image Processing Toolbox (bwdist) is used when its licence can be checked
% out; otherwise an exact base-MATLAB distance transform capped at the radius.
%
% Use as:
%   SKULL_BALLON = skull_rubber_wrap(parameters, BW, medium_masks, segmented_img)
%
% Input:
%   parameters    - (1,1) simulation configuration struct (grid.resolution_mm,
%                   headmodel.skull_wrap_radius_mm, layers, medium_properties)
%   BW            - [Nx x Ny x Nz] binary skull mask
%   medium_masks  - [Nx x Ny x Nz] medium label array
%   segmented_img - [Nx x Ny x Nz] SimNIBS tissue label volume
%
% Output:
%   SKULL_BALLON - [Nx x Ny x Nz] filled binary skull mask (logical)
%
% See also: SKULL_FILL_HOLES, PCT_SKULLEXPAND

    disp('Using rubber wrap to fill potential local holes in skull layer...')
    t_start = tic;

    sz = size(BW, 1:3);
    if ~isequal(size(medium_masks, 1:3), sz)
        error('skull_rubber_wrap:sizeMismatch', 'Size mismatch: skull %s vs tissues %s', ...
            mat2str(sz), mat2str(size(medium_masks)));
    end

    r_mm = 5;
    if isfield(parameters.headmodel, 'skull_wrap_radius_mm') && ~isempty(parameters.headmodel.skull_wrap_radius_mm)
        r_mm = parameters.headmodel.skull_wrap_radius_mm;
    end
    r = r_mm / parameters.grid.resolution_mm;   % [grid voxels]
    if ~(isscalar(r) && isfinite(r) && r > 0)
        error('skull_rubber_wrap:badRadius', ...
            'headmodel.skull_wrap_radius_mm must be a positive number (got %s).', mat2str(r_mm));
    end
    use_ipt = ipt_available();

    %% Largest connected skull component
    BW1 = largest_component(logical(BW), use_ipt);

    %% Voxels next to both brain and skin (face neighbours) that are not brain
    BRAIN = layer_mask(parameters, medium_masks, 'brain');
    SKIN  = layer_mask(parameters, medium_masks, 'skin');
    SKIN_touch = dilate6(SKIN);
    C = dilate6(BRAIN) & SKIN_touch & ~BRAIN;

    %% Rubber wrap: closing of the skull, evaluated at those voxels
    n_blocks = 0;
    B_out = false(sz);
    if any(C(:))
        [B_out, n_blocks] = close_at(BW1, C, r, use_ipt);
    end
    clear C

    %% CSF voxels face-adjacent to skin, within 2 voxels of the wrapped voxels
    if any(B_out(:))
        near = find(dilate_sparse(B_out, 2) & SKIN_touch);
        csf = near(ismember(segmented_img(near), charm_seg_labels().csf));
        B_out(csf) = true;
        B_out = B_out & ~BRAIN;
    end
    clear SKIN_touch

    SKULL = BW1;
    B0 = BW1 | B_out;
    clear B_out

    %% Single-voxel gaps along z with skull in front and behind
    ZADDED = false(sz);
    ZADDED(:, :, 2:end-1) = ~B0(:, :, 2:end-1) & B0(:, :, 1:end-2) & B0(:, :, 3:end);

    BALLOON = B0 & ~SKULL;          % voxels the wrap adds
    SKULL_BALLON = B0 | ZADDED;     % final filled skull mask (logical)

    if use_ipt, how = 'bwdist'; else, how = 'base-MATLAB distance transform'; end
    fprintf(['Rubber wrap: radius %g mm (%g voxels) on a [%s] grid; %d blocks with brain next to skin; ' ...
        '%d voxels added (%s, %.1f s)\n'], r_mm, r, num2str(sz), n_blocks, ...
        nnz(SKULL_BALLON & ~SKULL), how, toc(t_start));

    %% [DEBUG] Save as .nii.gz
    if debug_on(parameters, 'preproc', 'volumes')
        outNii = fullfile(parameters.io.dir_debug_preproc, ...
            sprintf('balloon_mask_added%s.nii', parameters.io.output_affix));
        niftiwrite(uint8(SKULL_BALLON & ~SKULL), outNii, 'Compressed', true);
        debug_register(parameters, 'preproc', outNii, ...
            'Rubber wrap: voxels added to the skull (grid voxels, no header)', 'volume');

        outNii = fullfile(parameters.io.dir_debug_preproc, ...
            sprintf('balloon_mask_final%s.nii', parameters.io.output_affix));
        niftiwrite(uint8(SKULL_BALLON), outNii, 'Compressed', true);
        debug_register(parameters, 'preproc', outNii, ...
            'Rubber wrap: final filled skull mask (grid voxels, no header)', 'volume');
    end

    %% [DEBUG] Visualize the expanded skull
    % Computationally heavy on the HPC; off by default (headmodel.skull_wrap_visualize)
    if isfield(parameters.headmodel, 'skull_wrap_visualize') && parameters.headmodel.skull_wrap_visualize == 1
        skull_rubber_wrap_visualize(parameters, SKULL, ZADDED, BALLOON);
    end

end

%% Helpers ---------------------------------------------------------------------

function tf = ipt_available()
% Image Processing Toolbox installed and its licence obtainable
    tf = license('test', 'image_toolbox') && ~isempty(which('bwdist'));
    if tf
        try
            tf = logical(license('checkout', 'image_toolbox'));
        catch
            tf = false;
        end
    end
end

function M = layer_mask(parameters, medium_masks, name)
% Voxels of a medium layer, or none when the layer is not modelled
    i = find(strcmp(fieldnames(parameters.medium_properties), name));
    if any(strcmp(fieldnames(parameters.layers), name)) && ~isempty(i)
        M = (medium_masks == i);
    else
        M = false(size(medium_masks, 1:3));
    end
end

function L = largest_component(BW, use_ipt)
% Largest 26-connected component, labelled within the bounding box of BW
    if ~any(BW(:))
        error('skull_rubber_wrap:emptySkull', 'Skull mask is empty.');
    end
    idx = bbox_idx(BW, 0);
    sub = BW(idx{:});
    if use_ipt
        CC = bwconncomp(sub, 26);
    else
        CC = bwconncomp_3d(sub, 26);
    end
    [~, iBig] = max(cellfun(@numel, CC.PixelIdxList));
    Lsub = false(size(sub));
    Lsub(CC.PixelIdxList{iBig}) = true;
    L = false(size(BW));
    L(idx{:}) = Lsub;
end

function idx = bbox_idx(M, pad)
% Index ranges of the bounding box of M, padded and clipped to the array
    sz = size(M, 1:3);
    idx = cell(1, 3);
    for d = 1:3
        other = setdiff(1:3, d);
        prof = any(any(M, other(1)), other(2));
        prof = prof(:);
        idx{d} = max(find(prof, 1, 'first') - pad, 1) : min(find(prof, 1, 'last') + pad, sz(d));
    end
end

function D = dilate6(M)
% Dilation by the voxel and its 6 face neighbours (no wrap at the edges)
    D = M;
    D(2:end, :, :)   = D(2:end, :, :)   | M(1:end-1, :, :);
    D(1:end-1, :, :) = D(1:end-1, :, :) | M(2:end, :, :);
    D(:, 2:end, :)   = D(:, 2:end, :)   | M(:, 1:end-1, :);
    D(:, 1:end-1, :) = D(:, 1:end-1, :) | M(:, 2:end, :);
    D(:, :, 2:end)   = D(:, :, 2:end)   | M(:, :, 1:end-1);
    D(:, :, 1:end-1) = D(:, :, 1:end-1) | M(:, :, 2:end);
end

function M = dilate_sparse(M, rad)
% Dilation of a sparse mask by a sphere of radius rad [voxels]
    sz = size(M, 1:3);
    rr = ceil(rad);
    [kx, ky, kz] = ndgrid(-rr:rr);
    keep = kx.^2 + ky.^2 + kz.^2 <= rad^2;
    off = [kx(keep), ky(keep), kz(keep)];
    [i, j, k] = ind2sub(sz, find(M));
    for o = 1:size(off, 1)
        ii = i + off(o, 1); jj = j + off(o, 2); kk = k + off(o, 3);
        v = ii >= 1 & ii <= sz(1) & jj >= 1 & jj <= sz(2) & kk >= 1 & kk <= sz(3);
        M(sub2ind(sz, ii(v), jj(v), kk(v))) = true;
    end
end

function [out, n_blocks] = close_at(S, C, r, use_ipt)
% Closing of S by a sphere of radius r [voxels], evaluated at the voxels of C
%
% The grid is split into blocks; each block holding a voxel of C is closed
% with a margin of 2*floor(r)+1 voxels, the reach of the closing (dilation
% then erosion, each by r), so the result inside the block equals the
% closing of the whole grid.
    sz = size(S, 1:3);
    pad = 2 * floor(r) + 1;
    B = 64;                                        % block edge [voxels]
    [ci, cj, ck] = ind2sub(sz, find(C));
    blk = unique([ceil(ci / B), ceil(cj / B), ceil(ck / B)], 'rows');
    out = false(sz);
    n_blocks = 0;
    for b = 1:size(blk, 1)
        core = cell(1, 3); ext = cell(1, 3); off = cell(1, 3);
        for d = 1:3
            core{d} = ((blk(b, d) - 1) * B + 1) : min(blk(b, d) * B, sz(d));
            ext{d}  = max(core{d}(1) - pad, 1) : min(core{d}(end) + pad, sz(d));
            off{d}  = core{d} - ext{d}(1) + 1;
        end
        Sx = S(ext{:});
        if ~any(Sx(:))
            continue                               % no skull within reach: nothing closes
        end
        n_blocks = n_blocks + 1;
        closed = morph_close_sphere(Sx, r, use_ipt);
        out(core{:}) = closed(off{:}) & C(core{:});
    end
end
