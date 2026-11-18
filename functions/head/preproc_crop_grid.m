function [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, ...
         focus_pos_final, translation_matrix] = preproc_crop_grid( ...
         parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img, all_trans_pos_grid, all_focus_pos_grid)
% PREPROC_CROP_GRID  Crop 3-D head model to the simulation-relevant region
%
% Crops the full-head model to a CSF-guided, transducer-inclusive bounding
% box, adds symmetric user padding, optimises dimensions for FFT, and
% applies conditional pre/post padding to prevent out-of-bounds access.
% The PML is always applied outside the returned grid (PMLInside=false).
% Symmetric padding (parameters.headmodel.head_pad_mm) is applied BEFORE
% transducer_setup to prevent out-of-bounds access. Grid dimensions are
% then adjusted to be FFT-efficient. Skull edges are computed to derive
% the crop bounding box.
%
% Use as:
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, ...
%    focus_pos_final, translation_matrix] = preproc_crop_grid( ...
%       parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid)
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, ...
%    focus_pos_final, translation_matrix] = preproc_crop_grid( ...
%       parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img)
%   [...] = preproc_crop_grid(..., pseudoCT_img, all_trans_pos_grid, all_focus_pos_grid)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   medium_masks   - [Nx x Ny x Nz] medium label array
%   segmentation   - [Nx x Ny x Nz] tissue segmentation
%   bone_mask_img  - [Nx x Ny x Nz] binary skull mask
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   pseudoCT_img   - [Nx x Ny x Nz] Hounsfield-unit skull image (optional, [] if unused)
%   all_trans_pos_grid - [N x 3] positions of all transducers, kept inside the
%                        crop (optional, default: zeros(0,3))
%   all_focus_pos_grid - [N x 3] focus positions of all transducers (optional)
%
% Output:
%   medium_masks       - cropped medium label array
%   segmentation_crop  - cropped tissue segmentation
%   bone_mask_crop     - cropped binary skull mask
%   pseudoCT_crop      - cropped pseudoCT image ([] when pseudoCT_img is not provided)
%   parameters         - updated struct with .grid.dims = [Nx Ny Nz]
%   trans_pos_final    - [1x3] transducer position in cropped grid
%   focus_pos_final    - [1x3] focus position in cropped grid
%   translation_matrix - [4x4] homogeneous translation matrix
%
% See also: PREPROC_CROP_ECSF, HEAD_SMOOTH_AND_CROP

if nargin < 7; pseudoCT_img = []; end
if nargin < 8; all_trans_pos_grid = zeros(0,3); end
if nargin < 9; all_focus_pos_grid = zeros(0,3); end

% === USER SYMMETRIC PADDING (BOTH SIDES) ===
if ~isfield(parameters.headmodel, 'head_pad_mm') || isempty(parameters.headmodel.head_pad_mm)
    parameters.headmodel.head_pad_mm = 0;  % mm
end
pad_voxels = round(parameters.headmodel.head_pad_mm / parameters.grid.resolution_mm);
pad_pre_post = [pad_voxels, pad_voxels, pad_voxels];  % Symmetric BOTH sides

% Track TOTAL pre-padding offset for translation matrix
total_pre_offset = zeros(1,3);

% Apply symmetric padding BEFORE transducer_setup
if any(pad_pre_post > 0)
    segmentation    = pad3_both(segmentation,  pad_pre_post, 0);
    medium_masks    = pad3_both(medium_masks,  pad_pre_post, 0);
    bone_mask_img   = pad3_both(bone_mask_img, pad_pre_post, 0);
    if ~isempty(pseudoCT_img)
        pseudoCT_img = pad3_both(pseudoCT_img, pad_pre_post, NaN);
    end
    total_pre_offset = total_pre_offset + pad_pre_post;  % User padding contribution
    % Positions shift by PRE-padding amount ('both' adds pre first)
    trans_pos_grid = trans_pos_grid + pad_pre_post;
    focus_pos_grid = focus_pos_grid + pad_pre_post;
    if ~isempty(all_trans_pos_grid)
        all_trans_pos_grid = all_trans_pos_grid + pad_pre_post;
        all_focus_pos_grid = all_focus_pos_grid + pad_pre_post;
    end
end

% Expand CSF mask to crop the layered medium
[medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmentation, trans_pos_grid, focus_pos_grid);

% Unlabelled voxels (padding, resampling fill) are water, as in MEDIUM_SETUP,
% so they do not widen the crop
i_water = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
medium_masks(medium_masks == 0) = i_water;

% Minimal fixed crop margin — PML is added outside the grid by k-Wave (PMLInside=false)
crop_margin = 1;

% Include transducer bowl geometry for all transducers (safe with padding).
% For multi-transducer setups, using only transducer(1) leaves the contralateral
% side over-extended: the CSF-expanded head tissue reaches the far scalp + 40 mm,
% while the opposite transducer bowl is never included in the crop mask.

% Use all_trans_pos_grid when provided (Nt×3), otherwise fall back to trans(1) only.
if ~isempty(all_trans_pos_grid)
    n_bowls = size(all_trans_pos_grid, 1);
else
    n_bowls = 1;
end

% === AUTO-PAD: guarantee every transducer centre is inside the volume ===
% makeBowl hard-errors when the centre falls outside the grid.  This can
% happen when head_pad_mm is insufficient for a given anatomy/placement.
% Compute the minimum extra padding needed per face and apply it now.
all_centres = zeros(n_bowls, 3);
for i_t_chk = 1:n_bowls
    if ~isempty(all_trans_pos_grid)
        all_centres(i_t_chk,:) = all_trans_pos_grid(i_t_chk,:);
    else
        all_centres(i_t_chk,:) = trans_pos_grid;
    end
end
cur_dims  = size(medium_masks);
need_pre  = max(0, ceil(1 - min(all_centres, [], 1)));   % shortfall at low faces
need_post = max(0, ceil(max(all_centres, [], 1) - cur_dims));  % shortfall at high faces
if any(need_pre > 0) || any(need_post > 0)
    need_pre_mm  = need_pre  * parameters.grid.resolution_mm;
    need_post_mm = need_post * parameters.grid.resolution_mm;
    warning(['preproc_crop_grid: transducer centre(s) outside volume after head_pad_mm=%.1f mm. '  ...
             'Auto-adding [%d %d %d] vox pre / [%d %d %d] vox post to fit bowl. '                  ...
             '(%.1f / %.1f / %.1f mm pre; %.1f / %.1f / %.1f mm post). '                           ...
             'Consider increasing headmodel.head_pad_mm by ~%.0f mm.'], ...
            parameters.headmodel.head_pad_mm, ...
            need_pre(1), need_pre(2), need_pre(3), ...
            need_post(1), need_post(2), need_post(3), ...
            need_pre_mm(1), need_pre_mm(2), need_pre_mm(3), ...
            need_post_mm(1), need_post_mm(2), need_post_mm(3), ...
            max([need_pre_mm, need_post_mm]));
    if any(need_pre > 0)
        segmentation    = pad3_pre(segmentation,    need_pre, 0);
        medium_masks    = pad3_pre(medium_masks,    need_pre, i_water);
        bone_mask_img   = pad3_pre(bone_mask_img,   need_pre, 0);
        if ~isempty(pseudoCT_img)
            pseudoCT_img = pad3_pre(pseudoCT_img,   need_pre, NaN);
        end
        total_pre_offset    = total_pre_offset + need_pre;
        trans_pos_grid      = trans_pos_grid  + need_pre;
        focus_pos_grid      = focus_pos_grid  + need_pre;
        if ~isempty(all_trans_pos_grid)
            all_trans_pos_grid  = all_trans_pos_grid  + need_pre;
            all_focus_pos_grid  = all_focus_pos_grid  + need_pre;
        end
    end
    if any(need_post > 0)
        segmentation    = pad3_post(segmentation,    need_post, 0);
        medium_masks    = pad3_post(medium_masks,    need_post, i_water);
        bone_mask_img   = pad3_post(bone_mask_img,   need_post, 0);
        if ~isempty(pseudoCT_img)
            pseudoCT_img = pad3_post(pseudoCT_img,   need_post, NaN);
        end
    end
end
% === END AUTO-PAD ===

combinedmask = (medium_masks ~= i_water);
for i_t = 1:n_bowls
    if ~isempty(all_trans_pos_grid)
        tp_i = all_trans_pos_grid(i_t,:);
        fp_i = all_focus_pos_grid(i_t,:);
    else
        tp_i = trans_pos_grid;
        fp_i = focus_pos_grid;
    end
    bowl_i = transducer_setup(parameters.transducer(i_t), tp_i, fp_i, ...
        size(segmentation), parameters.grid.resolution_mm, parameters);
    combinedmask = combinedmask | bowl_i;
end

% Compute crop bounds
orig_dims = size(medium_masks);
[min_dims, max_dims, new_grid_dims] = get_crop_dims(double(combinedmask), crop_margin);
clear combinedmask;

% === CONDITIONAL PRE-PADDING if min_dims < 1 ===
if any(min_dims < 1)
    pad_amount = abs(min(min_dims, [1 1 1]));
    segmentation  = pad3_pre(segmentation,  pad_amount, 0);
    medium_masks  = pad3_pre(medium_masks,  pad_amount, i_water);
    bone_mask_img = pad3_pre(bone_mask_img, pad_amount, 0);
    if ~isempty(pseudoCT_img)
        pseudoCT_img = pad3_pre(pseudoCT_img, pad_amount, NaN);
    end
    total_pre_offset = total_pre_offset + pad_amount;  % Add conditional padding
    min_dims = max(min_dims, [1 1 1]);
    max_dims = max_dims + pad_amount;
    new_grid_dims = max_dims - min_dims + 1;
    trans_pos_grid = trans_pos_grid + pad_amount;
    focus_pos_grid = focus_pos_grid + pad_amount;
end

% Optimize dimensions for FFT
new_grid_dims(1) = find_min_factor(new_grid_dims(1), new_grid_dims(1) + parameters.grid.max_expand);
new_grid_dims(2) = find_min_factor(new_grid_dims(2), new_grid_dims(2) + parameters.grid.max_expand);
new_grid_dims(3) = find_min_factor(new_grid_dims(3), new_grid_dims(3) + parameters.grid.max_expand);
max_dims = min_dims + new_grid_dims - 1;

% === CONDITIONAL POST-PADDING if FFT expansion exceeds padding ===
if any(max_dims > size(medium_masks))
    pad_post_amount = max(0, max_dims - size(medium_masks));
    segmentation  = pad3_post(segmentation,  pad_post_amount, 0);
    medium_masks  = pad3_post(medium_masks,  pad_post_amount, i_water);
    bone_mask_img = pad3_post(bone_mask_img, pad_post_amount, 0);
    if ~isempty(pseudoCT_img)
        pseudoCT_img = pad3_post(pseudoCT_img, pad_post_amount, NaN);
    end
    fprintf('Post-padding applied: [%d %d %d] voxels\n', pad_post_amount);
end

% Apply final crop
medium_masks    = medium_masks(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
segmentation_crop = segmentation(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
bone_mask_crop  = bone_mask_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
if ~isempty(pseudoCT_img)
    pseudoCT_crop = pseudoCT_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
else
    pseudoCT_crop = [];
end

% Update parameters
parameters.grid.dims = size(medium_masks);

% Final positions in cropped grid coordinates
trans_pos_final = trans_pos_grid - min_dims + 1;
focus_pos_final = focus_pos_grid - min_dims + 1;

% === COMPLETE TRANSLATION MATRIX ===
% Maps: original_world → final_cropped_world
% Offset = total_pre_padding + crop_offset
translation_matrix = makehgtform('translate', total_pre_offset + (1-min_dims));

seg_qc('add', 'grid_crop', [], medium_masks, focus_pos_final, struct( ...
    'title', 'Grid crop', ...
    'description', sprintf(['Cropped from [%d %d %d] to [%d %d %d] voxels around the head and all ' ...
    'transducers (padding %g mm, then enlarged to FFT-efficient sizes). This is the medium the ' ...
    'simulation uses.'], orig_dims, parameters.grid.dims, parameters.headmodel.head_pad_mm), ...
    'grid', seg_qc('grid_note', parameters, 'cropped'), 'label_space', 'medium'));

% Display summary
fprintf('Grid: [%dx%dx%d] → [%dx%dx%d]\n', orig_dims, parameters.grid.dims);
fprintf('  User padding: %.1fmm = [%d %d %d] voxels\n', ...
        parameters.headmodel.head_pad_mm, pad_pre_post);
fprintf('  Total pre-offset: [%d %d %d] voxels\n', total_pre_offset);
fprintf('  Crop bounds: min=[%d %d %d], max=[%d %d %d]\n', min_dims, max_dims);

end

function B = pad3_both(A, p, val)
new_sz = size(A) + 2*p;
B = fill_val(new_sz, val, class(A));
B(p(1)+1:end-p(1), p(2)+1:end-p(2), p(3)+1:end-p(3)) = A;
end

function B = pad3_pre(A, p, val)
new_sz = size(A) + p;
B = fill_val(new_sz, val, class(A));
B(p(1)+1:end, p(2)+1:end, p(3)+1:end) = A;
end

function B = pad3_post(A, p, val)
new_sz = size(A) + p;
B = fill_val(new_sz, val, class(A));
B(1:end-p(1), 1:end-p(2), 1:end-p(3)) = A;
end

function B = fill_val(sz, val, cls)
if isnan(val)
    B = NaN(sz, cls);
else
    B = repmat(cast(val, cls), sz);
end
end