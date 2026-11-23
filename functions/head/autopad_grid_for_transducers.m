function [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = ...
    autopad_grid_for_transducers(parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg)
% AUTOPAD_GRID_FOR_TRANSDUCERS  Expand grid if any transducer/focus position falls outside
%
% After GRID_TISSUE_SETUP produces a cropped grid (potentially sized for a
% single transducer), this function checks whether all transducer and focus
% positions would land within [1, grid.dims] after applying planimg.transf.
% If any position exceeds the grid in any dimension, the minimum required
% padding is computed and applied symmetrically, the volumetric arrays are
% expanded, and planimg.transf / planimg.inv_transf are updated accordingly.
%
% A note is printed whenever padding is applied, recommending that the user
% set headmodel.head_pad_mm to avoid this auto-correction on future runs.
%
% Only active for 'layered' simulations (where planimg.transf is defined).
% For non-layered media the function is a no-op.
%
% Use as:
%   [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = ...
%       autopad_grid_for_transducers(parameters, medium_masks, segmentation, ...
%                                    bone_mask, pseudoCT, planimg)
%
% See also: PREPROC_CROP_GRID, GRID_TRANSDUCER_LOCATION, PRESTUS_PIPELINE

arguments
    parameters  (1,1) struct
    medium_masks
    segmentation
    bone_mask
    pseudoCT
    planimg     (1,1) struct
end

% Only layered simulations use an affine transform; skip others.
if ~contains(parameters.simulation.medium, 'layered')
    return
end
if ~isfield(planimg, 'transf') || isempty(planimg.transf)
    return
end

N = numel(parameters.transducer);
if N == 0; return; end

% Collect all T1-space positions that grid_transducer_location will transform.
all_pts = zeros(2*N, 3);
for ti = 1:N
    tr = parameters.transducer(ti);
    if ~isfield(tr,'trans_pos') || isempty(tr.trans_pos) || ...
       ~isfield(tr,'focus_pos') || isempty(tr.focus_pos)
        continue
    end
    all_pts(2*ti-1, :) = tr.trans_pos(:)';
    all_pts(2*ti,   :) = tr.focus_pos(:)';
end

% Compute proposed sim-grid positions.
proposed = round(affine_apply_pts(all_pts, planimg.transf));

dims = parameters.grid.dims(:)';  % [1×3]

% Per-dimension padding needed:
%   pre_pad(d)  — positions < 1 require this many extra voxels on the low side
%   post_pad(d) — positions > dims(d) require this many extra on the high side
pre_pad  = max(0, 1        - min(proposed, [], 1));   % 1×3
post_pad = max(0, max(proposed, [], 1) - dims);       % 1×3

if ~any(pre_pad) && ~any(post_pad)
    return  % all positions within bounds — nothing to do
end

% Build a per-transducer description for the note.
msgs = {};
for ti = 1:N
    tp_prop = proposed(2*ti-1, :);
    fp_prop = proposed(2*ti,   :);
    tp_over = (tp_prop < 1) | (tp_prop > dims);
    fp_over = (fp_prop < 1) | (fp_prop > dims);
    if any(tp_over)
        msgs{end+1} = sprintf('  Transducer %d trans_pos %s → grid [%d %d %d] (overflow dims: %s)', ...
            ti, mat2str(tp_prop), dims(1), dims(2), dims(3), ...
            mat2str(find(tp_over)));  %#ok<AGROW>
    end
    if any(fp_over)
        msgs{end+1} = sprintf('  Transducer %d focus_pos %s → grid [%d %d %d] (overflow dims: %s)', ...
            ti, mat2str(fp_prop), dims(1), dims(2), dims(3), ...
            mat2str(find(fp_over)));  %#ok<AGROW>
    end
end

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  NOTE [autopad]: transducer position(s) exceed grid bounds   ║\n');
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
for k = 1:numel(msgs)
    fprintf('║  %s\n', msgs{k});
end
fprintf('║\n');
fprintf('║  Auto-applying padding:\n');
fprintf('║    pre-pad  [%d %d %d] voxels (%.1f mm)\n', pre_pad,  max(pre_pad)  * parameters.grid.resolution_mm);
fprintf('║    post-pad [%d %d %d] voxels (%.1f mm)\n', post_pad, max(post_pad) * parameters.grid.resolution_mm);
total_pad_mm = (max(pre_pad) + max(post_pad)) * parameters.grid.resolution_mm;
fprintf('║\n');
fprintf('║  To avoid auto-padding, set headmodel.head_pad_mm ≥ %.0f mm\n', total_pad_mm);
fprintf('║  in your config (currently: %.0f mm).\n', ...
    get_head_pad_mm(parameters));
fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

% Determine the water label (pad value for medium_masks).
medium_names = fieldnames(parameters.medium_properties);
i_water = find(strcmp(medium_names, 'water'));
if isempty(i_water); i_water = 0; end

% Apply pre-padding (low side of each dimension).
if any(pre_pad > 0)
    medium_masks = pad3_pre(medium_masks, pre_pad, i_water);
    segmentation = pad3_pre(segmentation, pre_pad, 0);
    bone_mask    = pad3_pre(bone_mask,    pre_pad, 0);
    if ~isempty(pseudoCT)
        pseudoCT = pad3_pre(pseudoCT, pre_pad, NaN);
    end
    % Shift the translation component of planimg.transf.
    % Row-vector convention: [out 1] = [in 1] * M  →  translation in M(4, 1:3).
    planimg.transf(4, 1:3) = planimg.transf(4, 1:3) + pre_pad;
end

% Apply post-padding (high side of each dimension).
if any(post_pad > 0)
    medium_masks = pad3_post(medium_masks, post_pad, i_water);
    segmentation = pad3_post(segmentation, post_pad, 0);
    bone_mask    = pad3_post(bone_mask,    post_pad, 0);
    if ~isempty(pseudoCT)
        pseudoCT = pad3_post(pseudoCT, post_pad, NaN);
    end
    % Post-padding does not shift existing voxels, so no transform update needed.
end

% Update grid dimensions and inverse transform.
parameters.grid.dims = size(medium_masks);
planimg.inv_transf   = affine_inverse(planimg.transf);

fprintf('[autopad] Grid expanded: [%d %d %d] → [%d %d %d]\n', ...
    dims, parameters.grid.dims);

end

% =========================================================================
%% Helpers
% =========================================================================

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

function mm = get_head_pad_mm(parameters)
if isfield(parameters, 'headmodel') && isfield(parameters.headmodel, 'head_pad_mm') && ...
        ~isempty(parameters.headmodel.head_pad_mm)
    mm = parameters.headmodel.head_pad_mm;
else
    mm = 0;
end
end
