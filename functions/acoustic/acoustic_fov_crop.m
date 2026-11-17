function [parameters, kwave_medium, medium_masks, segmentation, trans_pos, focus_pos, fov_state] = ...
    acoustic_fov_crop(parameters, kwave_medium, medium_masks, segmentation, trans_pos, focus_pos)
% ACOUSTIC_FOV_CROP  Resolve FOV settings and crop the acoustic grid if active.
%
% Wraps resolve_acoustic_fov_params (auto-resolution of 'auto' sentinel values),
% jwave PML expansion, and the acoustic_grid_fov crop into a single call.
% Returns a fov_state struct that downstream stages (thermal_stage,
% nifti_acoustic) use to back-project from FOV-grid to full-grid coordinates.
%
% fov_state fields:
%   .active            - true when FOV crop was applied
%   .offset            - [1x3] voxel offset of FOV within full simulation grid
%   .full_dims         - full simulation grid dimensions before crop
%   .kwave_medium_full - kwave_medium before crop (for thermal restore)
%   .medium_masks_full - medium_masks before crop (for thermal restore)
%   .segmentation_full - segmentation before crop (for thermal restore)
%   .trans_pos_full    - trans_pos in full sim-grid voxels
%   .focus_pos_full    - focus_pos in full sim-grid voxels
%   .transducer_pos_full - struct array (trans_pos, focus_pos) of every
%                        transducer in full sim-grid voxels; restore with
%                        FOV_RESTORE_TRANSDUCERS
%
% Use as:
%   [parameters, kwave_medium, medium_masks, segmentation, ...
%    trans_pos, focus_pos, fov_state] = ...
%       acoustic_fov_crop(parameters, kwave_medium, medium_masks, ...
%                         segmentation, trans_pos, focus_pos)

arguments
    parameters  (1,1) struct
    kwave_medium (1,1) struct
    medium_masks
    segmentation
    trans_pos
    focus_pos
end

fov_state = struct('active', false, 'offset', [], 'full_dims', [], ...
    'kwave_medium_full', [], 'medium_masks_full', [], 'segmentation_full', [], ...
    'trans_pos_full',    [], 'focus_pos_full',    [], 'transducer_pos_full', []);

% No transducer/focus geometry (e.g. run_grid_setup=0 with no cache): an FOV
% crop is impossible (acoustic_grid_fov centres on the focus) and
% resolve_acoustic_fov_params rejects empty positions.  Leave the grid
% untouched, mirroring the original inline behaviour of skipping the crop.
if isempty(trans_pos) || isempty(focus_pos)
    return;
end

% Resolve 'auto' sentinel values; determine whether FOV crop is requested.
[parameters, has_acoustic_fov] = resolve_acoustic_fov_params(parameters, trans_pos, focus_pos);

if ~has_acoustic_fov
    return;
end

% jwave embeds its PML inside the domain. Expand the requested FOV by one
% PML thickness on each side so the physical FOV occupies the non-PML region.
% The thickness is resolved as jwave_wrapper does (resolve_jwave_pml_size,
% sized from the maximum sound speed); the full-grid maximum is used, which
% equals the cropped one unless the fastest voxel lies outside the FOV.
use_jwave = isfield(parameters, 'acoustic_backend') && ...
    strcmpi(parameters.acoustic_backend, 'jwave');
if use_jwave
    if isfield(parameters.grid, 'pml_size')
        jw_pml_raw = parameters.grid.pml_size;
    else
        jw_pml_raw = 'auto';
    end
    jw_pml_eff = max(resolve_jwave_pml_size(jw_pml_raw, parameters.transducer(1).freq_hz, ...
        max(kwave_medium.sound_speed(:)), parameters.grid.resolution_mm / 1e3));
    if (ischar(jw_pml_raw) || isstring(jw_pml_raw)) && strcmp(jw_pml_raw, 'auto')
        fprintf('[pipeline] jwave pfov PML expansion: pml_size=''auto'' resolved to %d voxels.\n', jw_pml_eff);
    end
    jw_pml_mm = jw_pml_eff * parameters.grid.resolution_mm;
    parameters.grid.acoustic_fov_diameter_mm = parameters.grid.acoustic_fov_diameter_mm + 2 * jw_pml_mm;
    parameters.grid.acoustic_fov_length_mm   = parameters.grid.acoustic_fov_length_mm   + 2 * jw_pml_mm;
    if isfield(parameters, 'state') && isfield(parameters.state, 'acoustic_fov_per_tx')
        parameters.state.acoustic_fov_per_tx.diameter_mm = parameters.state.acoustic_fov_per_tx.diameter_mm + 2 * jw_pml_mm;
        parameters.state.acoustic_fov_per_tx.length_mm   = parameters.state.acoustic_fov_per_tx.length_mm   + 2 * jw_pml_mm;
    end
    fprintf('[pipeline] jwave pfov PML expansion (+%.0f mm each side): diameter=%.1f mm, length=%.1f mm\n', ...
        jw_pml_mm, parameters.grid.acoustic_fov_diameter_mm, parameters.grid.acoustic_fov_length_mm);
end

% Save full-grid state before crop.
fov_state.full_dims         = parameters.grid.dims;
fov_state.kwave_medium_full = kwave_medium;
fov_state.medium_masks_full = medium_masks;
fov_state.segmentation_full = segmentation;
fov_state.trans_pos_full    = trans_pos;
fov_state.focus_pos_full    = focus_pos;
fov_state.transducer_pos_full = transducer_positions(parameters);

% Apply FOV crop.
[parameters, kwave_medium, medium_masks, trans_pos, focus_pos, fov_offset] = ...
    acoustic_grid_fov(parameters, kwave_medium, medium_masks, trans_pos, focus_pos);
fov_state.offset = fov_offset;
fov_state.active = true;

fprintf('[pipeline] Acoustic FOV crop active: offset [%s], new dims [%s]\n', ...
    num2str(fov_offset), num2str(parameters.grid.dims));

% Crop segmentation to same FOV so analysis masks match sensor_data size.
fov_end = fov_offset + parameters.grid.dims - 1;
if ~isempty(segmentation) && isequal(size(segmentation), double(fov_state.full_dims))
    segmentation = segmentation( ...
        fov_offset(1):fov_end(1), ...
        fov_offset(2):fov_end(2), ...
        fov_offset(3):fov_end(3));
end

end

function pos = transducer_positions(parameters)
% trans_pos/focus_pos of every transducer (full sim-grid voxels, before the crop).
pos = struct('trans_pos', {}, 'focus_pos', {});
for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);
    pos(ti).trans_pos = [];
    pos(ti).focus_pos = [];
    if isfield(tr, 'trans_pos'), pos(ti).trans_pos = tr.trans_pos; end
    if isfield(tr, 'focus_pos'), pos(ti).focus_pos = tr.focus_pos; end
end
end
