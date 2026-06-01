function [parameters_fov, kwave_medium_fov, medium_masks_fov, ...
          trans_pos_fov, focus_pos_fov, fov_offset_ac, full_ac_dims] = ...
    acoustic_grid_fov(parameters, kwave_medium, medium_masks, trans_pos, focus_pos)
% ACOUSTIC_GRID_FOV  Crop the acoustic simulation domain to a beam-axis FOV box.
%
% Reduces the k-Wave grid to a rectangular box centred on the focus before
% acoustic simulation.  When grid.mode = 'transducer_axis' (default), the
% grid z-axis is the acoustic beam axis; x and y are lateral.  The FOV box
% dimensions are therefore:
%
%   lateral (dims 1 & 2):  grid.acoustic_fov_diameter_mm  [mm]
%   axial   (dim  3):      grid.acoustic_fov_length_mm    [mm]
%
% The crop offsets are stored in fov_offset_ac and full_ac_dims so that the
% resulting pressure field can be back-projected onto the full acoustic grid
% by ASSEMBLE_LIMITED_FOV_FIELDS.
%
% Required parameters.grid fields:
%   acoustic_fov_diameter_mm  — lateral FOV half-width [mm]
%   acoustic_fov_length_mm    — axial FOV depth [mm]
%
% Optional parameters.grid fields:
%   mode  — 'transducer_axis' (default) or 'ras_plus'.
%           In 'ras_plus' mode the beam axis is derived from trans_pos→focus_pos
%           and a warning is issued; FOV orientation may be approximate.
%
% Use as:
%   [parameters_fov, kwave_medium_fov, medium_masks_fov, ...
%    trans_pos_fov, focus_pos_fov, fov_offset_ac, full_ac_dims] = ...
%       acoustic_grid_fov(parameters, kwave_medium, medium_masks, trans_pos, focus_pos)
%
% Input:
%   parameters   - PRESTUS config with grid.acoustic_fov_diameter_mm and
%                  grid.acoustic_fov_length_mm set
%   kwave_medium - struct with 3-D field arrays on the full acoustic grid
%   medium_masks - [Nx × Ny × Nz] uint8/int label map on the full acoustic grid
%   trans_pos    - [1×3] transducer centre in full acoustic grid voxels
%   focus_pos    - [1×3] focal point in full acoustic grid voxels
%
% Output:
%   parameters_fov   - parameters copy with grid.dims, trans_pos, focus_pos
%                      updated to the cropped FOV grid
%   kwave_medium_fov - kwave_medium fields cropped to FOV (same resolution)
%   medium_masks_fov - medium_masks cropped to FOV (integer labels preserved)
%   trans_pos_fov    - transducer position in FOV grid voxels (clamped)
%   focus_pos_fov    - focus position in FOV grid voxels
%   fov_offset_ac    - [1×3] start voxel index in the FULL acoustic grid (1-based)
%   full_ac_dims     - [1×3] full acoustic grid dimensions
%
% See also: ASSEMBLE_LIMITED_FOV_FIELDS, THERMAL_GRID_SETUP, PRESTUS_PIPELINE

arguments
    parameters   (1,1) struct
    kwave_medium (1,1) struct
    medium_masks {mustBeNumericOrLogical}
    trans_pos    (1,:) {mustBeNumeric}
    focus_pos    (1,:) {mustBeNumeric}
end

% -------------------------------------------------------------------------
% Validate config fields
% -------------------------------------------------------------------------
if ~isfield(parameters.grid, 'acoustic_fov_diameter_mm') || ...
        isempty(parameters.grid.acoustic_fov_diameter_mm)
    error('acoustic_grid_fov:missingDiameter', ...
        'parameters.grid.acoustic_fov_diameter_mm must be set.');
end
if ~isfield(parameters.grid, 'acoustic_fov_length_mm') || ...
        isempty(parameters.grid.acoustic_fov_length_mm)
    error('acoustic_grid_fov:missingLength', ...
        'parameters.grid.acoustic_fov_length_mm must be set.');
end

fov_diameter_mm = parameters.grid.acoustic_fov_diameter_mm;
fov_length_mm   = parameters.grid.acoustic_fov_length_mm;
dx_mm           = parameters.grid.resolution_mm;
full_ac_dims    = parameters.grid.dims;   % [Nx Ny Nz]
ndim            = numel(full_ac_dims);

if ndim ~= 3
    error('acoustic_grid_fov:not3D', ...
        'acoustic_grid_fov only supports 3-D grids (ndim=%d).', ndim);
end

% -------------------------------------------------------------------------
% Determine beam axis
% -------------------------------------------------------------------------
ras_plus_mode = isfield(parameters.grid, 'mode') && ...
    strcmp(parameters.grid.mode, 'ras_plus');

if ras_plus_mode
    % Derive beam axis from transducer→focus vector in grid voxels.
    beam_vec = focus_pos - trans_pos;
    [~, beam_dim] = max(abs(beam_vec));
    warning('acoustic_grid_fov:rasPlusMode', ...
        ['grid.mode=''ras_plus'': beam axis inferred as dim %d from trans→focus vector. ' ...
         'FOV orientation may be approximate.'], beam_dim);
    lateral_dims = setdiff(1:3, beam_dim);
else
    % transducer_axis mode: z (dim 3) is the beam axis.
    beam_dim     = 3;
    lateral_dims = [1 2];
end

% -------------------------------------------------------------------------
% Build per-dimension FOV [mm]: diameter for lateral, length for axial
% -------------------------------------------------------------------------
fov_mm = zeros(1, 3);
fov_mm(lateral_dims) = fov_diameter_mm;
fov_mm(beam_dim)     = fov_length_mm;

ac_fov_mm = full_ac_dims * dx_mm;   % full acoustic FOV in mm

% Clamp requested FOV to full acoustic extent
fov_mm = min(fov_mm, ac_fov_mm);

% Centre FOV on focus_pos
focus_mm  = (focus_pos - 1) * dx_mm;   % [1×3]
fov_lo_mm = focus_mm - fov_mm / 2;
fov_hi_mm = focus_mm + fov_mm / 2;

% Clamp lo/hi, then nudge lo if hi was clamped
fov_hi_mm = min(fov_hi_mm, ac_fov_mm);
fov_lo_mm = max(fov_lo_mm, 0);
fov_lo_mm = min(fov_lo_mm, fov_hi_mm - fov_mm);
fov_lo_mm = max(fov_lo_mm, 0);

% Convert to 1-based voxel indices in the full grid
fov_start_ac = max(ones(1,3),         round(fov_lo_mm / dx_mm) + 1);
fov_end_ac   = min(full_ac_dims,      round(fov_hi_mm / dx_mm) + 1);
fov_dims     = fov_end_ac - fov_start_ac + 1;   % cropped dims

fov_offset_ac = fov_start_ac;   % 1-based start index in full grid

fprintf('[acoustic_grid_fov] Full acoustic grid: %.2f mm, [%s] vox\n', ...
    dx_mm, num2str(full_ac_dims));
fprintf('[acoustic_grid_fov] FOV box: lateral=%.1f mm, axial=%.1f mm  →  [%s] vox  (offset [%s])\n', ...
    fov_diameter_mm, fov_length_mm, num2str(fov_dims), num2str(fov_start_ac));

% -------------------------------------------------------------------------
% Crop kwave_medium fields
% -------------------------------------------------------------------------
kwave_medium_fov = kwave_medium;
field_names = fieldnames(kwave_medium);
for fi = 1:numel(field_names)
    fname = field_names{fi};
    vol   = kwave_medium.(fname);
    if ~isnumeric(vol) || ~isequal(size(vol), double(full_ac_dims))
        continue
    end
    kwave_medium_fov.(fname) = vol( ...
        fov_start_ac(1):fov_end_ac(1), ...
        fov_start_ac(2):fov_end_ac(2), ...
        fov_start_ac(3):fov_end_ac(3));
end

% -------------------------------------------------------------------------
% Crop medium_masks (nearest-neighbour: just array indexing)
% -------------------------------------------------------------------------
medium_masks_fov = medium_masks( ...
    fov_start_ac(1):fov_end_ac(1), ...
    fov_start_ac(2):fov_end_ac(2), ...
    fov_start_ac(3):fov_end_ac(3));

% -------------------------------------------------------------------------
% Update trans_pos and focus_pos to FOV grid voxel coordinates
% -------------------------------------------------------------------------
fov_origin_mm = (fov_start_ac - 1) * dx_mm;   % physical offset of FOV origin

% ac_vox → fov_vox:  fov_v = (ac_v - 1)*dx - fov_origin + 1  (in mm/dx units)
trans_pos_fov = round(((trans_pos  - 1) * dx_mm - fov_origin_mm) / dx_mm + 1);
focus_pos_fov = round(((focus_pos  - 1) * dx_mm - fov_origin_mm) / dx_mm + 1);

% Clamp (transducer may sit outside the FOV)
trans_pos_fov = max(1, min(trans_pos_fov, fov_dims));
focus_pos_fov = max(1, min(focus_pos_fov, fov_dims));

% -------------------------------------------------------------------------
% Clone parameters and update for FOV grid
% -------------------------------------------------------------------------
parameters_fov            = parameters;
parameters_fov.grid.dims  = fov_dims;
parameters_fov.transducer(1).trans_pos = trans_pos_fov;
if isfield(parameters_fov.transducer(1), 'focus_pos')
    parameters_fov.transducer(1).focus_pos = focus_pos_fov;
end

end
