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
% A beam more than 2 degrees off the grid axis gets a box around its own beam
% axis instead (see below).
%
% With several transducers (coherent coupling) the box is the union of one
% box per transducer, each enclosing a cylinder of the FOV diameter that runs
% the FOV length from that transducer's centre along its own beam axis, and
% every transducer's trans_pos/focus_pos is moved into FOV voxels.
%
% The crop offsets are stored in fov_offset_ac and full_ac_dims so that the
% resulting pressure field can be back-projected onto the full acoustic grid
% by ASSEMBLE_LIMITED_FOV_FIELDS.
%
% Required parameters.grid fields:
%   acoustic_fov_diameter_mm  — lateral FOV diameter [mm] (numeric; resolve 'auto' in caller)
%   acoustic_fov_length_mm    — axial FOV depth [mm]      (numeric; resolve 'auto' in caller)
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
%   parameters_fov   - parameters copy with grid.dims and every transducer's
%                      trans_pos, focus_pos updated to the cropped FOV grid
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
% Build FOV extents [mm]
%
% Axial  (beam_dim):    starts at the transducer centre, extends fov_length_mm
%                       forward along the beam axis.
% Lateral (lateral_dims): centred on the focus, width = fov_diameter_mm.
% -------------------------------------------------------------------------
ac_fov_mm   = full_ac_dims * dx_mm;   % full acoustic extent in mm
trans_mm    = (trans_pos  - 1) * dx_mm;   % [1×3]
focus_mm    = (focus_pos  - 1) * dx_mm;   % [1×3]

fov_lo_mm = zeros(1, 3);
fov_hi_mm = zeros(1, 3);

% Axial: transducer centre → focus direction for fov_length_mm.
% Extend in the direction the beam actually travels (sign of focus−trans),
% not always in +beam_dim, so the FOV contains the focus even when the
% transducer sits at a higher coordinate than the focal point.
if focus_mm(beam_dim) >= trans_mm(beam_dim)
    fov_lo_mm(beam_dim) = trans_mm(beam_dim);
    fov_hi_mm(beam_dim) = trans_mm(beam_dim) + fov_length_mm;
else
    fov_lo_mm(beam_dim) = trans_mm(beam_dim) - fov_length_mm;
    fov_hi_mm(beam_dim) = trans_mm(beam_dim);
end

% Lateral: centred on focus
fov_lo_mm(lateral_dims) = focus_mm(lateral_dims) - fov_diameter_mm / 2;
fov_hi_mm(lateral_dims) = focus_mm(lateral_dims) + fov_diameter_mm / 2;

% Several transducers (coherent coupling): the box must hold every beam, so
% it becomes the union of one box per transducer (see MULTI_TRANSDUCER_BOX).
% A single beam that is oblique to the grid axis (e.g. transducer 2 on an
% async shared grid aligned to transducer 1) gets the same beam-axis box:
% the box above, centred laterally on the focus, would cut off its bowl.
n_tx = numel(parameters.transducer);
u_beam  = (focus_pos - trans_pos) / max(norm(focus_pos - trans_pos), eps);
oblique = ~ras_plus_mode && acosd(min(1, abs(u_beam(beam_dim)))) > 2;
if n_tx > 1 || oblique
    [fov_lo_mm, fov_hi_mm] = multi_transducer_box(parameters, trans_pos, focus_pos, ...
        fov_diameter_mm, fov_length_mm, dx_mm);
end

% Clamp to full acoustic grid
fov_lo_mm = max(fov_lo_mm, 0);
fov_hi_mm = min(fov_hi_mm, ac_fov_mm);

% Convert to 1-based voxel indices in the full grid
fov_start_ac = max(ones(1,3),    round(fov_lo_mm / dx_mm) + 1);
fov_end_ac   = min(full_ac_dims, round(fov_hi_mm / dx_mm) + 1);
fov_dims     = fov_end_ac - fov_start_ac + 1;   % cropped dims

fov_offset_ac = fov_start_ac;   % 1-based start index in full grid

fprintf('[acoustic_grid_fov] Full acoustic grid: [%s] vox @ %.2f mm\n', ...
    num2str(full_ac_dims), dx_mm);
if n_tx > 1
    fprintf('[acoustic_grid_fov] FOV box = union of %d transducer beam boxes\n', n_tx);
elseif oblique
    fprintf('[acoustic_grid_fov] Beam %.0f deg off the grid axis: FOV box around the beam axis\n', ...
        acosd(min(1, abs(u_beam(beam_dim)))));
end
fprintf('[acoustic_grid_fov] FOV box: lateral=%.1f mm, axial=%.1f mm  ->  [%s] vox  (offset [%s])\n', ...
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

% Focus must be inside the FOV; warn if it was clamped (fov_length_mm too short).
focus_pos_fov_unclamped = focus_pos_fov;
trans_pos_fov = max(1, min(trans_pos_fov, fov_dims));
focus_pos_fov = max(1, min(focus_pos_fov, fov_dims));
if any(focus_pos_fov ~= focus_pos_fov_unclamped)
    warning('acoustic_grid_fov:focusOutsideFOV', ...
        ['Focus (full-grid vox [%s]) maps to FOV vox [%s], clamped from [%s]. ' ...
         'Increase acoustic_fov_length_mm so the FOV reaches beyond the focus.'], ...
        num2str(focus_pos), num2str(focus_pos_fov), num2str(focus_pos_fov_unclamped));
end

% -------------------------------------------------------------------------
% Clone parameters and update for FOV grid
% -------------------------------------------------------------------------
parameters_fov            = parameters;
parameters_fov.grid.dims  = fov_dims;
parameters_fov.transducer(1).trans_pos = trans_pos_fov;
if isfield(parameters_fov.transducer(1), 'focus_pos')
    parameters_fov.transducer(1).focus_pos = focus_pos_fov;
end

% Transducers 2..N: source_sensor_setup reads their positions from
% parameters.transducer(ti), so they must move into FOV voxels as well.
for ti = 2:n_tx
    for key = {'trans_pos', 'focus_pos'}
        pos = parameters.transducer(ti).(key{1});
        pos_fov = round(((pos(:)' - 1) * dx_mm - fov_origin_mm) / dx_mm + 1);
        if any(pos_fov < 1 | pos_fov > fov_dims)
            warning('acoustic_grid_fov:transducerOutsideFOV', ...
                ['Transducer %d %s (full-grid vox [%s]) lies outside the FOV box [%s] ' ...
                 '(offset [%s]); it is clamped to the box edge.'], ti, key{1}, ...
                num2str(pos(:)'), num2str(fov_dims), num2str(fov_start_ac));
            pos_fov = max(1, min(pos_fov, fov_dims));
        end
        parameters_fov.transducer(ti).(key{1}) = pos_fov;
    end
end

end

% -------------------------------------------------------------------------
function [lo_mm, hi_mm] = multi_transducer_box(parameters, trans_pos1, focus_pos1, ...
                                               diameter_mm, length_mm, dx_mm)
% Union of one box per transducer, in full-grid mm (0 = first voxel). Each
% transducer's box encloses a cylinder of diameter_mm that starts at the
% transducer centre and runs length_mm along its own trans->focus axis, so
% it does not rely on the grid z axis being any one beam's axis (focal_axis
% mode aligns z to the MEAN axis of all transducers). Per-transducer
% 'auto' values resolved by RESOLVE_ACOUSTIC_FOV_PARAMS, when present,
% replace the scalar diameter/length for that transducer.
n_tx = numel(parameters.transducer);
diam = repmat(diameter_mm, 1, n_tx);
len  = repmat(length_mm,   1, n_tx);
if isfield(parameters, 'state') && isfield(parameters.state, 'acoustic_fov_per_tx')
    per = parameters.state.acoustic_fov_per_tx;
    if numel(per.diameter_mm) == n_tx, diam = per.diameter_mm; end
    if numel(per.length_mm)   == n_tx, len  = per.length_mm;   end
end

lo_mm = inf(1, 3);
hi_mm = -inf(1, 3);
for ti = 1:n_tx
    if ti == 1
        t = trans_pos1(:)'; f = focus_pos1(:)';
    else
        t = parameters.transducer(ti).trans_pos(:)';
        f = parameters.transducer(ti).focus_pos(:)';
    end
    t_mm = (t - 1) * dx_mm;
    u    = (f - t) / max(norm(f - t), eps);
    e_mm = t_mm + len(ti) * u;
    % half-extent of a disc of radius r with normal u, per axis
    r_ax = diam(ti) / 2 * sqrt(max(0, 1 - u.^2));
    lo_mm = min(lo_mm, min(t_mm, e_mm) - r_ax);
    hi_mm = max(hi_mm, max(t_mm, e_mm) + r_ax);
end
end
