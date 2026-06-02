function [kgrid_th, kwave_medium_th, medium_masks_th, parameters_th, transf_th] = ...
    thermal_grid_setup(parameters, planimg, kwave_medium, medium_masks, trans_pos, focus_pos)
% THERMAL_GRID_SETUP  Build an independent k-Wave grid for thermal simulation
%
% Creates a kWaveGrid at a (potentially coarser) thermal resolution and
% resamples all medium properties and masks from the acoustic grid into it.
% When grid.thermal_resolution_mm equals grid.resolution_mm (or is unset),
% the returned grid is identical to the acoustic one — no resampling is done.
%
% The thermal FOV can be independent of the acoustic FOV:
%   - Default: same physical extent as the acoustic grid (full-domain coupling)
%   - With grid.thermal_fov_mm: a box of that size centred on focus_pos
%
% Use as:
%   [kgrid_th, kwave_medium_th, medium_masks_th, parameters_th, transf_th] = ...
%       thermal_grid_setup(parameters, planimg, kwave_medium, medium_masks, trans_pos, focus_pos)
%
% Input:
%   parameters   - PRESTUS config; may contain:
%                    grid.thermal_resolution_mm [mm]  (default: grid.resolution_mm)
%                    grid.thermal_fov_mm        [1x3 mm] (default: full acoustic FOV)
%   planimg      - struct with fields transf (T1→acoustic-grid affine, 4×4)
%   kwave_medium - struct with 3-D field arrays on the acoustic grid
%                  (sound_speed, density, alpha_coeff, alpha_power,
%                   thermal_conductivity, specific_heat, perfusion_coeff,
%                   absorption_fraction, temp_0)
%   medium_masks - 3-D integer label map on the acoustic grid
%   trans_pos    - [1x3] transducer position in acoustic grid voxels
%   focus_pos    - [1x3] focus position in acoustic grid voxels
%
% Output:
%   kgrid_th        - kWaveGrid at thermal resolution / FOV
%   kwave_medium_th - kwave_medium fields resampled to thermal grid
%   medium_masks_th - medium_masks resampled to thermal grid (nearest)
%   parameters_th   - parameters copy with grid.dims and transducer(1).trans_pos
%                     updated to thermal grid coordinates
%   transf_th       - 4×4 affine: T1 voxels → thermal grid voxels
%                     (for resampling NIfTIs into the thermal grid)
%
% See also: THERMAL_SIMULATION, SOURCE_SENSOR_SETUP

arguments
    parameters   (1,1) struct
    planimg      (1,1) struct
    kwave_medium (1,1) struct
    medium_masks {mustBeNumericOrLogical}
    trans_pos    (1,:) {mustBeNumeric}
    focus_pos    (1,:) {mustBeNumeric}
end

ac_dx_mm  = parameters.grid.resolution_mm;
ac_dims   = parameters.grid.dims;    % [N1 N2 N3]
ndim      = numel(ac_dims);

% -------------------------------------------------------------------------
% Thermal resolution
% -------------------------------------------------------------------------
if isfield(parameters.grid, 'thermal_resolution_mm') && ...
        ~isempty(parameters.grid.thermal_resolution_mm)
    th_dx_mm = parameters.grid.thermal_resolution_mm;
else
    th_dx_mm = ac_dx_mm;
end
th_dx_m = th_dx_mm / 1e3;
scale    = ac_dx_mm / th_dx_mm;   % >1 when coarsening

% -------------------------------------------------------------------------
% Thermal FOV: default = same physical extent as acoustic grid
% -------------------------------------------------------------------------
ac_fov_mm = ac_dims * ac_dx_mm;   % [1×ndim] full acoustic FOV in mm

if isfield(parameters.grid, 'thermal_fov_mm') && ...
        ~isempty(parameters.grid.thermal_fov_mm)
    th_fov_mm = parameters.grid.thermal_fov_mm(:)';   % [1×ndim]
    if numel(th_fov_mm) ~= ndim
        error('thermal_grid_setup:dimMismatch', ...
            'grid.thermal_fov_mm must have %d elements (one per grid dimension).', ndim);
    end
    % Clamp to acoustic FOV
    th_fov_mm = min(th_fov_mm, ac_fov_mm);

    % Centre the thermal FOV on the focus position (in acoustic grid voxels).
    % Convert focus to mm from acoustic grid origin.
    focus_mm  = (focus_pos - 1) * ac_dx_mm;          % [1×ndim]
    fov_lo_mm = max(0,          focus_mm - th_fov_mm/2);
    fov_hi_mm = min(ac_fov_mm,  focus_mm + th_fov_mm/2);
    % Re-centre if clamped
    fov_lo_mm = min(fov_lo_mm, fov_hi_mm - th_fov_mm);
    fov_lo_mm = max(fov_lo_mm, 0);

    % Convert FOV bounds to acoustic voxel indices (1-based)
    fov_start_ac = max(1,       round(fov_lo_mm / ac_dx_mm) + 1);
    fov_end_ac   = min(ac_dims, round(fov_hi_mm / ac_dx_mm) + 1);
    cropped_ac_dims = fov_end_ac - fov_start_ac + 1;
else
    % Full acoustic FOV
    fov_start_ac    = ones(1, ndim);
    fov_end_ac      = ac_dims;
    cropped_ac_dims = ac_dims;
    th_fov_mm       = ac_fov_mm;
end

% Thermal grid dimensions (at least 1 voxel per axis)
th_dims = max(1, round(cropped_ac_dims * scale));

% -------------------------------------------------------------------------
% Short-circuit: if grids are identical, skip all resampling
% -------------------------------------------------------------------------
same_grid = (th_dx_mm == ac_dx_mm) && all(fov_start_ac == 1) && all(th_dims == ac_dims);

if same_grid
    if ndim == 3
        kgrid_th = kWaveGrid(ac_dims(1), th_dx_m, ac_dims(2), th_dx_m, ac_dims(3), th_dx_m);
    else
        kgrid_th = kWaveGrid(ac_dims(1), th_dx_m, ac_dims(2), th_dx_m);
    end
    kwave_medium_th = kwave_medium;
    medium_masks_th = medium_masks;
    parameters_th   = parameters;
    transf_th       = planimg.transf;
    fprintf('[thermal_grid_setup] Thermal grid identical to acoustic grid (%s mm, %s vox).\n', ...
        num2str(th_dx_mm), num2str(th_dims));
    return
end

% -------------------------------------------------------------------------
% Build kWaveGrid at thermal resolution
% -------------------------------------------------------------------------
if ndim == 3
    kgrid_th = kWaveGrid(th_dims(1), th_dx_m, th_dims(2), th_dx_m, th_dims(3), th_dx_m);
else
    kgrid_th = kWaveGrid(th_dims(1), th_dx_m, th_dims(2), th_dx_m);
end

fprintf('[thermal_grid_setup] Acoustic grid: %.2f mm, [%s] vox\n', ...
    ac_dx_mm, num2str(ac_dims));
fprintf('[thermal_grid_setup] Thermal  grid: %.2f mm, [%s] vox  (FOV crop: vox %s–%s)\n', ...
    th_dx_mm, num2str(th_dims), num2str(fov_start_ac), num2str(fov_end_ac));

% -------------------------------------------------------------------------
% Resample medium fields (linear) and masks (nearest)
% -------------------------------------------------------------------------
kwave_medium_th = kwave_medium;
field_names = fieldnames(kwave_medium);
for fi = 1:numel(field_names)
    fname = field_names{fi};
    vol   = kwave_medium.(fname);
    if ~isnumeric(vol) || ~isequal(size(vol), double(ac_dims))
        continue   % skip scalars and non-grid fields
    end
    % 1. Crop to thermal FOV
    if ndim == 3
        vol_crop = vol(fov_start_ac(1):fov_end_ac(1), ...
                       fov_start_ac(2):fov_end_ac(2), ...
                       fov_start_ac(3):fov_end_ac(3));
    else
        vol_crop = vol(fov_start_ac(1):fov_end_ac(1), ...
                       fov_start_ac(2):fov_end_ac(2));
    end
    % 2. Resample to thermal dims
    if th_dx_mm ~= ac_dx_mm
        kwave_medium_th.(fname) = imresize3(single(vol_crop), th_dims, 'linear');
    else
        kwave_medium_th.(fname) = vol_crop;
    end
end

% Resample masks with nearest-neighbour to preserve integer labels
if ndim == 3
    masks_crop = medium_masks(fov_start_ac(1):fov_end_ac(1), ...
                               fov_start_ac(2):fov_end_ac(2), ...
                               fov_start_ac(3):fov_end_ac(3));
else
    masks_crop = medium_masks(fov_start_ac(1):fov_end_ac(1), ...
                               fov_start_ac(2):fov_end_ac(2));
end
if th_dx_mm ~= ac_dx_mm
    medium_masks_th = imresize3(int16(masks_crop), th_dims, 'nearest');
else
    medium_masks_th = masks_crop;
end

% -------------------------------------------------------------------------
% Update transducer / focus positions in thermal grid voxel coordinates
% -------------------------------------------------------------------------
% Physical offset of the thermal FOV origin from the acoustic grid origin [mm]
fov_origin_mm = (fov_start_ac - 1) * ac_dx_mm;

% From acoustic grid voxels to thermal grid voxels:
%   th_vox = ((ac_vox - 1) * ac_dx - fov_origin_mm) / th_dx + 1
trans_pos_th  = round(((trans_pos  - 1) * ac_dx_mm - fov_origin_mm) / th_dx_mm + 1);
focus_pos_th  = round(((focus_pos  - 1) * ac_dx_mm - fov_origin_mm) / th_dx_mm + 1);

% Clamp to thermal grid (trans_pos may be outside when using cropped FOV)
trans_pos_th = max(1, min(trans_pos_th, th_dims));
focus_pos_th = max(1, min(focus_pos_th, th_dims));

% -------------------------------------------------------------------------
% Compute T1 → thermal grid affine (transf_th)
% -------------------------------------------------------------------------
% transf (planimg.transf) maps T1 voxels → acoustic grid voxels (row-vector convention).
% The acoustic→thermal mapping is:
%   th_vox_i = (ac_vox_i - 1) * scale - fov_origin_mm/th_dx_mm + 1
%             = ac_vox_i * scale + (1 - scale - fov_origin_mm/th_dx_mm)
%
% In 4×4 homogeneous (row-vector) form:
%   [th 1] = [ac 1] * M_ac2th
M_ac2th = diag([scale, scale, scale, 1]);
offset = 1 - scale - fov_origin_mm / th_dx_mm;   % [1×ndim]
if ndim == 3
    M_ac2th(4, 1:3) = offset;
else
    M_ac2th(3, 1:2) = offset(1:2);
end
transf_th = planimg.transf * M_ac2th;

% -------------------------------------------------------------------------
% Clone parameters and patch for thermal grid
% -------------------------------------------------------------------------
parameters_th = parameters;
parameters_th.grid.dims = th_dims;
parameters_th.transducer(1).trans_pos = trans_pos_th;
if isfield(parameters_th.transducer(1), 'focus_pos')
    parameters_th.transducer(1).focus_pos = focus_pos_th;
end

end
