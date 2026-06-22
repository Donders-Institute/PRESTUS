function [kgrid_th, kwave_medium_th, medium_masks_th, parameters_th, transf_th, sensor_th, source_th] = ...
    resolve_thermal_grid(parameters, has_external_acoustic, has_acoustic_fov, ...
                         use_babelbrain_acoustic, planimg, kwave_medium, medium_masks, ...
                         trans_pos, focus_pos, kgrid, sensor, source)
% RESOLVE_THERMAL_GRID  Set up the thermal simulation grid and placeholders.
%
% Decides whether to call thermal_grid_setup (independent thermal grid) or to
% reuse the acoustic grid objects directly, then builds sensor_th and source_th
% placeholders sized to whichever grid is used.
%
% An independent thermal grid is used when any of:
%   - An external acoustic NIfTI was provided (no acoustic kgrid exists)
%   - An acoustic FOV crop was active (kgrid must be rebuilt for correct domain)
%   - BabelBrain acoustic backend (never produces kgrid/source/sensor)
%   - A coarser thermal resolution is requested (grid.thermal_resolution_mm)
%   - An independent thermal FOV is specified (grid.thermal_fov_mm)
%
% Use as:
%   [kgrid_th, kwave_medium_th, medium_masks_th, parameters_th, transf_th, ...
%    sensor_th, source_th] = ...
%       resolve_thermal_grid(parameters, has_external_acoustic, has_acoustic_fov, ...
%                            use_babelbrain_acoustic, planimg, kwave_medium, medium_masks, ...
%                            trans_pos, focus_pos, kgrid, sensor, source)

use_thermal_grid = has_external_acoustic || has_acoustic_fov || use_babelbrain_acoustic || ...
    (isfield(parameters.grid, 'thermal_resolution_mm') && ...
     ~isempty(parameters.grid.thermal_resolution_mm) && ...
     parameters.grid.thermal_resolution_mm ~= parameters.grid.resolution_mm) || ...
    (isfield(parameters.grid, 'thermal_fov_mm') && ...
     ~isempty(parameters.grid.thermal_fov_mm));

if use_thermal_grid
    [kgrid_th, kwave_medium_th, medium_masks_th, parameters_th, transf_th] = ...
        thermal_grid_setup(parameters, planimg, kwave_medium, medium_masks, trans_pos, focus_pos);
    sensor_th      = struct();
    sensor_th.mask = ones(parameters_th.grid.dims);
    source_th      = struct();
    source_th.Q    = zeros(parameters_th.grid.dims);
else
    kgrid_th        = kgrid;
    kwave_medium_th = kwave_medium;
    medium_masks_th = medium_masks;
    parameters_th   = parameters;
    assert(isfield(planimg, 'transf') && ~isempty(planimg.transf), ...
        ['prestus_pipeline:missingPlanimg — planimg.transf is required for thermal ' ...
         'NIfTI export but was not found. ' ...
         'Re-run with run_grid_setup=1 to regenerate the grid cache and planimg sidecar.']);
    transf_th = planimg.transf;
    if ~isempty(sensor) && isstruct(sensor) && isfield(sensor, 'mask')
        sensor_th = sensor;
        source_th = source;
    else
        sensor_th      = struct();
        sensor_th.mask = ones(parameters_th.grid.dims);
        source_th      = struct();
        source_th.Q    = zeros(parameters_th.grid.dims);
    end
end

end
