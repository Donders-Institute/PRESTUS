function [dims, info] = water_grid_dims(parameters)
% WATER_GRID_DIMS  Free-water grid size from the transducer aperture and focal depth
%
% Resolves grid.default_dims = 'auto' for water simulations. The free-water
% grid needs to hold the transducer and the beam up to a little beyond the
% focus, as the acoustic FOV of a head simulation does
% (RESOLVE_ACOUSTIC_FOV_PARAMS):
%
%   lateral = aperture diameter + 2 * grid.water_margin_mm      (default 5 mm)
%   axial   = focal_distance_bowl + grid.z_beyond_focus_mm      (default 40 mm)
%
% The transducer sits at the near face of the grid (GRID_TRANSDUCER_LOCATION),
% so the axial length is measured from the bowl apex. The PML is added
% outside these dimensions by k-Wave. With grid.axisymmetric the result is
% [Nlateral, Naxial] (full lateral width; GRID_AXISYMMETRY keeps one half),
% else [Nlateral, Nlateral, Naxial].
%
% Use as:
%   dims = water_grid_dims(parameters)
%   [dims, info] = water_grid_dims(parameters)
%
% Input:
%   parameters - PRESTUS config with grid.resolution_mm and transducer(1)
%                (annular: annular.elem_od_mm; matrix: matrix.outer_diameter_mm;
%                focal_distance_ep or focal_distance_bowl)
%
% Output:
%   dims - [1x2] or [1x3] grid points
%   info - .aperture_mm, .lateral_mm, .axial_mm, .focal_distance_bowl_mm
%
% See also: GRID_TISSUE_SETUP, RESOLVE_ACOUSTIC_FOV_PARAMS, KWAVE_CALIBRATE_DRIVE

    res = parameters.grid.resolution_mm;
    margin_mm = get_num(parameters.grid, 'water_margin_mm', 5);
    beyond_mm = get_num(parameters.grid, 'z_beyond_focus_mm', 40);

    tr = parameters.transducer(1);
    aperture_mm = [];
    if isfield(tr, 'annular') && isfield(tr.annular, 'elem_od_mm') && ~isempty(tr.annular.elem_od_mm)
        aperture_mm = max(tr.annular.elem_od_mm);
    elseif isfield(tr, 'elem_od_mm') && ~isempty(tr.elem_od_mm)
        aperture_mm = max(tr.elem_od_mm);
    elseif isfield(tr, 'matrix') && isfield(tr.matrix, 'outer_diameter_mm') && ~isempty(tr.matrix.outer_diameter_mm)
        aperture_mm = tr.matrix.outer_diameter_mm;
    end
    if isempty(aperture_mm)
        error('water_grid_dims:aperture', ['grid.default_dims ''auto'' needs the transducer aperture ' ...
            '(annular.elem_od_mm or matrix.outer_diameter_mm). Set grid.default_dims explicitly instead.']);
    end

    p = focal_distance_calculation(parameters);
    fd_bowl = p.transducer(1).focal_distance_bowl;
    if isempty(fd_bowl) || ~isfinite(fd_bowl)
        error('water_grid_dims:focus', ['grid.default_dims ''auto'' needs transducer.focal_distance_ep ' ...
            '(or focal_distance_bowl). Set grid.default_dims explicitly instead.']);
    end

    lateral_mm = aperture_mm + 2 * margin_mm;
    axial_mm   = fd_bowl + beyond_mm;
    n_lat = 2 * ceil(lateral_mm / res / 2) + 1;      % odd: the beam axis is a grid line
    n_ax  = ceil(axial_mm / res) + 2;                % transducer at voxel 2
    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric
        dims = [n_lat, n_ax];
    else
        dims = [n_lat, n_lat, n_ax];
    end
    info = struct('aperture_mm', aperture_mm, 'lateral_mm', lateral_mm, 'axial_mm', axial_mm, ...
        'focal_distance_bowl_mm', fd_bowl);
end

function v = get_num(s, f, d)
    if isfield(s, f) && isnumeric(s.(f)) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
