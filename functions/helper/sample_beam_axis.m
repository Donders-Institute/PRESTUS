function [dist_mm, values] = sample_beam_axis(field, trans_pos, focus_pos, resolution_mm)
% SAMPLE_BEAM_AXIS  Sample a field along the beam axis from the transducer bowl to the grid edge
%
% Linearly interpolates a 2D or 3D field along the line starting at
% trans_pos and passing through focus_pos, in half-voxel steps, up to the
% last point inside the grid. Distances are measured from trans_pos (the
% bowl apex), so the result is independent of the grid frame (e.g.,
% axisymmetric vs. expanded grids, or an acoustic FOV crop).
%
% Use as:
%   [dist_mm, values] = sample_beam_axis(field, trans_pos, focus_pos, resolution_mm)
%
% Input:
%   field         - 2D or 3D numeric array on the simulation grid
%   trans_pos     - transducer (bowl apex) position [voxels]
%   focus_pos     - intended focus position [voxels]
%   resolution_mm - grid resolution [mm]
%
% Output:
%   dist_mm - [1 x N] distance from trans_pos along the beam axis [mm]
%   values  - [1 x N] interpolated field values (NaN outside the grid)
%
% See also: WATER_BASELINE, PLOT_AXIAL_PROFILE

arguments
    field         {mustBeNumeric}
    trans_pos     {mustBeNumeric}
    focus_pos     {mustBeNumeric}
    resolution_mm (1,1) {mustBeNumeric, mustBePositive}
end

    n_dims    = ndims(field);
    trans_pos = double(trans_pos(1:n_dims));
    focus_pos = double(focus_pos(1:n_dims));

    beam_vec = focus_pos - trans_pos;
    if norm(beam_vec) == 0
        error('sample_beam_axis: trans_pos equals focus_pos; beam axis undefined.');
    end
    beam_dir = beam_vec / norm(beam_vec);

    step_vox = 0.5;
    d_vox = 0:step_vox:norm(size(field));
    pts = trans_pos(:) + beam_dir(:) * d_vox;   % [n_dims x N]

    grid_vecs = arrayfun(@(d) 1:size(field, d), 1:n_dims, 'UniformOutput', false);
    pts_cell  = num2cell(pts, 2);
    values = interpn(grid_vecs{:}, double(gather(field)), pts_cell{:}, 'linear', NaN);
    values = values(:)';

    inside  = find(~isnan(values), 1, 'last');
    values  = values(1:inside);
    dist_mm = d_vox(1:inside) * resolution_mm;
end
