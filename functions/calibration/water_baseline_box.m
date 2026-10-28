function [kgrid_b, source_b, sensor_b, offset] = water_baseline_box(parameters, kgrid, source, sensor)
% WATER_BASELINE_BOX  Grid, source and sensor of the free-water baseline cropped to the beam
%
% The free-water baseline (WATER_BASELINE) only needs the source and the
% beam: in homogeneous water nothing outside reflects back. This crops the
% grid to the bounding box of the source points and of each transducer's
% beam axis from the bowl apex to 40 mm beyond the deeper of the focus and
% the geometric focus, widened by 10 mm, rounded up to a multiple of 8
% voxels per axis (within the grid). The source signals are kept as they
% are: k-Wave orders source points by linear index, and that order is the
% same inside any sub-box. Same time axis, so the steady state is reached.
%
% Returns the inputs unchanged (offset [0 0 0]) when the box would not save
% at least 20 % of the voxels, for axisymmetric grids, and for sources other
% than an additive pressure source on a mask.
%
% Use as:
%   [kgrid_b, source_b, sensor_b, offset] = water_baseline_box(parameters, kgrid, source, sensor)
%
% Output:
%   offset - [1 x n_dims] voxel offset of the box in the grid (box voxel v
%            is grid voxel v + offset)
%
% See also: WATER_BASELINE

    nd = kgrid.dim;
    dims = [kgrid.Nx, kgrid.Ny, kgrid.Nz]; dims = dims(1:nd);
    kgrid_b = kgrid; source_b = source; sensor_b = sensor; offset = zeros(1, nd);
    axisym = isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1;
    other_source = any(isfield(source, {'u_mask', 'ux', 'uy', 'uz', 'p0', 's_mask'}));
    if axisym || ~isfield(source, 'p_mask') || other_source || nd < 2
        return
    end

    res = parameters.grid.resolution_mm;
    margin = ceil(10 / res);
    sub = cell(1, nd);
    [sub{:}] = ind2sub(dims, find(source.p_mask));
    pts = cell2mat(sub);
    for ti = 1:numel(parameters.transducer)
        tr = parameters.transducer(ti);
        if ~isfield(tr, 'trans_pos') || numel(tr.trans_pos) ~= nd || ~isfield(tr, 'focus_pos') || numel(tr.focus_pos) ~= nd
            return                                      % positions unknown: keep the full grid
        end
        t = double(tr.trans_pos(:)'); f = double(tr.focus_pos(:)');
        u = (f - t) / max(norm(f - t), eps);
        reach = norm(f - t);
        if isfield(tr, char(tr.type)) && isfield(tr.(char(tr.type)), 'curv_radius_mm') && ...
                isfinite(tr.(char(tr.type)).curv_radius_mm)
            reach = max(reach, tr.(char(tr.type)).curv_radius_mm / res);
        end
        s = (0:1:ceil(reach + 40 / res))';
        pts = [pts; t + s .* u]; %#ok<AGROW>
    end
    lo = max(1, floor(min(pts, [], 1)) - margin);
    hi = min(dims, ceil(max(pts, [], 1)) + margin);
    n  = hi - lo + 1;
    n  = min(dims, ceil(n / 8) * 8);                    % FFT-friendlier, within the grid
    hi = min(dims, lo + n - 1); lo = max(1, hi - n + 1);
    if prod(hi - lo + 1) > 0.8 * prod(dims)
        return
    end

    idx = arrayfun(@(d) lo(d):hi(d), 1:nd, 'UniformOutput', false);
    nb = hi - lo + 1;
    dx = kgrid.dx;
    if nd == 3
        kgrid_b = kWaveGrid(nb(1), dx, nb(2), dx, nb(3), dx);
    else
        kgrid_b = kWaveGrid(nb(1), dx, nb(2), dx);
    end
    kgrid_b.setTime(kgrid.Nt, kgrid.dt);
    source_b.p_mask = source.p_mask(idx{:});
    if nnz(source_b.p_mask) ~= nnz(source.p_mask)
        error('PRESTUS:waterBaselineBox:source', 'The baseline box cuts the source (%d of %d points).', ...
            nnz(source_b.p_mask), nnz(source.p_mask));
    end
    if isfield(sensor, 'mask') && isequal(size(sensor.mask), dims)
        sensor_b.mask = sensor.mask(idx{:});
        if ~any(sensor_b.mask(:)), sensor_b.mask = true(nb); end
    end
    offset = lo - 1;
    fprintf('[water_baseline] beam box [%s] of grid [%s] (%.0f%% of the voxels)\n', ...
        num2str(nb), num2str(dims), 100 * prod(nb) / prod(dims));
end
