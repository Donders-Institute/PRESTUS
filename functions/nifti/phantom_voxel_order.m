function [data, transform] = phantom_voxel_order(data, parameters, planimg)
% PHANTOM_VOXEL_ORDER  Return a non-layered output volume to the voxel order of its input
%
% Undoes what the grid stage did to the input before simulating, so that an
% output NIfTI overlays the phantom it was simulated on:
%   - axisymmetric runs: a 3D expansion is cut through the beam axis
%     (transducer(1).trans_pos_bilateral(2)); the [Nlateral x Naxial] plane is
%     transposed back when GRID_AXISYMMETRY found the input axial-first
%     (grid.axisym_transposed = false);
%   - a 2D phantom transposed to grid.default_dims (planimg.phantom_transposed)
%     is transposed back;
%   - a plane cut from a 3D phantom (planimg.phantom_slice) gets the phantom's
%     affine for that plane.
%
% Use as:
%   [data, transform] = phantom_voxel_order(data, parameters, planimg)
%
% Input:
%   data       - simulation-grid output volume (2D or 3D)
%   parameters - PRESTUS config (grid.axisymmetric, grid.axisym_bilateral_dims,
%                grid.axisym_transposed, transducer(1).trans_pos_bilateral)
%   planimg    - planning image struct (phantom_header, phantom_transposed,
%                phantom_slice; all optional)
%
% Output:
%   data      - volume in the phantom file's voxel order
%   transform - affine3d for the output NIfTI header (the phantom's, or the
%               phantom's restricted to the cut plane); [] without a phantom
%               header or when the phantom has no transform (TransformName 'None')
%
% See also: NIFTI_TO_T1W, GRID_TISSUE_SETUP, GRID_AXISYMMETRY

arguments
    data
    parameters (1,1) struct
    planimg    (1,1) struct
end

    transform = [];
    grid = parameters.grid;
    axisym = isfield(grid, 'axisymmetric') && grid.axisymmetric == 1;
    has_snapshot = isfield(grid, 'axisym_bilateral_dims') && ...
        isfield(parameters.transducer(1), 'trans_pos_bilateral') && ...
        ~isempty(parameters.transducer(1).trans_pos_bilateral);

    in_plane = false;   % data is the [Nlateral x Naxial] bilateral plane
    if axisym && ndims(data) == 3
        ax = round(size(data, 1) / 2);
        if has_snapshot && size(data, 1) == grid.axisym_bilateral_dims(2) && ...
                size(data, 2) == grid.axisym_bilateral_dims(2)
            ax = parameters.transducer(1).trans_pos_bilateral(2);   % beam axis
        end
        data = reshape(data(ax, :, :), size(data, 2), size(data, 3));
        in_plane = has_snapshot && isequal(size(data), fliplr(grid.axisym_bilateral_dims));
    elseif axisym && has_snapshot && ismatrix(data)
        in_plane = isequal(size(data), fliplr(grid.axisym_bilateral_dims));
    end
    % Axial-first input: GRID_AXISYMMETRY did not transpose it
    if in_plane && isfield(grid, 'axisym_transposed') && ~grid.axisym_transposed
        data = data.';
    end

    if ~isfield(planimg, 'phantom_header') || isempty(planimg.phantom_header)
        return
    end
    % A phantom without a spatial transform gives none to reuse
    hdr = planimg.phantom_header;
    if ~(isfield(hdr, 'TransformName') && strcmpi(hdr.TransformName, 'None'))
        transform = hdr.Transform;
    end

    if ismatrix(data) && isfield(planimg, 'phantom_transposed') && ...
            ~isempty(planimg.phantom_transposed) && planimg.phantom_transposed
        data = data.';
    end

    if ismatrix(data) && isfield(planimg, 'phantom_slice') && ~isempty(planimg.phantom_slice) && ...
            ~isempty(transform)
        % Plane at index k of dim d of the 3D phantom: its two in-plane voxel
        % axes, the cut axis as the (singleton) third, origin moved to the plane
        % (row-vector convention, 0-based voxels: [ras 1] = [i j k 1] * T)
        T = affine_matrix(transform);
        d = planimg.phantom_slice.dim;
        k = planimg.phantom_slice.index;
        keep = setdiff(1:3, d);
        T2 = eye(4);
        T2(1, :) = T(keep(1), :);
        T2(2, :) = T(keep(2), :);
        T2(3, :) = T(d, :);
        T2(4, :) = T(4, :) + (k - 1) * T(d, :);
        T2(1:3, 4) = 0; T2(4, 4) = 1;
        transform = affine3d(T2);
    end
end

function T = affine_matrix(tf)
% 4x4 row-vector matrix of an affine3d (T) or affinetform3d (A, column convention)
    if isprop(tf, 'T')
        T = tf.T;
    else
        T = tf.A';
    end
end
