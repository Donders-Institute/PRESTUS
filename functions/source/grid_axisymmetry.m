function [parameters, segmentation, bone_mask, pseudoCT, medium_masks] = ...
    grid_axisymmetry(parameters, segmentation, bone_mask, pseudoCT, medium_masks)
% GRID_AXISYMMETRY  Adapt grid and tissue masks for k-Wave axisymmetric simulation
%
% When parameters.grid.axisymmetric == 1 and a 2D focus position is specified,
% the full bilateral grid is halved along the radial axis so that column 1
% of the returned arrays corresponds to r = 0 (the axis of symmetry), as
% required by kspaceFirstOrderAS. Transducer and focus positions are updated
% accordingly. Has no effect when axisymmetric mode is not requested.
%
% --- Axisymmetric coordinate pipeline ---
%
%  Step 1 — grid_axisymmetry (this function)
%    Input:  bilateral grid [Nlateral x Naxial] or [Naxial x Nlateral];
%            transducer at the axis column (trans_pos(lateral) = ax_center).
%    Action: if dims(2) > dims(1), transpose so the layout is [Naxial x Nlateral]
%            (axial along rows, lateral along columns).
%            Save the bilateral snapshot:
%              parameters.grid.axisym_bilateral_dims     = [Naxial, Nlateral]
%              parameters.transducer(1).trans_pos_bilateral = [axial, lateral]
%              parameters.transducer(1).focus_pos_bilateral = [axial, lateral]
%            Halve: keep columns ax_center:end → [Naxial x Nr], Nr = Nlateral-ax_center+1.
%            Set trans_pos = [axial, 1], focus_pos = [axial_focus, 1].
%    Output: half-grid [Naxial x Nr], grid.dims = [Naxial, Nr].
%
%  Step 2 — kgrid, segmentation, medium_masks are now Nr columns wide.
%    kspaceFirstOrderAS normalises y_vec = kgrid.y_vec - kgrid.y_vec(1),
%    so column 1 is always r = 0 regardless of the absolute y coordinates.
%
%  Step 3 — source_create (kWaveArray mirrored-grid approach for annular)
%    Builds kgrid_mirrored with 2*Nr-1 columns (odd; exact centre at y=0).
%    Computes element weights on the full bilateral mirrored grid.
%    Crops columns Nr:end → [Naxial x Nr] source.p_mask where col 1 = r=0.
%
%  Step 4 — kspaceFirstOrderAS simulation on [Naxial x Nr] grid.
%    PML chosen with WSWA symmetry (getOptimalPMLSize(..., 'WSWA')).
%    sensor_data.p_max_all: [Naxial x Nr] half-plane result.
%
%  Step 5 — convert_axisymmetric_to_2d  OR  convert_axisymmetric_to_3d
%    Restores the original bilateral coordinate space from the snapshot:
%      2D: mirrors half-plane → [Naxial x Nlateral], then transposes to
%          [Nlateral x Naxial]; sets grid.dims = [Nlateral, Naxial].
%      3D: radially expands → [Nlateral x Nlateral x Naxial];
%          sets grid.dims = [Nlateral, Nlateral, Naxial].
%    trans_pos and focus_pos are restored from the bilateral snapshot.
%    All downstream code (thermal, reporting, NIfTI export) sees the
%    original grid as if axisymmetric mode was never used.
%
% Requirement: grid.mode must be 'transducer_axis' (default); the beam axis
% must be aligned with the grid z-axis.  'ras_plus' mode is incompatible.
%
% Use as:
%   [parameters, segmentation, bone_mask, pseudoCT, medium_masks] = ...
%       grid_axisymmetry(parameters, segmentation, bone_mask, pseudoCT, medium_masks)
%
% Input:
%   parameters   - (struct) PRESTUS config; must contain grid.dims [1x2],
%                    grid.axisymmetric, and transducer(1).trans_pos / focus_pos
%   segmentation - [Nx x Ny] numeric, tissue label map
%   bone_mask    - [Nx x Ny] logical, binary skull mask
%   pseudoCT     - [Nx x Ny] numeric, Hounsfield-unit skull image ([] when unused)
%   medium_masks - [Nx x Ny] numeric, medium layer label map
%
% Output:
%   parameters   - updated: grid.dims = [Naxial, Nr] where Nr = Nlateral-ax_center+1;
%                    trans_pos = [axial, 1]; focus_pos = [axial_focus, 1];
%                    axisym_bilateral_dims, trans_pos_bilateral, focus_pos_bilateral stored
%   segmentation - [Naxial x Nr] numeric, halved along radial axis
%   bone_mask    - [Naxial x Nr] logical, halved along radial axis
%   pseudoCT     - [Naxial x Nr] numeric, halved along radial axis ([] when unused)
%   medium_masks - [Naxial x Nr] numeric, halved along radial axis
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, CONVERT_AXISYMMETRIC_TO_2D, SOURCE_CREATE

arguments
    parameters   (1,1) struct
    segmentation {mustBeNumericOrLogical}
    bone_mask    {mustBeNumericOrLogical}
    pseudoCT     {mustBeNumericOrLogical}
    medium_masks {mustBeNumericOrLogical}
end
    if isfield(parameters.grid, 'mode') && strcmp(parameters.grid.mode, 'ras_plus') && ...
            isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        error('PRESTUS:gridMode:incompatible', ...
            ['grid.mode=''ras_plus'' is incompatible with grid.axisymmetric=1. ' ...
             'Axisymmetric mode requires the beam axis to be aligned with the grid z-axis, ' ...
             'which is only guaranteed in the default ''transducer_axis'' mode.']);
    end

    if numel(parameters.transducer(1).focus_pos) == 2 && ...
            isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        if numel(parameters.transducer) > 1
            error('Axisymmetric simulations with multiple transducers are not supported (only a single transducer is allowed when axisymmetric == 1)');
        end
        trans_pos = parameters.transducer(1).trans_pos;
        focus_pos = parameters.transducer(1).focus_pos;
        % ensure that radial(y) dim is shorter than axial (x) dim
        if parameters.grid.dims(2) > parameters.grid.dims(1)
            parameters.grid.dims = fliplr(parameters.grid.dims);
            trans_pos    = fliplr(trans_pos);
            focus_pos    = fliplr(focus_pos);
            segmentation = segmentation';
            bone_mask    = bone_mask';
            if ~isempty(pseudoCT); pseudoCT = pseudoCT'; end
            medium_masks = medium_masks';
        end
        % Snapshot the bilateral state (post-flip, pre-halve) so that the
        % expansion functions can restore the original grid dimensions and
        % positions.  The axisymmetric half-grid is an internal detail;
        % all downstream outputs should be in the original coordinate space.
        parameters.grid.axisym_bilateral_dims     = parameters.grid.dims;          % [Naxial, Nlateral]
        parameters.transducer(1).trans_pos_bilateral = trans_pos;                  % [axial, lateral]
        parameters.transducer(1).focus_pos_bilateral = focus_pos;
        % halve the grid along the radial axis so that column 1 = r = 0.
        % kspaceFirstOrderAS always normalises y_vec to start at 0
        % (y_vec = kgrid.y_vec - kgrid.y_vec(1)), so whatever ends up in
        % column 1 is treated as the axis.  The transducer is on the axis,
        % so the correct split point is trans_pos(2): include that column
        % as column 1 of the half-grid and discard everything to its left.
        % Using floor(dims(2)/2)+1 instead would silently drop the axis
        % column and shift every phantom feature by one voxel radially.
        center_col = trans_pos(2);
        segmentation = segmentation(:, center_col:end);
        bone_mask    = bone_mask(:, center_col:end);
        if ~isempty(pseudoCT); pseudoCT = pseudoCT(:, center_col:end); end
        medium_masks = medium_masks(:, center_col:end);
        parameters.grid.dims(2) = size(segmentation, 2);
        % set transducer and focus position to the radial midline
        trans_pos(2) = 1;
        focus_pos(2) = 1;
        % Retain transducer and focus positions after grid manipulations
        parameters.transducer(1).trans_pos = trans_pos;
        parameters.transducer(1).focus_pos = focus_pos;
    end