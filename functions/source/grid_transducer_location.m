function [parameters] = grid_transducer_location(parameters, planimg)
% GRID_TRANSDUCER_LOCATION  Map transducer and focus positions into the simulation grid
%
% For layered simulations, applies the affine transform stored in
% planimg.transf to map T1-space coordinates of every transducer to
% simulation-grid indices. For non-layered (water / phantom) simulations,
% positions are either taken directly from parameters or placed
% automatically: the transducer at the near face of the grid, the focus at
% the expected focal distance along the axial dimension. All positions are
% validated to lie within [1, grid.dims] (the PML is added outside by
% k-Wave so no internal clearance is required).
%
% Use as:
%   [parameters] = grid_transducer_location(parameters, planimg)
%
% Input:
%   parameters - PRESTUS config; must contain simulation.medium, grid.dims,
%                grid.resolution_mm [mm], transducer(i).trans_pos, transducer(i).focus_pos
%   planimg    - planning image info from GRID_TISSUE_SETUP; must contain
%                planimg.transf (affine matrix) for layered simulations;
%                may be empty for water/phantom
%
% Output:
%   parameters - updated: transducer(i).trans_pos and transducer(i).focus_pos
%                set to integer grid indices
%
% See also: GRID_TISSUE_SETUP, GRID_AXISYMMETRY, FOCAL_DISTANCE_CALCULATION

arguments
    parameters (1,1) struct
    planimg    (1,1) struct
end

    if contains(parameters.simulation.medium, {'layered'})
        % map all transducers from T1 grid to sim grid using the same transform
        for ti = 1:numel(parameters.transducer)
            tr = parameters.transducer(ti);

            % require T1-grid positions for each transducer
            assert(isfield(tr, 'trans_pos') && isfield(tr, 'focus_pos') && ...
                   ~isempty(tr.trans_pos) && ~isempty(tr.focus_pos), ...
                   'trans_pos and focus_pos must be defined for each transducer');

            pts_t1 = [tr.trans_pos(:).'; tr.focus_pos(:).'];  % 2×3
            pts_sim = round(tformfwd(pts_t1, maketform('affine', planimg.transf))); % 2×3

            parameters.transducer(ti).trans_pos  = pts_sim(1,:);
            parameters.transducer(ti).focus_pos  = pts_sim(2,:);
        end
    else
        % non-layered media currently only support a single transducer
        % only the first specified transducer will be modeled
        % this is a known limitation

        if numel(parameters.transducer)>1
            warning("Non-skull and/or non-layered simulations currently only support a single transducer. Only the first specified transducer will be retained...");
            parameters.transducer = parameters.transducer(1);
        end

        if strcmp(parameters.transducer.type, 'annular')
            % for water medium remove potential position specifications
            % the grid has an arbitrary size that does not necessarily map onto the planning image
            if strcmp(parameters.simulation.medium, 'water')
                parameters.transducer.trans_pos = [];
                parameters.transducer.focus_pos = [];
            end
        end

        if (~isfield(parameters.transducer, 'trans_pos') || isempty(parameters.transducer.trans_pos)) ...
                || (~isfield(parameters.transducer, 'focus_pos')|| isempty(parameters.transducer.focus_pos))
            disp('Either grid or focus position is not set, positioning them arbitrarily based on the focal distance')
            % note that the focus position matters only for the orientation of the transducer
        end
        % set transducer position in grid
        if ~isfield(parameters.transducer, 'trans_pos') || isempty(parameters.transducer.trans_pos)
            % transducer positioned arbitrarily (2D only)
            % axial: first interior voxel (PML is added outside the grid by k-Wave)
            % lateral: centred
            trans_pos = round(...
                [parameters.grid.dims(1:(numel(parameters.grid.dims)-1))/2, ...
                2]);
        else
            trans_pos = parameters.transducer.trans_pos;
            % Adjust if the positions are transposed
            if size(trans_pos,1)>size(trans_pos, 2)
                warning('Specified transducer position appears transposed...adjusting');
                trans_pos = trans_pos';
            end
        end
        % set focus position in grid
        if ~isfield(parameters.transducer, 'focus_pos') || isempty(parameters.transducer.focus_pos)
            % no focus point specified
            % position focus at expected distance from transducer
            % index dimension depends on 2D/3D
            parameters = focal_distance_calculation(parameters);
            focus_pos = trans_pos;
            focus_pos(numel(parameters.grid.dims)) = ...
                round(focus_pos(numel(parameters.grid.dims)) + ...
                parameters.transducer(1).focal_distance_bowl/parameters.grid.resolution_mm);
        else
            focus_pos = parameters.transducer.focus_pos;
            % Adjust if the positions are transposed (2D only)
            % In 2D, we expect the second index as the axial dimension
            if numel(parameters.grid.dims) == 2 && focus_pos(1)>focus_pos(2)
                warning('Specified focus position appears transposed...adjusting');
                focus_pos = focus_pos';
            end
        end
        % Retain transducer and focus positions
        parameters.transducer(1).trans_pos = trans_pos;
        parameters.transducer(1).focus_pos = focus_pos;
    end
    
    % Verify that transducer and focus positions lie within the interior grid.
    % The PML is added outside the grid by k-Wave (PMLInside=false), so positions
    % only need to be in [1, dims] — no in-grid PML clearance is required.
    for ti = 1:numel(parameters.transducer)
        tp = parameters.transducer(ti).trans_pos;
        fp = parameters.transducer(ti).focus_pos;
        assert(all(tp >= 1) && all(tp <= parameters.grid.dims), ...
            sprintf('Transducer %i position %s is outside the grid [1, %s]. Adjust trans_pos.', ...
            ti, mat2str(tp), mat2str(parameters.grid.dims)))
        assert(all(fp >= 1) && all(fp <= parameters.grid.dims), ...
            sprintf('Focus position of transducer %i (%s) is outside the grid [1, %s]. Adjust focus_pos.', ...
            ti, mat2str(fp), mat2str(parameters.grid.dims)))
    end