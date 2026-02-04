function [parameters] = grid_transducer_location(parameters, planimg)
% Position transducer(s) in the grid

    if contains(parameters.simulation_medium, {'layered'})
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

        if (~isfield(parameters.transducer, 'trans_pos') || isempty(parameters.transducer.trans_pos)) ...
                || (~isfield(parameters.transducer, 'focus_pos')|| isempty(parameters.transducer.focus_pos))
            disp('Either grid or focus position is not set, positioning them arbitrarily based on the focal distance')
            % note that the focus position matters only for the orientation of the transducer
        end
        % set transducer position in grid
        if ~isfield(parameters.transducer, 'trans_pos') || isempty(parameters.transducer.trans_pos)
            % transducer positioned arbitrarily (2D only)
            % y: first position beyond pml layer
            % x: halfway
            trans_pos = round(...
                [parameters.grid_dims(1:(parameters.n_sim_dims-1))/2, ...
                parameters.pml_size+1]);
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
            % this already accounts for PML size
            focus_pos = trans_pos;
            focus_pos(parameters.n_sim_dims) = ...
                round(focus_pos(parameters.n_sim_dims) + ...
                parameters.expected_focal_distance_mm/parameters.grid_step_mm);
        else
            focus_pos = parameters.transducer.focus_pos;
            % Adjust if the positions are transposed (2D only)
            % In 2D, we expect the second index as the axial dimension
            if parameters.n_sim_dims == 2 && focus_pos(1)>focus_pos(2)
                warning('Specified focus position appears transposed...adjusting');
                focus_pos = focus_pos';
            end
        end
        % Retain transducer and focus positions
        parameters.transducer(1).trans_pos = trans_pos;
        parameters.transducer(1).focus_pos = focus_pos;
    end
    
    % If a PML layer is used to absorb waves reaching the edge of the grid,
    % this will check if there is enough room for a PML layer between the
    % transducers and the edge of the grid
    for ti = 1:numel(parameters.transducer)
        tp = parameters.transducer(ti).trans_pos;
        fp = parameters.transducer(ti).focus_pos;
        assert(min(abs([repmat(0, 1, numel(parameters.grid_dims));parameters.grid_dims]-...
            tp ),[],'all') > parameters.pml_size, ...
            sprintf('The minimal distance between the transducer %i and the simulation grid boundary should be larger than the PML size. Adjust transducer position or the PML size', ti))
        assert(min(abs([repmat(0, 1, numel(parameters.grid_dims));parameters.grid_dims]-...
            fp ),[],'all') > parameters.pml_size, ...
            sprintf('The minimal distance between the focus position of transducer %i and the simulation grid boundary should be larger than the PML size. Adjust transducer position or the PML size', ti))
    end