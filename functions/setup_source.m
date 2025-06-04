function [source, source_labels, transducer_pars] = setup_source(parameters, kgrid, trans_pos, focus_pos)

% SETUP_SOURCE Creates a source for k-Wave simulations based on transducer parameters.
%
% This function generates a source structure for k-Wave simulations using the provided 
% transducer parameters, computational grid, and transducer/focus positions. It supports 
% both custom element geometries and kWaveArray-based transducer setups.
%
% Input:
%   parameters     - Struct containing simulation parameters (e.g., transducer settings).
%   kgrid          - Struct representing the k-Wave computational grid (e.g., `kWaveGrid`).
%   trans_pos      - [1x3] or [3x1] array specifying the transducer position in grid coordinates.
%   focus_pos      - [1x3] or [3x1] array specifying the geometric focus position in grid coordinates.
%
% Output:
%   source         - Struct containing the source mask and pressure signals for k-Wave simulations.
%   source_labels  - 3D matrix labeling each transducer element in the computational grid.
%   transducer_pars- Updated struct with additional fields (e.g., element diameters in grid points).

    % Validate input dimensions
    if ~isequal(size(focus_pos), size(trans_pos))
        error('Transducer and focus positions should be arrays of equal size');
    end

    if isequal(size(focus_pos), [parameters.n_sim_dims 1])
        focus_pos = focus_pos';
        trans_pos = trans_pos';
    elseif ~isequal(size(focus_pos), [1 parameters.n_sim_dims])
        error('Transducer and focus positions should have size [N 1] or [1 N], where N is the number of simulation dimensions.');
    end

    % Convert element diameters from mm to grid points
    transducer_pars = parameters.transducer;
    grid_step_mm = parameters.grid_step_mm;

    transducer_pars.Elements_OD = 2 * floor(transducer_pars.Elements_OD_mm / grid_step_mm / 2) + 1; % Outer diameter in grid points
    transducer_pars.Elements_ID = 2 * floor(transducer_pars.Elements_ID_mm / grid_step_mm / 2) + 1; % Inner diameter in grid points
    transducer_pars.Elements_ID(transducer_pars.Elements_ID_mm == 0) = 0;

    % Convert curvature radius from mm to grid points
    transducer_pars.radius_grid = round(transducer_pars.curv_radius_mm / grid_step_mm);

    %% Setup source signal
    
    cw_signal = createCWSignals(...
        kgrid.t_array, ...
        parameters.transducer.source_freq_hz, ...
        parameters.transducer.source_amp, ...
        parameters.transducer.source_phase_rad);

    source = struct();

    %% Custom element geometry (non-kWaveArray setup)
    if parameters.use_kWaveArray == 0
        grid_dims = parameters.grid_dims;
        transducer_mask = zeros(grid_dims);
        source_labels = zeros(grid_dims);

        % Create element bowls one by one
        for el_i = 1:transducer_pars.n_elements
            if parameters.n_sim_dims == 3
                bowl = makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, ...
                                transducer_pars.Elements_OD(el_i), focus_pos);
            else
                bowl = makeArc(grid_dims, trans_pos, transducer_pars.radius_grid, ...
                               transducer_pars.Elements_OD(el_i), focus_pos);
            end

            % Subtract inner bowl geometry if applicable
            if transducer_pars.Elements_ID(el_i) > 0
                if length(grid_dims) == 3
                    bowl = bowl - makeBowl(grid_dims, trans_pos, ...
                                           transducer_pars.radius_grid, ...
                                           transducer_pars.Elements_ID(el_i), focus_pos);
                else
                    bowl = bowl - makeArc(grid_dims, trans_pos, ...
                                          transducer_pars.radius_grid, ...
                                          transducer_pars.Elements_ID(el_i), focus_pos);
                end
            end

            % Define the binary source mask and label elements
            transducer_mask = transducer_mask + bowl;
            source_labels = source_labels + el_i * bowl;
        end

        % Assign pressure signals to each source point
        p_mask_source_p = reshape(source_labels(source_labels > 0), [], 1);
        source.p = zeros(length(p_mask_source_p), length(cw_signal));

        for ii = 1:length(p_mask_source_p)
            source.p(ii, :) = cw_signal(p_mask_source_p(ii), :);
        end

        source.p_mask = transducer_mask;

    %% kWaveArray-based setup (for advanced simulations)
    else
        disp('Setting up kWaveArray (might take a bit of time)');

        % Determine if axisymmetric mode should be enabled
        if numel(trans_pos) == 2 && ...
                isfield(parameters, 'axisymmetric') && ...
                parameters.axisymmetric == 1
            axisymmetric = true;
            disp("Using axisymmetric setup for 2D input...")
        else
            axisymmetric = false;
        end
        
        % Create kWaveArray object (always without axisymmetry)
        karray = kWaveArray('Axisymmetric', false, ...
                           'BLITolerance', 0.1, ...
                           'UpsamplingRate', 10, ...
                           'BLIType', 'sinc');
        
        % Set focus position and transducer position vectors in physical coordinates
        if numel(trans_pos) == 3
            % 3D annular array
            pos_vec = [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))];
            focus_vec = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))];
            
            karray.addAnnularArray(pos_vec, ...
                                   transducer_pars.curv_radius_mm * 1e-3, ...
                                   [transducer_pars.Elements_ID_mm; transducer_pars.Elements_OD_mm] * 1e-3, ...
                                   focus_vec);
        elseif numel(trans_pos) == 2 && axisymmetric == false

            % 2D arc-shaped element
            pos_vec = [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2))];
            focus_vec = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2))];
            
            karray.addArcElement(pos_vec, ...
                                 transducer_pars.curv_radius_mm * 1e-3, ...
                                 transducer_pars.Elements_OD_mm * 1e-3, ...
                                 focus_vec);
        end
        
        if axisymmetric
            
            % Approach: create a double sized grid and crop the source

            % make a new grid with odd number of pts. 
            kgrid_mirrored = kWaveGrid(kgrid.Nx, kgrid.dx, 2*kgrid.Ny - 1, kgrid.dy);
            karray_full = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 100);
            % position of the central point on the bowl surface
            x_offset = (trans_pos(1)-1).*(1/parameters.grid_step_mm); % account for the plm (mm offset)
            position_full = [kgrid.x_vec(1) + x_offset*kgrid.dx, 0+eps];
            % lateral grid point
%             focus_pos_half = [0, kgrid.y_vec(1)];
            focus_pos_full = [0, 0+eps]; % eps =  smallest step in MATLAB
            % add the arc element (2D bowl)
            karray_full.addArcElement(position_full, ...
                                 transducer_pars.curv_radius_mm * 1e-3, ...
                                 transducer_pars.Elements_OD_mm * 1e-3, ...
                                 focus_pos_full);
            % Get binary mask for the array elements on the grid
            source.p_mask = karray_full.getArrayBinaryMask(kgrid_mirrored);
            % Get the source weighting for each point in the mask
            grid_weights = karray_full.getElementGridWeights(kgrid_mirrored, 1);
            % Get distributed source signal weighted by element geometry
            % consider half of the y-axis space for axisymmetry now
            source.p_mask = source.p_mask(:,kgrid.Ny:end);
            grid_weights = grid_weights(:, kgrid.Ny:end);
            tmp = grid_weights(source.p_mask);
            source.p = bsxfun(@times, cw_signal, tmp(:));
            % binarize the source mask based on weight threshold
            source_labels = (source.p_mask.*grid_weights)>.25;
            
        else
            % --- Non-axisymmetric mode: compute weights and distribute signals manually ---
            
            % Initialize 4D array for element weights: [elements, Nx, Ny, Nz]
            grid_weights_4d = zeros(transducer_pars.n_elements, kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));
            
            for ind = 1:transducer_pars.n_elements
                fprintf('Computing weights for element %i...', ind);
                grid_weights_4d(ind,:,:,:) = karray.getElementGridWeights(kgrid, ind);
                fprintf(' done\n');
            end
            
            % Sum weights over elements and create binary mask of active source points
            binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;
            
            Nt = size(cw_signal, 2);
            mask_ind = find(binary_mask);
            num_source_points = sum(binary_mask(:));
            
            distributed_source_signal = zeros(num_source_points, Nt);
            source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));
            
            if canUseGPU
                cw_signal = gpuArray(cw_signal);
                distributed_source_signal = gpuArray(distributed_source_signal);
            end
            
            % Distribute signals weighted by element weights
            for ind = 1:transducer_pars.n_elements
                source_weights = squeeze(grid_weights_4d(ind,:,:,:));
                el_binary_mask = source_weights ~= 0;
                
                if canUseGPU
                    source_weights = gpuArray(source_weights);
                end
                
                element_mask_ind = find(el_binary_mask);
                local_ind = ismember(mask_ind, element_mask_ind);
                
                distributed_source_signal(local_ind,:) = ...
                    distributed_source_signal(local_ind,:) + ...
                    source_weights(element_mask_ind) * cw_signal(ind,:);
                
                source_labels = source_labels + ind * el_binary_mask;
            end
            
            % Assign outputs
            source.p_mask = binary_mask;
            source.p = distributed_source_signal;
        end


    end

end
