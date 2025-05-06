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
        
        % Create empty kWaveArray object
        karray = kWaveArray('BLITolerance', 0.1, 'UpsamplingRate', 10, 'BLIType', 'sinc');

        % Add annular array elements to the kWaveArray object
        karray.addAnnularArray([kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))], ...
                               transducer_pars.curv_radius_mm * 1e-3, ...
                               [transducer_pars.Elements_ID_mm; transducer_pars.Elements_OD_mm] * 1e-3, ...
                               [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]);

        % Compute weights for each element in the computational grid
        grid_weights_4d = zeros(transducer_pars.n_elements, kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
        
        for ind = 1:transducer_pars.n_elements
            fprintf('Computing weights for element %i...', ind);
            grid_weights_4d(ind,:,:,:) = karray.getElementGridWeights(kgrid, ind);                
            fprintf(' done\n');
        end

        % Generate binary mask and distributed source signal
        binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;
        Nt = size(cw_signal, 2);

        mask_ind = find(binary_mask);
        num_source_points = sum(binary_mask(:));

        distributed_source_signal = zeros(num_source_points, Nt);
        
        source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
        
        if canUseGPU
            cw_signal = gpuArray(cw_signal);
            distributed_source_signal = gpuArray(distributed_source_signal);
        end

        % Assign signals to elements based on weights
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

        % Assign binary mask and distributed signals to the source structure
        source.p_mask = binary_mask;
        source.p = distributed_source_signal;
    end

end
