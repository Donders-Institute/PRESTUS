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
    if strcmp(parameters.transducer.element_shape, 'annular') || ~strcmp(parameters.transducer.element_shape.matrix.matrix_shape, 'dolphin')
        transducer_pars.Elements_OD = 2 * floor(transducer_pars.Elements_OD_mm / grid_step_mm / 2) + 1; % Outer diameter in grid points
        transducer_pars.Elements_ID = 2 * floor(transducer_pars.Elements_ID_mm / grid_step_mm / 2) + 1; % Inner diameter in grid points
        transducer_pars.Elements_ID(transducer_pars.Elements_ID_mm == 0) = 0;
    
        % Convert curvature radius from mm to grid points
        transducer_pars.radius_grid = round(transducer_pars.curv_radius_mm / grid_step_mm);
    end

    %% Custom element geometry (non-kWaveArray setup)
    if parameters.use_kWaveArray == 0
        %% Setup source signal
        cw_signal = createCWSignals(...
        kgrid.t_array, ...
        transducer_pars.source_freq_hz, ...
        transducer_pars.source_amp, ...
        transducer_pars.source_phase_rad);

        source = struct();

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

        if strcmp(parameters.transducer.element_shape.type, 'annular')
            % Add annular array elements to the kWaveArray object
            karray.addAnnularArray([kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))], ...
                                   transducer_pars.curv_radius_mm * 1e-3, ...
                                   [transducer_pars.Elements_ID_mm; transducer_pars.Elements_OD_mm] * 1e-3, ...
                                   [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]);

        %% matrix designs
        else
            matrix_pars = parameters.transducer.element_shape.matrix;
            if strcmp(matrix_pars.matrix_shape.type, 'define_here')
                is_rect_grid = false;
                if strcmp(matrix_pars.matrix_shape.define_here.grid_shape.type, 'rect')
                    is_rect_grid = true;
                end
                
                is_curved = matrix_pars.matrix_shape.define_here.is_curved;
                is_clover_module_cut = matrix_pars.matrix_shape.define_here.clover_module_cut;
                [transducer_pars, karray] = setup_circular_matrix_transducer(kgrid, trans_pos, transducer_pars, focus_pos, parameters, karray, is_rect_grid, is_curved, is_clover_module_cut);

            elseif strcmp(matrix_pars.matrix_shape.type, 'stanford')
                [transducer_pars, karray] = setup_stanford_transducer(kgrid, trans_pos, transducer_pars, focus_pos, parameters, karray);

            elseif strcmp(matrix_pars.matrix_shape.type, 'dolphin')
                [transducer_pars, karray] = setup_dolphin_transducer_design(transducer_pars, kgrid, trans_pos, focus_pos, karray, parameters.medium.water.sound_speed, parameters.n_sim_dims);
            else
                error('ERROR: %s is an unknown transducer design.', matrix_pars.matrix_shape.type);
            end
            
            % Update Elements OD and ID in the case part of elements have
            % been removed
            transducer_pars.Elements_OD = 2 * floor(transducer_pars.Elements_OD_mm / grid_step_mm / 2) + 1; % Outer diameter in grid points
            transducer_pars.Elements_ID = 2 * floor(transducer_pars.Elements_ID_mm / grid_step_mm / 2) + 1; % Inner diameter in grid points
            transducer_pars.Elements_ID(transducer_pars.Elements_ID_mm == 0) = 0;
        end
 
        % Define signal
        %% Setup source signal
        cw_signal = createCWSignals(...
        kgrid.t_array, ...
        transducer_pars.source_freq_hz, ...
        transducer_pars.source_amp, ...
        transducer_pars.source_phase_rad);

        source = struct();

        binary_mask = false(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
        source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
        Nt = size(cw_signal, 2);

        % Pre-allocate temporary arrays
        if canUseGPU
            cw_signal = gpuArray(cw_signal);
        end

        for ind = 1:transducer_pars.n_elements
            fprintf('Computing weights for element %i...', ind);
            element_weights = karray.getElementGridWeights(kgrid, ind);  
            el_binary_mask = element_weights ~= 0;
            binary_mask = binary_mask | el_binary_mask;

            source_labels = source_labels + ind * el_binary_mask;

            fprintf(' done\n');
        end

        mask_ind = find(binary_mask);
        num_source_points = length(mask_ind);
        distributed_source_signal = zeros(num_source_points, Nt);
        if canUseGPU
            distributed_source_signal = gpuArray(distributed_source_signal);
        end

        % Create a mapping from 3D indices to linear indices in the mask
        [grid_i, grid_j, grid_k] = ind2sub([kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1)], mask_ind);
        mask_mapping = sparse(grid_i, grid_j + kgrid.Ny * (grid_k - 1), 1:num_source_points, ...
                             kgrid.Nx, kgrid.Ny * max(kgrid.Nz, 1));


        % Assign signals to elements based on weights
        for ind = 1:transducer_pars.n_elements
            fprintf('Assigning signals for element %i...', ind);

            element_weights = karray.getElementGridWeights(kgrid, ind); 
            el_binary_mask = element_weights ~= 0;

            if canUseGPU
                element_weights = gpuArray(element_weights);
            end

            % Find indices where this element contributes
            [el_i, el_j, el_k] = ind2sub([kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1)], find(el_binary_mask));
            
            % Map to positions in distributed_source_signal
            for p = 1:length(el_i)
                idx = full(mask_mapping(el_i(p), el_j(p) + kgrid.Ny * (el_k(p) - 1)));
                if idx > 0 % Should always be true if masks are constructed properly
                    distributed_source_signal(idx, :) = distributed_source_signal(idx, :) + ...
                        element_weights(el_i(p), el_j(p), el_k(p)) * cw_signal(ind, :);
                end
            end
            fprintf(' done\n');
        end

        % Assign binary mask and distributed signals to the source structure
        source.p_mask = binary_mask;
        source.p = distributed_source_signal;

        % Show the 3D binary mask of the transducer in the grid
        source_mask = karray.getArrayBinaryMask(kgrid);
        if parameters.n_sim_dims == 3
            h = figure;
            hold on;
            
            % Extract isosurface
            p = patch(isosurface(source_mask, 0.5));  % 0.5 is the threshold for binary mask
            isonormals(source_mask, p);  % Add surface normals for better visualization
            set(p, 'FaceColor', 'r', 'EdgeColor', 'none');  % Set color and edges
            camlight; lighting gouraud;  % Improve lighting
            
            axis equal;
            xlabel('x-grid'); ylabel('y-grid'); zlabel('z-grid');
            title('3D Visualization of Binary Mask');

            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_binary_transducer_mask%s.fig', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'fig')

            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_binary_transducer_mask%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);

        else
            h = figure('Name', 'Transducer in Simulation Grid');
            % 2D visualization
            imagesc(source_mask);
            axis equal tight;
            colormap(hot);
            title('Transducer Binary Mask in Simulation Grid');
            colorbar;
    
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_binary_transducer_mask%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
        end
    end
end
