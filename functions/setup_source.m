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
    if parameters.unique_tran_design == 1
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

        if parameters.unique_tran_design == 1  
            % Add annular array elements to the kWaveArray object
            karray.addAnnularArray([kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))], ...
                                   transducer_pars.curv_radius_mm * 1e-3, ...
                                   [transducer_pars.Elements_ID_mm; transducer_pars.Elements_OD_mm] * 1e-3, ...
                                   [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]);

        %% new design
        elseif parameters.unique_tran_design == 2
            [transducer_pars, karray] = setup_dolphin_transducer_design(transducer_pars, kgrid, trans_pos, focus_pos, karray, parameters.medium.water.sound_speed, parameters.n_sim_dims, parameters.display_tran);
        elseif parameters.unique_tran_design == 3    
            tran_width = (transducer_pars.n_elem_row * transducer_pars.elem_width) + ((transducer_pars.n_elem_row - 1) * transducer_pars.elem_spacing_width);   % [m]
            tran_height = (transducer_pars.n_elem_col * transducer_pars.elem_height) + ((transducer_pars.n_elem_col - 1) * transducer_pars.elem_spacing_height); % [m]
            
            x_vec = linspace(-tran_width/2 + transducer_pars.elem_width/2, tran_width/2 - transducer_pars.elem_width/2, transducer_pars.n_elem_row);
            y_vec = linspace(-tran_height/2 + transducer_pars.elem_height/2, tran_height/2 - transducer_pars.elem_height/2, transducer_pars.n_elem_col);
            
            % Create meshgrid of center points
            [X, Y] = meshgrid(x_vec, y_vec);

            % Add offset and flatten
            X = X + kgrid.x_vec(trans_pos(1));
            Y = Y + kgrid.y_vec(trans_pos(2));
            Z = kgrid.z_vec(trans_pos(3)) * ones(size(X)); % For 2D array in XY plane

            % Combine into Nx3 array of center points
            rect_pos = [X(:), Y(:), Z(:)]';

            h = figure;
            scatter(X(:), Y(:), 60, 'filled'); % 60 is marker size
            axis equal;
            xlabel('X [m]');
            ylabel('Y [m]');
            title('Grid - Element Center Positions');
            grid on;

            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_transducer_grid%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);

            % Initialize an empty array to store element positions
            elem_pos = [];
            
            radius = transducer_pars.Elements_OD_mm(end) * 1e-3 / 2;
            % Loop over the grid points and check if they fall within the circular region
            for i = 1:size(rect_pos, 2)
                % Shifted distance from new center
                dx = rect_pos(1,i) - kgrid.x_vec(trans_pos(1));
                dy = rect_pos(2,i) - kgrid.y_vec(trans_pos(2));
                distance = sqrt(dx^2 + dy^2);
                
                % If the point is within the circle's radius, keep it
                if distance <= radius
                    elem_pos = [elem_pos; rect_pos(1,i), rect_pos(2,i), rect_pos(3,i)];
                end

            end

            h = figure;
            scatter(elem_pos(:, 1)*1e3, elem_pos(:, 2)*1e3, 60, 'filled'); % 60 is marker size
            axis equal;
            xlabel('X [mm]');
            ylabel('Y [mm]');
            title('Circular Aperture - Element Center Positions');
            grid on;

            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_circular_transducer_grid%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);

             % Define focus in meters
            focus_pos_m = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]';
            
            % Calculate wavelength and wavenumber for phase calculation in
            % water medium
            lambda = parameters.medium.water.sound_speed / transducer_pars.source_freq_hz; % [m]
            k = 2 * pi / lambda; % [rad/m]

            source_phase_deg = zeros(1, transducer_pars.n_elements);
            elem_pos = elem_pos';
            for ind = 1:size(elem_pos, 2)
                distance = sqrt((elem_pos(1, ind) - focus_pos_m(1))^2 + (elem_pos(2, ind) - focus_pos_m(2))^2 + (elem_pos(3, ind) - focus_pos_m(3))^2);
                            
                source_phase_deg(ind) = mod(k * distance, 2 * pi);
                
                karray.addRectElement(elem_pos(:, ind), transducer_pars.elem_height, transducer_pars.elem_width, [0, 0, 0])
            end
           
            transducer_pars.source_phase_deg = source_phase_deg;
            transducer_pars.source_phase_rad = deg2rad(source_phase_deg);

        elseif parameters.unique_tran_design == 4
            elem_pos = makeCartBowl( ...
                [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))], ...
                transducer_pars.curv_radius_mm * 1e-3, ...
                transducer_pars.Elements_OD_mm(end) * 1e-3, ...
                [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))], ...
                transducer_pars.n_elements, ...
                true);
            
            % Define focus in meters
            focus_pos_m = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]';

            % Calculate wavelength and wavenumber for phase calculation in
            % water medium
            lambda = parameters.medium.water.sound_speed / transducer_pars.source_freq_hz; % [m]
            k = 2 * pi / lambda; % [rad/m]

            source_phase_deg = zeros(1, transducer_pars.n_elements);

            % Initialize tx, ty, tz
            scaled_vectors = zeros(3, size(elem_pos, 2));
            tx = zeros(1, size(elem_pos, 2));
            ty = zeros(1, size(elem_pos, 2));
            tz = zeros(1, size(elem_pos, 2));

            % Add elements to the array
            for ind = 1:transducer_pars.n_elements
                el_pos = elem_pos(:, ind);
                scaled_vec = focus_pos_m - el_pos;  % vector from element to focus
                norm_vec = scaled_vec / norm(scaled_vec);            % normalize to unit vector

                scaled_vectors(:, ind) = scaled_vec;

                % Convert direction vector to Euler angles
                % Assumes initial normal of element is along +z
                % Compute rotation that aligns [0; 0; 1] to vec
                z0 = [0; 0; 1];
                v = cross(z0, norm_vec);
                s = norm(v);
                c = dot(z0, norm_vec);
                vx = [  0    -v(3)  v(2);
                       v(3)   0    -v(1);
                      -v(2)  v(1)   0  ];
                R = eye(3) + vx + vx^2 * ((1 - c) / (s^2 + eps));  % rotation matrix
                
                % Convert rotation matrix to ZYX Euler angles (extrinsic)
                yaw   = atan2d(R(2,1), R(1,1));  % around z
                pitch = atan2d(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));  % around y
                roll  = atan2d(R(3,2), R(3,3));  % around x
              
                tx(ind) = roll;
                ty(ind) = pitch;
                tz(ind) = yaw;

                karray.addRectElement(el_pos, transducer_pars.elem_height, transducer_pars.elem_width, [roll, pitch, yaw])

                distance = sqrt((elem_pos(1, ind) - focus_pos_m(1))^2 + (elem_pos(2, ind) - focus_pos_m(2))^2 + (elem_pos(3, ind) - focus_pos_m(3))^2);
                            
                source_phase_deg(ind) = mod(k * distance, 2 * pi);

            end

            transducer_pars.source_phase_deg = source_phase_deg;
            transducer_pars.source_phase_rad = deg2rad(source_phase_deg);

            h = figure;
            hold on;
            axis equal;
            
            scatter3(elem_pos(1,:), elem_pos(2,:), elem_pos(3,:), 'b.');
            plot3(focus_pos_m(1), focus_pos_m(2), focus_pos_m(3), 'go', 'MarkerFaceColor', 'g');
            
            for ind = 1:transducer_pars.n_elements
                el_pos = elem_pos(:, ind);
            
                % Ground truth direction (focus - element)
                vec_truth = focus_pos_m - el_pos;
                vec_truth = vec_truth / norm(vec_truth);
            
                % Get Euler angles (in degrees)
                roll  = tx(ind);
                pitch = ty(ind);
                yaw   = tz(ind);
            
                % Convert to radians
                r = deg2rad(roll);
                p = deg2rad(pitch);
                y = deg2rad(yaw);
            
                % Compute rotation matrix (extrinsic ZYX)
                Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];
                Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
                Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
            
                R = Rz * Ry * Rx;
            
                % Rotate default normal to get predicted direction
                vec_rot = R * [0; 0; 1];

                scale_factor = focus_pos_m(3) - elem_pos(3, ind); 

                % Scale the vectors for visualization
                vec_rot_scaled = vec_rot * scale_factor;
                vec_truth_scaled = vec_truth * scale_factor;
            
                % Plot predicted orientation (from Euler angles)
                quiver3(el_pos(1), el_pos(2), el_pos(3), ...
                        vec_rot_scaled(1), vec_rot_scaled(2), vec_rot_scaled(3), 0, 'r', 'LineWidth', 1.5);
            
                % Plot ground truth orientation (focus - element)
                quiver3(el_pos(1), el_pos(2), el_pos(3), ...
                        vec_truth_scaled(1), vec_truth_scaled(2), vec_truth_scaled(3), 0, 'g--', 'LineWidth', 1.2);
            end
            
            view(0, 0);
            xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
            title('Validation of Element Orientations: Red = Computed, Green = Ground Truth');
            legend('Element Position', 'Focus', 'From Euler Angles', 'Ground Truth');

            output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_transducer_element_orientation%s.fig', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'fig')
            
            output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_transducer_element_orientation%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')

            close(h);

        else
            fprintf("Transducer design number %i hasn't been implemented. \n", parameters.unique_tran_design)
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
            xlabel('x'); ylabel('y'); zlabel('z');
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
