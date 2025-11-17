function [transducer_pars, karray] = setup_flat_circular_matrix_transducer(transducer_pars, kgrid, trans_pos, parameters, focus_pos, karray)
            tran_width = (transducer_pars.n_elem_row * transducer_pars.elem_width) + ((transducer_pars.n_elem_row - 1) * transducer_pars.elem_spacing_width);   % [m]
            tran_height = (transducer_pars.n_elem_col * transducer_pars.elem_height) + ((transducer_pars.n_elem_col - 1) * transducer_pars.elem_spacing_height); % [m]
            
            x_vec = linspace(-tran_width/2 + transducer_pars.elem_width/2, tran_width/2 - transducer_pars.elem_width/2, transducer_pars.n_elem_row);
            y_vec = linspace(-tran_height/2 + transducer_pars.elem_height/2, tran_height/2 - transducer_pars.elem_height/2, transducer_pars.n_elem_col);
            
            % Create meshgrid of center points
            [X, Y] = meshgrid(x_vec, y_vec);

            % Add offset and flatten
            X = X + kgrid.x_vec(trans_pos(1));
            Y = Y + kgrid.y_vec(trans_pos(2));
            

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

            dx = X(:) - kgrid.x_vec(trans_pos(1));
            dy = Y(:) - kgrid.y_vec(trans_pos(2));
            distances = sqrt(dx.^2 + dy.^2);
            
            % Logical mask for elements within aperture
            radius = transducer_pars.Elements_OD_mm(end) * 1e-3 / 2 ;
            mask = distances <= radius;

            % Select only elements within the circular aperture
            elem_pos = [X(:), Y(:), Z(:)];
            elem_pos = elem_pos(mask, :);  % Efficient filtering

            transducer_pars.n_elements = size(elem_pos, 1);
            transducer_pars.source_amp = transducer_pars.source_amp(1) * ones(1, transducer_pars.n_elements); 

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
            
            h = figure('Position', [100, 100, 800, 600]);
            hold on;
            
            source_phase_rad = zeros(1, transducer_pars.n_elements);
            elem_pos = elem_pos';
            for ind = 1:size(elem_pos, 2)
                distance = sqrt((elem_pos(1, ind) - focus_pos_m(1))^2 + (elem_pos(2, ind) - focus_pos_m(2))^2 + (elem_pos(3, ind) - focus_pos_m(3))^2);
                            
                source_phase_rad(ind) = mod(k * distance, 2 * pi);
                
                karray.addRectElement(elem_pos(:, ind), transducer_pars.elem_height, transducer_pars.elem_width, [0, 0, 0])
                
                 % Get the center position of this element
                center_x = elem_pos(1, ind);
                center_y = elem_pos(2, ind);
                
                % Calculate the corner positions (x,y) for the rectangle
                x_min = center_x - transducer_pars.elem_width/2;
                y_min = center_y - transducer_pars.elem_height/2;

                % Draw the rectangle
                rectangle('Position', [x_min, y_min, transducer_pars.elem_width, transducer_pars.elem_height], ...
                         'EdgeColor', 'blue', 'FaceColor', [0.8, 0.8, 1], 'LineWidth', 1);
            end

            % Then, plot the center points
            scatter(elem_pos(1,:), elem_pos(2,:), 20, 'red', 'filled');

            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_transducer_elements_and_spacing%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
           
            transducer_pars.source_phase_rad = source_phase_rad;
            transducer_pars.source_phase_deg = rad2deg(source_phase_rad);
end