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
end