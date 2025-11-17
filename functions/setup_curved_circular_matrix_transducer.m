function [transducer_pars, karray] = setup_curved_circular_matrix_transducer(kgrid, trans_pos, transducer_pars, focus_pos, parameters, karray, is_rect_grid, is_cut)
            if is_rect_grid
                % Step 1: Compute total width and height of the array
                tran_width = (transducer_pars.n_elem_row * transducer_pars.elem_width) + ...
                             ((transducer_pars.n_elem_row - 1) * transducer_pars.elem_spacing_width);
                
                tran_height = (transducer_pars.n_elem_col * transducer_pars.elem_height) + ...
                              ((transducer_pars.n_elem_col - 1) * transducer_pars.elem_spacing_height);
                
                % Step 2: Define x and y coordinates of element centers
                x_vec = linspace(-tran_width/2 + transducer_pars.elem_width/2, ...
                                  tran_width/2 - transducer_pars.elem_width/2, transducer_pars.n_elem_row);
                
                y_vec = linspace(-tran_height/2 + transducer_pars.elem_height/2, ...
                                  tran_height/2 - transducer_pars.elem_height/2, transducer_pars.n_elem_col);
                
                [X, Y] = meshgrid(x_vec, y_vec);
                
                % Step 4: Apply translation to position the transducer in the simulation grid
                X = X + kgrid.x_vec(trans_pos(1));
                Y = Y + kgrid.y_vec(trans_pos(2));
                Z = kgrid.z_vec(trans_pos(3)) * ones(size(X)); % For 2D array in XY plane 
    
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
                    sprintf('sub-%03d_%s_transducer_grid%s.png', ...
                    parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
                saveas(h, output_plot_filename, 'png')
                close(h);
    
                % Compute Z coordinates based on spherical surface (assuming a spherical cap)
                ROC = transducer_pars.curv_radius_mm / 1000;
                R2 = ROC^2;
                
                % Compute the sag (z-offset from flat plane)
                sagitta_term = R2 - elem_pos(:, 1).^2 - elem_pos(:, 2).^2;
                
                % Safety check to prevent imaginary numbers
                if any(sagitta_term < 0)
                    error('Some elements fall outside the spherical cap surface.');
                end
                
                % Apply curvature (spherical cap)
                elem_pos(:, 3) = elem_pos(:, 3) + ROC - sqrt(sagitta_term);
                elem_pos = elem_pos';
            else
                elem_pos = makeCartBowl( ...
                                [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))], ...
                                transducer_pars.curv_radius_mm * 1e-3, ...
                                transducer_pars.Elements_OD_mm(end) * 1e-3, ...
                                [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))], ...
                                transducer_pars.n_elements, ...
                                true);
            end



            h = figure;
            scatter3(elem_pos(1, :)*1e3, elem_pos(2, :)*1e3, elem_pos(3, :)*1e3, 60, 'filled');
            xlabel('X [mm]');
            ylabel('Y [mm]');
            zlabel('Z [mm]');
            title('3D Curved Circular Aperture');
            axis equal;
            grid on;

            output_plot_filename = fullfile(parameters.debug_dir, ...
                    sprintf('sub-%03d_%s_circular_transducer_grid%s.png', ...
                    parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
            
            if is_cut
            % Backup before cutting
                elem_pos_full = elem_pos;
    
                % Remove two parts of module to be able to combine modules
                x = elem_pos(1, :);
                y = elem_pos(2, :);
                
                % Remove left lobe
                cut1 = x < kgrid.x_vec(trans_pos(1)) - 23.695/1000;
                
                % Remove lower wedge-shaped region
                cx = kgrid.x_vec(trans_pos(1));
                cy = kgrid.y_vec(trans_pos(2));

                chord_len = 27.915e-3;
                center_dist = 23.695e-3;
                start_dist = 27.5e-3;
                chord_start = [cx, cy + start_dist];
                
                angle = asin(center_dist/start_dist);
                
                x_dist = sin(angle) * chord_len;
                y_dist = cos(angle) * chord_len;

                chord_end = chord_start + [-x_dist, y_dist];
                
                % Define the line vector
                chord_vec = chord_end - chord_start;
                
                % For each point, determine if it's on the "keep" side of the line
                % The cross product of (point - chord_start) and chord_vec will be positive
                % if the point is on one side and negative if on the other side
                cut2 = zeros(size(x));
                for i = 1:length(x)
                    point = [x(i), y(i)];
                    point_vec = point + chord_start;
                    % In 2D, the cross product is essentially: 
                    cross_prod = point_vec(1)*chord_vec(2) - point_vec(2)*chord_vec(1);
                    
                    % Negative cross product means the point is on the "remove" side (bottom left)
                    cut2(i) = (cross_prod < 0);
                end
                
                % Keep only the elements outside the clover cuts
                cut = cut1 | cut2;
                elem_pos_cut = elem_pos(:, ~cut);

                % === Plot Before and After ===
                h = figure;
                hold on;
                scatter3(elem_pos_full(1, :)*1e3, elem_pos_full(2, :)*1e3, elem_pos_full(3, :)*1e3, ...
                         60, [0.7 0.7 0.7], 'filled'); % full set in gray
                scatter3(elem_pos_cut(1, :)*1e3, elem_pos_cut(2, :)*1e3, elem_pos_cut(3, :)*1e3, ...
                         60, 'r', 'filled');           % cut set in red
                xlabel('X [mm]');
                ylabel('Y [mm]');
                zlabel('Z [mm]');
                legend({'Before Cut', 'After Cut'});
                title('Comparison: Full vs Cut 3D Curved Aperture');
                axis equal;
                grid on;
                hold off;
                
                output_plot_filename = fullfile(parameters.debug_dir, ...
                    sprintf('sub-%03d_%s_compare_transducer_grid%s.png', ...
                    parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
                saveas(h, output_plot_filename, 'png');
                close(h);
                
                % Update the main variable to the cut version
                elem_pos = elem_pos_cut;

                transducer_pars.n_elements = size(elem_pos, 2);
                transducer_pars.source_amp = transducer_pars.source_amp(1) * ones(1, transducer_pars.n_elements); 
            end

            % Define focus in meters
            focus_pos_m = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]';

            % Calculate wavelength and wavenumber for phase calculation in
            % water medium
            lambda = parameters.medium.water.sound_speed / transducer_pars.source_freq_hz; % [m]
            k = 2 * pi / lambda; % [rad/m]

            source_phase_rad = zeros(1, transducer_pars.n_elements);

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
                            
                source_phase_rad(ind) = mod(k * distance, 2 * pi);

            end

            transducer_pars.source_phase_rad = source_phase_rad;
            transducer_pars.source_phase_deg = rad2deg(source_phase_rad);

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

                scale_factor = norm(focus_pos_m - el_pos);

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
end