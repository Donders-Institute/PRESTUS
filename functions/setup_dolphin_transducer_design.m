function [transducer_pars, karray] = setup_dolphin_transducer_design(transducer_pars, kgrid, trans_pos, focus_pos, karray, water_c0, n_sim_dims)
            [totalMask] = setup_dolphin_mask(transducer_pars);
            totalMaskCols = size(totalMask, 2);  % Get the number of columns in the new mask
                        
            % Count active elements and their indices
            active_elements = find(totalMask);
            num_elems = length(active_elements);

            % Calculate wavelength and wavenumber for phase calculation in
            % water medium
            lambda = water_c0 / transducer_pars.source_freq_hz; % [m]
            k = 2 * pi / lambda; % [rad/m]

            % Convert element positions to grid indices
            tran_center_m = [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))];

            focus_m = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))];

            % Calculate the updated center X position for the totalMask
            totalMaskCenterX = (totalMaskCols + 1) / 2;
            
            % Element counter for numbering
            element_counter = 1;
            
            phases_rad = zeros(1, num_elems);
            % Process each element in the grid using the total mask dimensions
            for row = 1:transducer_pars.n_elem_row
                for col = 1:totalMaskCols
                    % Check if this element is active according to the mask
                    if totalMask(row, col)       
                        % Calculate physical position (from center of array)
                        phys_x = tran_center_m(1) + (col - totalMaskCenterX) * (transducer_pars.elem_width + transducer_pars.elem_spacing_width);  % [m]
                        phys_y = tran_center_m(2) + (row - centerY) * (transducer_pars.elem_height + transducer_pars.elem_spacing_height); % [m]
                        phys_z = tran_center_m(3);      
            
                        distance = sqrt((phys_x - focus_m(1))^2 + (phys_y - focus_m(2))^2 + (phys_z - focus_m(3))^2);
                            
                        phaseDifference = mod(k * distance, 2 * pi);
                        
                        phases_rad(element_counter) = phaseDifference;

                        % Add rectangular element to array with correct position
                        % Parameters: position, width, height, depth, normal direction
                        karray.addRectElement([phys_x, phys_y, phys_z], transducer_pars.elem_width, transducer_pars.elem_height, [0, 0, 0]);
                    end
                end
            end
 
            % plot geometry
            % 1. Visualize the mask pattern
            figure('Name', 'Transducer Mask Pattern');
            imagesc(totalMask);
            title('Element Layout Pattern');
            colormap(gray);
            axis equal tight;
            xlabel('Column Index');
            ylabel('Row Index');
            colorbar;
            
            % % 3. Visualize individual element positions in 3D space
            % figure('Name', 'Element Positions');
            % 
            % % Get total number of elements
            % n_elements = transducer_pars.n_elements;
            % 
            % % Extract element positions (you would need to store these when creating elements)
            % positions = zeros(n_elements, 3);
            % for i = 1:n_elements
            %     % This assumes you've stored positions when creating elements
            %     % You might need to adapt this based on how you're tracking element positions
            %     positions(i,:) = karray.getElementPosition(i);
            % end
            % 
            % % Plot element positions
            % scatter3(positions(:,1), positions(:,2), positions(:,3), 50, 1:n_elements, 'filled');
            % axis equal;
            % xlabel('X [m]');
            % ylabel('Y [m]');
            % zlabel('Z [m]');
            % title('Element Positions');
            % colorbar('Title', 'Element Index');
            % grid on;
            % 
            % % 4. Visualize phase distribution (if available)
            % if isfield(transducer_pars, 'source_phase_rad')
            %     figure('Name', 'Element Phases');
            %     phases = transducer_pars.source_phase_rad;
            % 
            %     % Plot phases as a colormap on elements
            %     subplot(2,1,1);
            %     scatter3(positions(:,1), positions(:,2), positions(:,3), 50, phases, 'filled');
            %     axis equal;
            %     xlabel('X [m]');
            %     ylabel('Y [m]');
            %     zlabel('Z [m]');
            %     title('Element Phases');
            %     colorbar('Title', 'Phase [rad]');
            %     colormap(hsv);
            %     grid on;
            % 
            %     % Also show as a histogram
            %     subplot(2,1,2);
            %     histogram(phases, 20);
            %     xlabel('Phase [rad]');
            %     ylabel('Count');
            %     title('Phase Distribution');
            %end
            
            % 5. Test beam forming by visualizing the pressure field at a reference plane
            % This would typically be done separately in a simulation
            disp('To verify focusing, run a test simulation with these parameters and visualize the pressure field.');
        
            % calculate phases
            transducer_pars.source_phase_rad = phases_rad;
            transducer_pars.source_phase_deg = rad2deg(phases_rad);
end