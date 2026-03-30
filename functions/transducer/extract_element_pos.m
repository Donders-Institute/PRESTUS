function elem_pos_m = extract_element_pos(parameters, tp, trans_pos_m) 
% EXTRACT_ELEMENT_POS Extracts matrix transducer element positions from file.
%
% This function reads element position coordinates from an external file
% and converts them into expected element position coordinates. The 
% positions are defined in physical space (typically in millimetres) and
% are transformed into metres relative to the simulation grid.
%
% The resulting element coordinates are shifted to the desired transducer
% location within the k-Wave simulation grid.
%
% Input:
%   parameters  Global simulation parameters
%   tp        - Struct containing transducer parameters.
%   trans_pos_m  [3x1] transducer position in meters
%
% Output:
%   elem_pos_m - 3 x N matrix containing element positions in metres,
%                expressed in simulation coordinates.

    % Extract matrix transducer configuration from the transducer parameters
    matrix_tp = tp.array_shape.matrix;

    % Extract file-based matrix shape parameters for reading element positions
    file_ext = matrix_tp.matrix_shape.extract_from_file;

    % ---------------------------------------------------------------------
    % Read element positions from file
    % ---------------------------------------------------------------------

    tran_info = readtable(file_ext.file_path);

    row_start = file_ext.start_row;
    row_end   = row_start + file_ext.n_elements - 1;

    col_start = file_ext.start_col;
    
    phys_positions_mm = table2array( ...
    tran_info(row_start:row_end, col_start:col_start+2) ...
    );
    
    % ---------------------------------------------------------------------
    % Optional: randomly select subset of elements
    % ---------------------------------------------------------------------
    if file_ext.select_random_subset

        if file_ext.subset.random_seed
            rng('shuffle'); % seed RNG based on current time
        end

        phys_positions_mm = datasample(phys_positions_mm, ...
            file_ext.subset.subset_n_elements, ...
            1, ...
            'Replace', ...
            false);
    end
    
    % ---------------------------------------------------------------------
    % Convert positions from mm to metres and center coordinates
    % ---------------------------------------------------------------------
    ROC = matrix_tp.curved.curv_radius_mm;

    % Translate origin from [0,0,ROC] to [0,0,0], convert mm→m, and flip Z-axis
    phys_positions_m = (phys_positions_mm - [0, 0, ROC]) .* [1, 1, -1] / 1000;
    
    % ---------------------------------------------------------------------
    % Optional: project elements onto new radius of curvature
    % ---------------------------------------------------------------------

    if file_ext.project_on_new_ROC

        new_ROC_mm = file_ext.ROC_projection.new_ROC_mm;

        r = hypot(phys_positions_m(:,1), phys_positions_m(:,2));
        r_max = max(r);

        sag_old = ROC - sqrt(ROC^2 - r_max^2);
        sag_new = new_ROC_mm - sqrt(new_ROC_mm^2 - r_max^2);

        scale_factor = sag_new / sag_old;

        z_ref = max(phys_positions_m(:,3));
        z_offset = phys_positions_m(:,3) - z_ref;

        new_z = z_ref + z_offset * scale_factor;

        phys_positions_m(:,3) = new_z;

        % ---------------------------------------------------------------------
        % Optional: projection sanity check
        % ---------------------------------------------------------------------

        % Extract X, Y, Z for convenience
        x = phys_positions_m(:,1);
        y = phys_positions_m(:,2);
        z = phys_positions_m(:,3);
    
        % Solve least-squares sphere fit (optional, purely for information)
        A = [2*x, 2*y, 2*z, ones(size(x))];
        b = x.^2 + y.^2 + z.^2;
        params = A \ b;
    
        xc = params(1); yc = params(2); zc = params(3); c = params(4);
        R_fit = sqrt(xc^2 + yc^2 + zc^2 + c);
        fprintf('Best-fit radius = %.6f m\n', R_fit);
    
        % [DEBUG] visualize grid
        if parameters.debug == 1
            h = figure('Name', 'Transducer Element Distribution');
            hold on; grid on;

            % Plot original element positions
            scatter3(x, y, z, 50, 'b', 'filled', 'DisplayName', 'Original positions');

            % Plot projected element positions
            scatter3(x, y, new_z, 50, 'r', 'filled', 'DisplayName', 'Projected positions');

            axis equal;
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            title('Transducer Element Distribution (3D)');
            legend('Location', 'best');
            view([90 0]); % Top-down view

            % Save figure
            fig_name = sprintf('sub-%03d_%s_transducer_redone_element_distribution%s', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix);

            saveas(h, fullfile(parameters.debug_dir, [fig_name '.fig']));
            saveas(h, fullfile(parameters.debug_dir, [fig_name '.png']));

            close(h);
        end

    end

    % ---------------------------------------------------------------------
    % Translate positions into simulation grid coordinates
    % ---------------------------------------------------------------------
    elem_pos_m = phys_positions_m + trans_pos_m';

    elem_pos_m = elem_pos_m';

end