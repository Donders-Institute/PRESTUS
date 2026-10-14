function [elem_pos_m, tr] = extract_element_pos(parameters, tr, trans_pos_m)
% EXTRACT_ELEMENT_POS  Read matrix transducer element positions from file
%
% Reads element position coordinates from an external file and converts
% them to metres relative to the simulation grid. Positions are
% optionally projected onto a new radius of curvature before being
% shifted to the transducer location in the k-Wave grid.
%
% Coordinate frame of the file: for a curved array (is_curved true), the
% origin is the centre of curvature and z points towards the array (the
% apex sits at z = curv_radius_mm). For a flat array (is_curved false),
% the origin is the array centre and z points along the beam; elements
% normally lie at z = 0.
%
% With extract_from_file.rot_col set, that column holds each element's
% rotation about its own normal [deg] (e.g. 90 for a rectangular element
% turned on its side); it is returned in tr.matrix.elem_rot_deg.
%
% Use as:
%   [elem_pos_m, tr] = extract_element_pos(parameters, tr, trans_pos_m)
%
% Input:
%   parameters  - (1,1) simulation parameters struct
%   tr          - (1,1) transducer parameter struct
%   trans_pos_m - [3x1] transducer position [m]
%
% Output:
%   elem_pos_m - [3xN] element positions in simulation coordinates [m]
%   tr         - transducer struct; tr.matrix.elem_rot_deg [1xN] set when
%                extract_from_file.rot_col is given
%
% See also: CONVERT_TO_ELEMENT_POS, TRANSDUCER_SETUP

    % Extract file-based matrix shape parameters for reading element positions
    file_ext = tr.matrix.matrix_shape.extract_from_file;

    % ---------------------------------------------------------------------
    % Read element positions from file
    % ---------------------------------------------------------------------

    tran_info = readtable(file_ext.file_path);

    row_start = file_ext.start_row;
    row_end   = row_start + file_ext.elem_n - 1;

    col_start = file_ext.start_col;
    
    phys_positions_mm = table2array( ...
    tran_info(row_start:row_end, col_start:col_start+2) ...
    );

    % Optional per-element rotation about the element normal [deg]
    elem_rot_deg = [];
    if isfield(file_ext, 'rot_col') && ~isempty(file_ext.rot_col)
        elem_rot_deg = table2array(tran_info(row_start:row_end, file_ext.rot_col));
    end
    
    % ---------------------------------------------------------------------
    % Optional: randomly select subset of elements
    % ---------------------------------------------------------------------
    if file_ext.select_random_subset

        if file_ext.subset.random_seed
            rng('shuffle'); % seed RNG based on current time
        end

        [phys_positions_mm, keep] = datasample(phys_positions_mm, ...
            file_ext.subset.subset_n_elements, ...
            1, ...
            'Replace', ...
            false);
        if ~isempty(elem_rot_deg); elem_rot_deg = elem_rot_deg(keep); end
    end
    if ~isempty(elem_rot_deg); tr.matrix.elem_rot_deg = elem_rot_deg(:)'; end
    
    % ---------------------------------------------------------------------
    % Convert positions from mm to metres and center coordinates
    % ---------------------------------------------------------------------
    ROC = tr.matrix.curv_radius_mm;

    if is_flat_matrix(tr.matrix)
        % Flat array: the file is already in the transducer frame (origin at
        % the array centre, beam along +z); only mm -> m
        phys_positions_m = phys_positions_mm / 1000;
        if isfield(file_ext, 'project_on_new_ROC') && isequal(file_ext.project_on_new_ROC, true)
            error('extract_element_pos:flatProjection', ...
                'project_on_new_ROC needs a curved array (is_curved true).');
        end
    else
        % Translate origin from [0,0,ROC] to [0,0,0], convert mm→m, and flip Z-axis
        phys_positions_m = (phys_positions_mm - [0, 0, ROC]) .* [1, 1, -1] / 1000;
    end
    
    % ---------------------------------------------------------------------
    % Optional: project elements onto new radius of curvature
    % ---------------------------------------------------------------------

    if ~is_flat_matrix(tr.matrix) && isequal(file_ext.project_on_new_ROC, true)

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
        if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'debug') && parameters.simulation.debug == 1
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
                parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix);

            saveas(h, fullfile(parameters.io.dir_debug_source, [fig_name '.fig']));
            saveas(h, fullfile(parameters.io.dir_debug_source, [fig_name '.png']));

            close(h);
        end

    end

    % ---------------------------------------------------------------------
    % Translate positions into simulation grid coordinates
    % ---------------------------------------------------------------------
    elem_pos_m = phys_positions_m + trans_pos_m';

    elem_pos_m = elem_pos_m';

end