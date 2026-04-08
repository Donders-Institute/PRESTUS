function [elem_pos_m, tr] = convert_to_element_pos(parameters, tr, trans_pos_m, focus_pos_m)
%CONVERT_TO_ELEMENT_POS Generate element positions for matrix transducers
%
% This function generates the 3-D element coordinates of a matrix
% array transducer based on the configuration defined in the parameter
% structure. Supported grid layouts are:
%
%   • rectangular grid
%   • Fibonacci spiral grid
%   • Fermat spiral grid (via makeCartBowl)
%
% The resulting element coordinates are optionally projected onto a
% spherical cap when the transducer is defined as curved.
%
% INPUTS
%   parameters  Global simulation parameters
%   tr          Transducer parameter structure
%   trans_pos_m  [3x1] transducer position in meters
%   focus_pos_m  [3x1] focus position in meters
%
% OUTPUTS
%   elem_pos_m  Element center coordinates [m] (3 × N)
%   tr          Updated transducer parameter structure
%
% NOTES
%   • Element coordinates are returned in k-Wave format (3 × N)
%   • Internal calculations are performed in millimeters unless noted
%   • Curvature is applied only when tr.matrix.is_curved = true

    % ----------------------------------------------------------------------
    % Extract configuration
    % -----------------------------------------------------------------------

    % Extract defined matrix shape parameters for reading element positions
    defined = tr.matrix.matrix_shape.define_here;
    grid_shape = defined.grid_shape;

    % Ensure column vectors
    trans_pos_m = trans_pos_m(:);
    focus_pos_m = focus_pos_m(:);

    switch grid_shape.type
        case 'fibonacci'

            fib = grid_shape.fibonacci;

            N = fib.elem_n;              % number of elements
            D = tr.matrix.outer_diameter_mm; % aperture diameter [mm]
            R = tr.matrix.curv_radius_mm;    % radius of curvature [mm]

            a = D/2;                         % aperture radius [mm]
            ga = pi*(3 - sqrt(5));           % golden angle

            elem_pos_mm = zeros(N,3);

            for i = 1:N

                % Uniform point distribution inside aperture
                r     = a * sqrt((i - 0.5)/N);
                theta = ga * (i - 1);

                x = r*cos(theta);
                y = r*sin(theta);

                % Project onto spherical cap if curved
                if tr.matrix.is_curved
                    z = R - sqrt(R^2 - x^2 - y^2);
                else
                    z = 0;
                end

                elem_pos_mm(i,:) = [x y z];

            end

            % --- Estimate average element pitch -------------------------------
            avg_dist = 0;
            for i = 1:N
                dists = vecnorm(elem_pos_mm - elem_pos_mm(i,:),2,2); % Euclidean distances
                dists(i) = inf; % Ignore self
                avg_dist = avg_dist + min(dists); % Nearest neighbor distance
            end

            avg_pitch = avg_dist / N; % Average center-to-center spacing

            kerf = fib.kerf_mm;                       % assumed kerf width [mm]
            element_size = avg_pitch - kerf;

            fprintf('Estimated avg pitch: %.3f mm\n', avg_pitch);
            fprintf('Recommended element size (kerf %.3f mm): %.3f mm\n', ...
                kerf, element_size);

            % Convert to meters + translate
            elem_pos_m = (elem_pos_mm / 1e3)' + trans_pos_m;

        case 'rect'

            rect = grid_shape.rect;

            % --- Compute array dimensions --------------------------------------
            tran_width = (rect.elem_n_row * tr.elem_width_mm) + ...
                ((rect.elem_n_row - 1) * rect.elem_spacing_width_mm);

            tran_height = (rect.elem_n_col * tr.elem_height_mm) + ...
                ((rect.elem_n_col - 1) * rect.elem_spacing_height_mm);

            % --- Element center coordinates ------------------------------------
            x_vec = linspace(-tran_width/2 + tr.elem_width_mm/2,...
                tran_width/2 - tr.elem_width_mm/2,...
                rect.elem_n_row);

            y_vec = linspace(-tran_height/2 + tr.elem_height_mm/2,...
                tran_height/2 - tr.elem_height_mm/2,...
                rect.elem_n_col);

            [X,Y] = meshgrid(x_vec,y_vec);

            % Flat array coordinates
            X = X + trans_pos_m(1);
            Y = Y + trans_pos_m(2);
            Z = trans_pos_m(3) * ones(size(X));

            % [DEBUG] visualize grid
            if parameters.debug == 1
                h = figure;
                scatter(X(:), Y(:), 60, 'filled')
                axis equal
                xlabel('X [m]')
                ylabel('Y [m]')
                title('Grid - Element Center Positions')
                grid on

                output_plot_filename = fullfile(parameters.io.debug_dir,...
                    sprintf('sub-%03d_%s_transducer_grid%s.png',...
                    parameters.subject_id,...
                    parameters.simulation.medium,...
                    parameters.io.output_affix));

                saveas(h,output_plot_filename,'png')
                close(h)
            end

            % Restrict to circular aperture
            dx = X(:) - trans_pos_m(1);
            dy = Y(:) - trans_pos_m(2);
            distances = sqrt(dx.^2 + dy.^2);

            radius = tr.matrix.elem_od_mm(end) * 1e-3 / 2;
            mask   = distances <= radius;

            elem_pos_m = [X(:) Y(:) Z(:)];
            elem_pos_m = elem_pos_m(mask,:);

            % Update number of elements
            tr.matrix.elem_n = size(elem_pos_m,1);

            % Duplicate source amplitude per element
            tr.elem_amp = tr.elem_amp(1) * ...
                ones(1,tr.matrix.elem_n);

            % [DEBUG] visualize circular aperture
            if parameters.debug == 1
                h = figure;
                scatter(elem_pos_m(:,1)*1e3,elem_pos_m(:,2)*1e3,60,'filled')
                axis equal
                xlabel('X [mm]')
                ylabel('Y [mm]')
                title('Circular Aperture - Element Center Positions')
                grid on

                saveas(h,output_plot_filename,'png')
                close(h)
            end

            % Apply curvature (spherical cap)
            if tr.matrix.is_curved

                ROC = tr.matrix.curv_radius_mm / 1000;
                R2  = ROC^2;

                sagitta_term = R2 - elem_pos_m(:,1).^2 - elem_pos_m(:,2).^2;

                if any(sagitta_term < 0)
                    error('Some elements fall outside the spherical cap.');
                end

                elem_pos_m(:,3) = elem_pos_m(:,3) + ROC - sqrt(sagitta_term);

            end

            elem_pos_m = elem_pos_m';

        case 'fermat'
            elem_pos_m = makeCartBowl( ...
                trans_pos_m', ...
                tr.matrix.curv_radius_mm * 1e-3, ...
                tr.matrix.outer_diameter_mm * 1e-3, ...
                focus_pos_m', ...
                tr.elem_n, ...
                true);

        otherwise
            error('Grid shape %s is unknown or not implemented.', grid_shape.type)
    end

end