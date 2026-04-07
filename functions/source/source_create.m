function [source, source_labels, transducer_pars] = source_create(parameters, kgrid, trans_pos, focus_pos)

% SOURCE_CREATE Creates a source for k-Wave simulations based on transducer parameters.
%
% Supports:
%   - multiple transducers, each with multiple elements (use_kWaveArray == 0)
%   - single transducer with kWaveArray-based setup (use_kWaveArray ~= 0)
%
% Input:
%   parameters      - struct containing simulation parameters (must contain parameters.transducer)
%   kgrid           - k-Wave grid struct (e.g., kWaveGrid)
%   trans_pos       - [transducer_N x n_sim_dims] or [1 x n_sim_dims] / [n_sim_dims x 1]
%   focus_pos       - same shape rules as trans_pos
%
% Output:
%   source          - struct with source.p_mask and source.p
%   source_labels   - grid of integer labels identifying active source regions
%   transducer_pars - struct array of transducer parameters with grid-based fields

    transducer_N = numel(parameters.transducer);

    if parameters.grid.use_kWaveArray ~= 0 && transducer_N > 1
        error(['Multiple transducers with kWaveArray (use_kWaveArray ~= 0) ' ...
               'are not implemented yet. Use use_kWaveArray == 0 or a single transducer.']);
    end

    %% Validate / broadcast positions

    if isequal(size(trans_pos), [numel(parameters.grid.dims) 1]), trans_pos = trans_pos'; end
    if isequal(size(focus_pos), [numel(parameters.grid.dims) 1]), focus_pos = focus_pos'; end

    if size(trans_pos,1) == 1,  trans_pos  = repmat(trans_pos,  transducer_N, 1); end
    if size(focus_pos,1) == 1,  focus_pos  = repmat(focus_pos, transducer_N, 1); end

    if size(trans_pos,1) ~= transducer_N || size(focus_pos,1) ~= transducer_N
        error('trans_pos and focus_pos must have one row per transducer.');
    end
    if size(trans_pos,2) ~= numel(parameters.grid.dims) || size(focus_pos,2) ~= numel(parameters.grid.dims)
        error('trans_pos and focus_pos columns must match numel(parameters.grid.dims).');
    end

    %% Convert element diameters from mm to grid points (for all transducers)

    transducer_pars = parameters.transducer;

    for it = 1:transducer_N
        tp = transducer_pars(it);

        if strcmp(tp.type, 'annular')
            tp.annular.Elements_OD = 2 * floor(tp.annular.Elements_OD_mm / parameters.grid.resolution_mm / 2) + 1;
            tp.annular.Elements_ID = 2 * floor(tp.annular.Elements_ID_mm / parameters.grid.resolution_mm / 2) + 1;
            tp.annular.Elements_ID(tp.annular.Elements_ID_mm == 0) = 0;
			
        end
        
        tp.(tp.type).radius_grid = round(tp.(tp.type).curv_radius_mm / parameters.grid.resolution_mm);

        if it == 1
            % initialise struct array with full field set of tp
            transducer_pars = repmat(tp, 1, transducer_N);
        end

        transducer_pars(it) = tp;
    end

    source = struct();

    %% Branch 1: Custom element geometry (non-kWaveArray setup, multi-transducer)

    if parameters.grid.use_kWaveArray == 0

        % --- per-transducer CW signals (per-element) ---

        freqs = [transducer_pars.(transducer_pars(1).type).source_freq_hz];
        if numel(unique(freqs)) > 1
            warning('Multiple source frequencies not yet supported. Using %i Hz.', freqs(1));
        end

        cw_signals_all    = cell(1, transducer_N);
        n_elements_per_T  = zeros(1, transducer_N);

        for it = 1:transducer_N
            tp = transducer_pars(it);

            % source_amp and source_phase_rad are per-element vectors
            cw_signals_all{it} = createCWSignals( ...
                kgrid.t_array, ...
                freqs(1), ...               % use canonical frequency
                tp.(tp.type).source_amp, ...
                tp.(tp.type).source_phase_rad);       % [n_elements x Nt]

            n_elements_per_T(it) = tp.n_elements;
        end

        % global element indexing: each (transducer, element) → unique index
        offsets          = [0, cumsum(n_elements_per_T(1:end-1))];  % length transducer_N
        n_elements_total = sum(n_elements_per_T);

        % --- geometry: per-element bowls, numeric source_labels = global element index ---

        grid_dims       = parameters.grid.dims;
        transducer_mask = false(grid_dims);
        source_labels   = zeros(grid_dims);
        
        for it = 1:transducer_N
            tp          = transducer_pars(it);
            trans_pos_i = transducer_pars(it).position.trans_pos;
            focus_pos_i = transducer_pars(it).position.focus_pos;
        
            for el_i = 1:tp.(tp.type).n_elements
                global_el = offsets(it) + el_i;  % 1..n_elements_total
        
                % outer element aperture
                if numel(parameters.grid.dims) == 3
                    bowl = makeBowl(grid_dims, trans_pos_i, tp.(tp.type).radius_grid, ...
                                    tp.(tp.type).Elements_OD(el_i), focus_pos_i);
                else
                    bowl = makeArc(grid_dims, trans_pos_i, tp.(tp.type).radius_grid, ...
                                   tp.(tp.type).Elements_OD(el_i), focus_pos_i);
                end

                % subtract inner aperture if applicable
                if tp.(tp.type).Elements_ID(el_i) > 0
                    if numel(parameters.grid.dims) == 3
                        bowl = bowl - makeBowl(grid_dims, trans_pos_i, ...
                                               tp.(tp.type).radius_grid, tp.(tp.type).Elements_ID(el_i), focus_pos_i);
                    else
                        bowl = bowl - makeArc(grid_dims, trans_pos_i, ...
                                              tp.(tp.type).radius_grid, tp.(tp.type).Elements_ID(el_i), focus_pos_i);
                    end
                end
        
                bmask = bowl > 0;
        
                % accumulate mask, but ASSIGN label instead of summing
                transducer_mask(bmask) = true;
                source_labels(bmask)   = global_el;
            end
        end

        % --- assemble global CW signal matrix [n_elements_total x Nt] ---

        Nt        = size(cw_signals_all{1}, 2);
        cw_global = zeros(n_elements_total, Nt);

        for it = 1:transducer_N
            idx0 = offsets(it) + 1;
            idx1 = offsets(it) + n_elements_per_T(it);
            cw_global(idx0:idx1, :) = cw_signals_all{it};  % rows = elements of transducer it
        end

        % --- assign pressure signals: one row per voxel, mapped by element index ---

        p_mask_source_p = source_labels(source_labels > 0);
        p_mask_source_p = reshape(p_mask_source_p, [], 1);

        source.p = zeros(length(p_mask_source_p), Nt);
        for ii = 1:length(p_mask_source_p)
            el_idx       = p_mask_source_p(ii);  % global element index
            source.p(ii, :) = cw_global(el_idx, :);
        end

        source.p_mask    = transducer_mask;

    %% Branch 2: kWaveArray-based setup (single transducer)

    else
        disp('Setting up kWaveArray (might take a bit of time)');

        % Note: Both `parameters` and `tp` are passed to this function to
        % allow support for multiple transducers in future implementations.
        tp = transducer_pars(1);  % use first (and only) transducer here

        % 3D/2D positions for kWaveArray are taken from the first transducer
        trans_pos_1  = tp.position.trans_pos;
        focus_pos_1  = tp.position.focus_pos;

        % Determine if axisymmetric mode should be enabled
        if numel(parameters.grid.dims) == 2 && ...
                isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
            axisymmetric = true;
            disp("Using axisymmetric setup for 2D input...");
        else
            axisymmetric = false;
        end

        % Create kWaveArray object (always without axisymmetry)
        karray = kWaveArray('Axisymmetric', false, ...
                           'BLITolerance', 0.1, ...
                           'UpsamplingRate', 10, ...
                           'BLIType', 'sinc');

        switch tp.type
            case 'annular'
                % Set focus position and transducer position vectors in physical coordinates
                if numel(parameters.grid.dims) == 3
                    % 3D annular array
                    pos_vec   = [kgrid.x_vec(trans_pos_1(1)), kgrid.y_vec(trans_pos_1(2)), kgrid.z_vec(trans_pos_1(3))];
                    focus_vec = [kgrid.x_vec(focus_pos_1(1)), kgrid.y_vec(focus_pos_1(2)), kgrid.z_vec(focus_pos_1(3))];

                    karray.addAnnularArray(pos_vec, ...
                                           tp.annular.curv_radius_mm * 1e-3, ...
                                           [tp.annular.Elements_ID_mm; tp.annular.Elements_OD_mm] * 1e-3, ...
                                           focus_vec);

                elseif numel(parameters.grid.dims) == 2 && axisymmetric == false

                    % 2D arc-shaped element
                    pos_vec   = [kgrid.x_vec(trans_pos_1(1)), kgrid.y_vec(trans_pos_1(2))];
                    focus_vec = [kgrid.x_vec(focus_pos_1(1)), kgrid.y_vec(focus_pos_1(2))];

                    karray.addArcElement(pos_vec, ...
                                         tp.annular.curv_radius_mm * 1e-3, ...
                                         tp.annular.Elements_OD_mm * 1e-3, ...
                                         focus_vec);
                end
            case 'matrix'
                matrix_tp = tp.matrix;
                
                trans_pos_m = [kgrid.x_vec(tp.position.trans_pos(1));
                    kgrid.y_vec(tp.position.trans_pos(2));
                    kgrid.z_vec(tp.position.trans_pos(3))];

                focus_pos_m = [kgrid.x_vec(tp.position.focus_pos(1));
                    kgrid.y_vec(tp.position.focus_pos(2));
                    kgrid.z_vec(tp.position.focus_pos(3))];

                switch matrix_tp.matrix_shape.type
                    case 'define_here'
                        [elem_pos_m, tp] = convert_to_element_pos(parameters, tp, trans_pos_m, focus_pos_m);
                    case 'extract_from_file'               
                        elem_pos_m = extract_element_pos(parameters, tp,  trans_pos_m);
                    otherwise 
                        error('Matrix shape %s is unknown or not implemented.', matrix_tp.matrix_shape.type)
                end

                % [DEBUG] visualize element distribution
                if parameters.simulation.debug == 1
                    % Convert positions to mm for plotting
                    elem_pos_mm = elem_pos_m' * 1e3;

                    h = figure;
                    scatter3(elem_pos_mm(:,1), elem_pos_mm(:,2), elem_pos_mm(:,3), 60, 'filled');
                    axis equal
                    xlabel('X [mm]')
                    ylabel('Y [mm]')
                    zlabel('Z [mm]')
                    view([0 90])
                    title('Transducer Element Distribution')
                    grid on;

                    % Build filenames
                    fig_filename = fullfile(parameters.io.debug_dir, ...
                        sprintf('sub-%03d_%s_transducer_element_distribution%s.fig', ...
                        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

                    png_filename = fullfile(parameters.io.debug_dir, ...
                        sprintf('sub-%03d_%s_transducer_element_distribution%s.png', ...
                        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

                    % Save outputs
                    saveas(h, fig_filename, 'fig')
                    saveas(h, png_filename, 'png')
                    close(h)
                end
                
                [karray, tp] = create_matrix_karray(kgrid, karray, parameters, tp, elem_pos_m, trans_pos, focus_pos);

            otherwise
                error('Array shape %s is unknown or not implemented.', tp.type)
        end

        % CW per-element for this transducer
        cw_signal = createCWSignals( ...
            kgrid.t_array, ...
            tp.(tp.type).source_freq_hz, ...
            tp.(tp.type).source_amp, ...
            tp.annular.source_phase_rad);   % [n_elements x Nt]

        if axisymmetric == true && strcmp(tp.type, 'annular')

            kgrid_mirrored = kWaveGrid(kgrid.Nx, kgrid.dx, 2*kgrid.Ny - 1, kgrid.dy);
            karray_full = kWaveArray('Axisymmetric', false, 'BLITolerance', 0.01, 'UpsamplingRate', 100);

            x_offset     = (trans_pos_1(1)-1) * (1/parameters.grid.resolution_mm);
            position_base = [kgrid.x_vec(1) + x_offset*kgrid.dx, 0+eps];
            focus_pos_full = [0, 0+eps];
            
            for el_i = 1:tp.annular.n_elements
                el_OD_m = tp.annular.Elements_OD_mm(el_i) * 1e-3;
                y_shift = (el_i - (tp.annular.n_elements+1)/2) * el_OD_m;
                element_pos = position_base + [0, y_shift];

                karray_full.addArcElement(element_pos, tp.curv_radius_mm * 1e-3, el_OD_m, focus_pos_full);
            end

            source_mask_full = false(kgrid_mirrored.Nx, kgrid_mirrored.Ny);
            grid_weights_3d  = zeros(tp.n_elements, kgrid_mirrored.Nx, kgrid_mirrored.Ny);

            for el_i = 1:tp.n_elements
                grid_weights_3d(el_i, :, :) = karray_full.getElementGridWeights(kgrid_mirrored, el_i);
                % Update the mask
                source_mask_full = source_mask_full | (squeeze(grid_weights_3d(el_i, :, :)) > 0.25);
            end

            source.p_mask = source_mask_full(:, kgrid.Ny:end);
            grid_weights_cropped = grid_weights_3d(:, :, kgrid.Ny:end);

            source_idx = find(source.p_mask);
            num_points = numel(source_idx);
            nTime      = size(cw_signal, 2);

            source.p = zeros(num_points, nTime);

            for el_i = 1:tp.n_elements
                weights_el = squeeze(grid_weights_cropped(el_i, :, :));
                weights_vec = weights_el(source_idx);
                source.p = source.p + bsxfun(@times, weights_vec, cw_signal(el_i, :));
            end

            % label sources at 50% weights
            source_labels = sum(grid_weights_cropped > 0.5, 1) > 0;
            source_labels = squeeze(source_labels);

        else
            % Non-axisymmetric: compute weights and distribute signals manually
            grid_weights_4d = zeros(tp.annular.n_elements, kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));

            for ind = 1:tp.annular.n_elements
                fprintf('Computing weights for element %i...', ind);
                grid_weights_4d(ind,:,:,:) = karray.getElementGridWeights(kgrid, ind);
                fprintf(' done\n');
            end

            binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;

            Nt        = size(cw_signal, 2);
            mask_ind  = find(binary_mask);
            num_src   = sum(binary_mask(:));

            distributed_source_signal = zeros(num_src, Nt);
            source_labels             = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));

            if canUseGPU
                cw_signal                 = gpuArray(cw_signal);
                distributed_source_signal = gpuArray(distributed_source_signal);
            end

            for ind = 1:tp.annular.n_elements
                source_weights = squeeze(grid_weights_4d(ind,:,:,:));
                el_binary_mask = source_weights ~= 0;

                if canUseGPU
                    source_weights = gpuArray(source_weights);
                end

                element_mask_ind = find(el_binary_mask);
                local_ind        = ismember(mask_ind, element_mask_ind);

                distributed_source_signal(local_ind,:) = ...
                    distributed_source_signal(local_ind,:) + ...
                    source_weights(element_mask_ind) * cw_signal(ind,:);

                source_labels = source_labels + ind * el_binary_mask;
            end

            source.p_mask = binary_mask;
            source.p      = distributed_source_signal;

            show_binary_mask_transducer(karray, kgrid, parameters);
        end
        transducer_pars(1) = tp;
    end
end
