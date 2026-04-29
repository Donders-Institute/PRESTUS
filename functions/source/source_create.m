function [source, source_labels, tr_arr] = source_create(parameters, kgrid, trans_pos, focus_pos)
% SOURCE_CREATE  Build a k-Wave source struct from PRESTUS transducer parameters
%
% Dispatches to one of two branches:
%   Branch 1 (use_kWaveArray == 0): Supports multiple transducers. Each
%     annular element is placed with makeBowl / makeArc. CW signals are
%     built with createCWSignals and assigned row-by-row to source.p.
%   Branch 2 (use_kWaveArray ~= 0): Single transducer only. kWaveArray
%     handles element geometry; grid weights distribute signals. For 2-D
%     axisymmetric annular arrays a mirrored grid approach is used.
% Element diameters are converted from mm to odd-integer grid points before
% geometry creation.
%
% Use as:
%   [source, source_labels, tr_arr] = source_create(parameters, kgrid, trans_pos, focus_pos)
%
% Input:
%   parameters - PRESTUS config; must contain grid.dims, grid.resolution_mm [mm],
%                grid.use_kWaveArray, and transducer (struct array)
%   kgrid      - kWaveGrid object
%   trans_pos  - [transducer_N x n_dims] or [1 x n_dims] transducer positions
%                in grid indices (broadcast to all transducers if 1-row)
%   focus_pos  - [transducer_N x n_dims] or [1 x n_dims] focus positions in grid indices
%
% Output:
%   source        - struct with fields source.p_mask (logical grid) and source.p
%   source_labels - grid array of integer element labels (0 = inactive)
%   tr_arr        - transducer struct array updated with grid-point element dimensions
%                   and phase/amplitude fields
%
% See also: SOURCE_SENSOR_SETUP, CREATE_MATRIX_KARRAY, CREATE_CLOVER_ARRAY

arguments
    parameters (1,1) struct
    kgrid      (1,1)
    trans_pos  (:,:) {mustBeNumeric}
    focus_pos  (:,:) {mustBeNumeric}
end

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

    tr_arr = parameters.transducer;

    for it = 1:transducer_N
        tr = tr_arr(it);

        if strcmp(tr.type, 'annular')
            tr.annular.Elements_OD = 2 * floor(tr.annular.elem_od_mm / parameters.grid.resolution_mm / 2) + 1;
            tr.annular.Elements_ID = 2 * floor(tr.annular.elem_id_mm / parameters.grid.resolution_mm / 2) + 1;
            tr.annular.Elements_ID(tr.annular.elem_id_mm == 0) = 0;
			
        end
        
        tr.(tr.type).radius_grid = round(tr.(tr.type).curv_radius_mm / parameters.grid.resolution_mm);

        if it == 1
            % initialise struct array with full field set of tr
            tr_arr = repmat(tr, 1, transducer_N);
        end

        tr_arr(it) = tr;
    end

    source = struct();

    %% Branch 1: Custom element geometry (non-kWaveArray setup, multi-transducer)

    if parameters.grid.use_kWaveArray == 0

        % --- per-transducer CW signals (per-element) ---

        freqs = [tr_arr.freq_hz];
        if numel(unique(freqs)) > 1
            warn('Multiple source frequencies not yet supported. Using %i Hz.', freqs(1));
        end

        cw_signals_all    = cell(1, transducer_N);
        n_elements_per_T  = zeros(1, transducer_N);

        for it = 1:transducer_N
            tr = tr_arr(it);

            % elem_amp and elem_phase_rad are per-element vectors
            cw_signals_all{it} = createCWSignals( ...
                kgrid.t_array, ...
                freqs(1), ...               % use canonical frequency
                tr.(tr.type).elem_amp, ...
                tr.(tr.type).elem_phase_rad);       % [elem_n x Nt]

            n_elements_per_T(it) = tr.(tr.type).elem_n;
        end

        % global element indexing: each (transducer, element) → unique index
        offsets          = [0, cumsum(n_elements_per_T(1:end-1))];  % length transducer_N
        n_elements_total = sum(n_elements_per_T);

        % --- geometry: per-element bowls, numeric source_labels = global element index ---

        grid_dims       = parameters.grid.dims;
        transducer_mask = false(grid_dims);
        source_labels   = zeros(grid_dims);
        
        for it = 1:transducer_N
            tr          = tr_arr(it);
            trans_pos_i = tr_arr(it).trans_pos;
            focus_pos_i = tr_arr(it).focus_pos;
        
            for el_i = 1:tr.(tr.type).elem_n
                global_el = offsets(it) + el_i;  % 1..n_elements_total
        
                % outer element aperture
                if numel(parameters.grid.dims) == 3
                    bowl = makeBowl(grid_dims, trans_pos_i, tr.(tr.type).radius_grid, ...
                                    tr.(tr.type).Elements_OD(el_i), focus_pos_i);
                else
                    bowl = makeArc(grid_dims, trans_pos_i, tr.(tr.type).radius_grid, ...
                                   tr.(tr.type).Elements_OD(el_i), focus_pos_i);
                end

                % subtract inner aperture if applicable
                if tr.(tr.type).Elements_ID(el_i) > 0
                    if numel(parameters.grid.dims) == 3
                        bowl = bowl - makeBowl(grid_dims, trans_pos_i, ...
                                               tr.(tr.type).radius_grid, tr.(tr.type).Elements_ID(el_i), focus_pos_i);
                    else
                        bowl = bowl - makeArc(grid_dims, trans_pos_i, ...
                                              tr.(tr.type).radius_grid, tr.(tr.type).Elements_ID(el_i), focus_pos_i);
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

        % Note: Both `parameters` and `tr` are passed to this function to
        % allow support for multiple transducers in future implementations.
        tr = tr_arr(1);  % use first (and only) transducer here

        % 3D/2D positions for kWaveArray are taken from the first transducer
        trans_pos_1  = tr.trans_pos;
        focus_pos_1  = tr.focus_pos;

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

        switch tr.type
            case 'annular'
                % Set focus position and transducer position vectors in physical coordinates
                if numel(parameters.grid.dims) == 3
                    % 3D annular array
                    pos_vec   = [kgrid.x_vec(trans_pos_1(1)), kgrid.y_vec(trans_pos_1(2)), kgrid.z_vec(trans_pos_1(3))];
                    focus_vec = [kgrid.x_vec(focus_pos_1(1)), kgrid.y_vec(focus_pos_1(2)), kgrid.z_vec(focus_pos_1(3))];

                    karray.addAnnularArray(pos_vec, ...
                                           tr.annular.curv_radius_mm * 1e-3, ...
                                           [tr.annular.elem_id_mm; tr.annular.elem_od_mm] * 1e-3, ...
                                           focus_vec);

                elseif numel(parameters.grid.dims) == 2 && axisymmetric == false

                    % 2D arc-shaped element
                    pos_vec   = [kgrid.x_vec(trans_pos_1(1)), kgrid.y_vec(trans_pos_1(2))];
                    focus_vec = [kgrid.x_vec(focus_pos_1(1)), kgrid.y_vec(focus_pos_1(2))];

                    karray.addArcElement(pos_vec, ...
                                         tr.annular.curv_radius_mm * 1e-3, ...
                                         tr.annular.elem_od_mm * 1e-3, ...
                                         focus_vec);
                end
            case 'matrix'
                trans_pos_m = [kgrid.x_vec(tr.trans_pos(1));
                    kgrid.y_vec(tr.trans_pos(2));
                    kgrid.z_vec(tr.trans_pos(3))];

                focus_pos_m = [kgrid.x_vec(tr.focus_pos(1));
                    kgrid.y_vec(tr.focus_pos(2));
                    kgrid.z_vec(tr.focus_pos(3))];

                switch tr.matrix.matrix_shape.type
                    case 'define_here'
                        [elem_pos_m, tr] = convert_to_element_pos(parameters, tr, trans_pos_m, focus_pos_m);
                    case 'extract_from_file'               
                        elem_pos_m = extract_element_pos(parameters, tr,  trans_pos_m);
                    otherwise 
                        error('Matrix shape %s is unknown or not implemented.', tr.matrix.matrix_shape.type)
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
                    fig_filename = fullfile(parameters.io.dir_debug_source, ...
                        sprintf('sub-%03d_%s_transducer_element_distribution%s.fig', ...
                        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

                    png_filename = fullfile(parameters.io.dir_debug_source, ...
                        sprintf('sub-%03d_%s_transducer_element_distribution%s.png', ...
                        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

                    % Save outputs
                    saveas(h, fig_filename, 'fig')
                    saveas(h, png_filename, 'png')
                    close(h)
                end
                
                [karray, tr] = create_matrix_karray(kgrid, karray, parameters, tr, elem_pos_m, trans_pos, focus_pos);

            otherwise
                error('Array shape %s is unknown or not implemented.', tr.type)
        end

        % CW per-element for this transducer
        cw_signal = createCWSignals( ...
            kgrid.t_array, ...
            tr.freq_hz, ...
            tr.(tr.type).elem_amp, ...
            tr.(tr.type).elem_phase_rad);   % [elem_n x Nt]

        if axisymmetric == true && strcmp(tr.type, 'annular')

            kgrid_mirrored = kWaveGrid(kgrid.Nx, kgrid.dx, 2*kgrid.Ny - 1, kgrid.dy);
            karray_full = kWaveArray('Axisymmetric', false, 'BLITolerance', 0.01, 'UpsamplingRate', 100);

            x_offset     = (trans_pos_1(1)-1) * (1/parameters.grid.resolution_mm);
            position_base = [kgrid.x_vec(1) + x_offset*kgrid.dx, 0+eps];
            focus_pos_full = [0, 0+eps];
            
            for el_i = 1:tr.annular.elem_n
                el_OD_m = tr.annular.elem_od_mm(el_i) * 1e-3;
                y_shift = (el_i - (tr.annular.elem_n+1)/2) * el_OD_m;
                element_pos = position_base + [0, y_shift];

                karray_full.addArcElement(element_pos, tr.annular.curv_radius_mm * 1e-3, el_OD_m, focus_pos_full);
            end

            source_mask_full = false(kgrid_mirrored.Nx, kgrid_mirrored.Ny);
            grid_weights_3d  = zeros(tr.(tr.type).elem_n, kgrid_mirrored.Nx, kgrid_mirrored.Ny);

            for el_i = 1:tr.(tr.type).elem_n
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

            for el_i = 1:tr.(tr.type).elem_n
                weights_el = squeeze(grid_weights_cropped(el_i, :, :));
                weights_vec = weights_el(source_idx);
                source.p = source.p + bsxfun(@times, weights_vec, cw_signal(el_i, :));
            end

            % label sources at 50% weights
            source_labels = sum(grid_weights_cropped > 0.5, 1) > 0;
            source_labels = squeeze(source_labels);

        else
            % Non-axisymmetric: compute weights and distribute signals manually
            grid_weights_4d = zeros(tr.(tr.type).elem_n, kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));

            for ind = 1:tr.(tr.type).elem_n
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

            for ind = 1:tr.(tr.type).elem_n
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

            if parameters.simulation.debug == 1
                show_binary_mask_transducer(karray, kgrid, parameters);
            end
        end
        tr_arr(1) = tr;
    end
end
