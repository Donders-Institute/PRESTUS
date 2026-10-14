function [source, source_labels, tr_arr] = source_create(parameters, kgrid, trans_pos, focus_pos)
% SOURCE_CREATE  Build a k-Wave source struct from PRESTUS transducer parameters
%
% Dispatches to one of two branches:
%   Branch 1 (use_kWaveArray == 0): Supports multiple transducers. Each
%     annular element is placed with makeBowl / makeArc. CW signals are
%     built with createCWSignals and assigned row-by-row to source.p.
%   Branch 2 (use_kWaveArray ~= 0): kWaveArray handles element geometry;
%     grid weights distribute signals. Multiple transducers supported in
%     non-axisymmetric mode; all transducers are added to one kWaveArray
%     with global element indexing. Axisymmetric mode is single-transducer only.
%
%   Axisymmetric mode (Branch 2, annular type, grid.axisymmetric == 1):
%     kWaveArray does not offer native axisymmetric source weighting, so a
%     mirrored-grid approximation is used:
%       1. A bilateral grid with 2*kgrid.Ny-1 columns is created (odd size,
%          so the centre column is exactly at y=0).
%       2. Each element is placed on the bilateral grid according to its
%          radial position: central disks (inner_diam == 0) are placed at
%          y=0 with diameter = outer_diam; annular rings (inner_diam > 0)
%          are the segment of the bowl's arc between r_inner and r_outer:
%          midpoint on the bowl surface at the mid polar angle (set back
%          from the apex plane by the sagitta), diameter = chord between the
%          ring edges. All elements are aimed at the bowl's centre of
%          curvature (curv_radius_mm along trans_pos -> focus_pos).
%       3. Element weights are computed on the full bilateral grid, then
%          cropped to columns kgrid.Ny:end → [Nx x Ny] where column 1 = r=0.
%
%   Flat annular elements (curv_radius_mm = Inf, e.g. a plane piston): a line
%     element (addLineElement) across the axis in axisymmetric and 2D mode, a
%     disc (addDiscElement) in 3D; the face lies at trans_pos, normal to the
%     trans_pos -> focus_pos direction. kWaveArray only (Branch 2); in 3D a
%     single element without an inner diameter.
%
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

        for it = 1:transducer_N
            if strcmp(tr_arr(it).type, 'annular') && any(isinf(tr_arr(it).annular.curv_radius_mm))
                error('PRESTUS:sourceCreate:flatNeedsKWaveArray', ...
                    ['Transducer %d is flat (curv_radius_mm = Inf), which needs the kWaveArray source ' ...
                     '(grid.use_kWaveArray = 1). makeBowl/makeArc only build curved elements.'], it);
            end
        end

        % --- per-transducer CW signals (per-element) ---

        freqs = [tr_arr.freq_hz];
        cw_signals_all    = cell(1, transducer_N);
        n_elements_per_T  = zeros(1, transducer_N);

        for it = 1:transducer_N
            tr = tr_arr(it);

            % elem_p_amp and elem_phase_rad are per-element vectors
            cw_signals_all{it} = createCWSignals( ...
                kgrid.t_array, ...
                freqs(it), ...              % per-transducer frequency
                tr.(tr.type).elem_p_amp, ...
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
            trans_pos_i = trans_pos(it,:);
            focus_pos_i = focus_pos(it,:);
        
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

    %% Branch 2: kWaveArray-based setup (multi-transducer supported in non-axisymmetric mode)

    else
        disp('Setting up kWaveArray (might take a bit of time)');

        % Determine if axisymmetric mode is requested.
        % Axisymmetric mode uses a mirrored-grid approximation and supports
        % only a single annular transducer (bilateral grid cannot represent
        % two physical transducers on the radial half-plane).
        if numel(parameters.grid.dims) == 2 && ...
                isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
            axisymmetric = true;
            disp("Using axisymmetric setup for 2D input...");
            if transducer_N > 1
                error('PRESTUS:sourceCreate:axisymMultiTxNotSupported', ...
                    ['Axisymmetric kWaveArray mode supports only a single transducer. ' ...
                     'Use use_kWaveArray == 0 or disable grid.axisymmetric for multiple transducers.']);
            end
            if ~strcmp(tr_arr(1).type, 'annular')
                error('PRESTUS:sourceCreate:axisymTypeNotSupported', ...
                    ['Axisymmetric mode only supports annular transducers. ' ...
                     'Transducer type ''%s'' requires a full 2D or 3D grid. ' ...
                     'Disable grid.axisymmetric or use a non-axisymmetric setup.'], tr_arr(1).type);
            end
        else
            axisymmetric = false;
        end

        % Band-limited interpolation of the source: grid.bli_tolerance and
        % grid.bli_upsampling_rate, else 0.1 / 10 (3D) and 0.01 / 100 (axisymmetric)
        if axisymmetric; bli = [0.01, 100]; else; bli = [0.1, 10]; end
        if isfield(parameters.grid, 'bli_tolerance') && ~isempty(parameters.grid.bli_tolerance)
            bli(1) = parameters.grid.bli_tolerance;
        end
        if isfield(parameters.grid, 'bli_upsampling_rate') && ~isempty(parameters.grid.bli_upsampling_rate)
            bli(2) = parameters.grid.bli_upsampling_rate;
        end

        % Create kWaveArray object (always without native axisymmetry)
        karray = kWaveArray('Axisymmetric', false, ...
                           'BLITolerance', bli(1), ...
                           'UpsamplingRate', bli(2), ...
                           'BLIType', 'sinc', ...
                           'SinglePrecision', strcmp(parameters.simulation.precision, 'single'));

        if axisymmetric

            % ---- Axisymmetric single-transducer path ----
            % Uses mirrored bilateral grid then crops to right half-plane.
            % See Branch 2 docstring in source_create for derivation details.

            tr          = tr_arr(1);
            trans_pos_1 = trans_pos(1,:);
            focus_pos_1 = focus_pos(1,:);

            % Add arc to main karray (not used below but keeps object consistent)
            pos_vec   = [kgrid.x_vec(trans_pos_1(1)), kgrid.y_vec(trans_pos_1(2))];
            focus_vec = [kgrid.x_vec(focus_pos_1(1)), kgrid.y_vec(focus_pos_1(2))];
            is_flat   = isinf(tr.annular.curv_radius_mm);
            if ~is_flat
                karray.addArcElement(pos_vec, tr.annular.curv_radius_mm * 1e-3, max(tr.annular.elem_od_mm) * 1e-3, focus_vec);   % full aperture (addArcElement takes one diameter)
            end

            cw_signal = createCWSignals(kgrid.t_array, tr.freq_hz, ...
                tr.(tr.type).elem_p_amp, tr.(tr.type).elem_phase_rad);

            kgrid_mirrored = kWaveGrid(kgrid.Nx, kgrid.dx, 2*kgrid.Ny - 1, kgrid.dy);
            karray_full = kWaveArray('Axisymmetric', false, 'BLITolerance', bli(1), 'UpsamplingRate', bli(2));

            x_trans_m      = kgrid.x_vec(trans_pos_1(1));
            ax_dir         = sign(kgrid.x_vec(focus_pos_1(1)) - x_trans_m);   % +1 / -1: axial direction to the focus
            curv_r_m       = tr.annular.curv_radius_mm * 1e-3;
            % Rings are aimed at the bowl's centre of curvature (geometric focus),
            % not at the (steered) focus: a ring off the axis aimed at the focus
            % is tilted away from the bowl surface
            geom_focus     = [x_trans_m + ax_dir * curv_r_m, 0+eps];

            for el_i = 1:tr.annular.elem_n
                inner_r_m = tr.annular.elem_id_mm(el_i) / 2 * 1e-3;
                outer_r_m = tr.annular.elem_od_mm(el_i) / 2 * 1e-3;
                if is_flat
                    % Flat face at the transducer plane, normal to the axis: a disc is
                    % the line r in [-R, R] of the mirrored grid, a ring r in [Ri, Ro]
                    % (only r >= 0 is kept by the crop below)
                    if inner_r_m == 0
                        karray_full.addLineElement([x_trans_m, -outer_r_m], [x_trans_m, outer_r_m]);
                    else
                        karray_full.addLineElement([x_trans_m, inner_r_m], [x_trans_m, outer_r_m]);
                    end
                    continue
                end
                if inner_r_m == 0
                    element_pos   = [x_trans_m, 0+eps];
                    el_diameter_m = 2 * outer_r_m;
                else
                    % Ring = the segment of the bowl's arc between polar angles
                    % asin(Ri/R) and asin(Ro/R) about the centre of curvature:
                    % midpoint on the bowl surface (set back from the apex plane by
                    % the sagitta, ~1-2 wavelengths for the outer CTX/IGT rings),
                    % diameter = chord between the ring edges
                    th_in         = asin(inner_r_m / curv_r_m);
                    th_out        = asin(outer_r_m / curv_r_m);
                    th_mid        = (th_in + th_out) / 2;
                    element_pos   = [x_trans_m + ax_dir * curv_r_m * (1 - cos(th_mid)), curv_r_m * sin(th_mid)];
                    el_diameter_m = 2 * curv_r_m * sin((th_out - th_in) / 2);
                end
                karray_full.addArcElement(element_pos, curv_r_m, el_diameter_m, geom_focus);
            end

            source_mask_full = false(kgrid_mirrored.Nx, kgrid_mirrored.Ny);
            grid_weights_3d  = zeros(tr.(tr.type).elem_n, kgrid_mirrored.Nx, kgrid_mirrored.Ny);

            for el_i = 1:tr.(tr.type).elem_n
                grid_weights_3d(el_i, :, :) = karray_full.getElementGridWeights(kgrid_mirrored, el_i);
                % Keep every grid point with a non-zero BLI weight, including the
                % negative sinc side lobes, as the 3D path (getElementBinaryMask)
                % and k-Wave do. A positive threshold (formerly > 0.25) drops part
                % of the band-limited source and raised the free-water focal
                % pressure by 6-11 % (ITRUSST BM1, 0.5-0.25 mm; see test_itrusst_benchmark).
                source_mask_full = source_mask_full | (squeeze(grid_weights_3d(el_i, :, :)) ~= 0);
            end

            source.p_mask        = source_mask_full(:, kgrid.Ny:end);
            grid_weights_cropped = grid_weights_3d(:, :, kgrid.Ny:end);

            source_idx = find(source.p_mask);
            num_points = numel(source_idx);
            nTime      = size(cw_signal, 2);

            if ismember(parameters.simulation.code_type, {'matlab_cpu', 'matlab_gpu'})
                % Lazy rank-2 CWSourceMatrix, as the non-axisymmetric path: the
                % dense num_points x Nt matrix reaches tens of GB at fine grids
                % with long runs (ITRUSST reference grid: 0.05 mm, up to 720k
                % steps). createCWSignals is A*sin(wt + phi) = Re(A e^{i(phi - pi/2)} e^{iwt}),
                % so the signal is that of the dense path below without its 4-period
                % start ramp (as in the non-axisymmetric lazy path); the steady state
                % is the same.
                C_src = zeros(num_points, 1);
                for el_i = 1:tr.(tr.type).elem_n
                    weights_el = squeeze(grid_weights_cropped(el_i, :, :));
                    C_src = C_src + weights_el(source_idx) * tr.annular.elem_p_amp(el_i) * ...
                        exp(1i * (tr.annular.elem_phase_rad(el_i) - pi/2));
                end
                source.p = CWSourceMatrix(real(C_src), imag(C_src), 2 * pi * tr.freq_hz, kgrid.dt, nTime);
                fprintf('[source_create] axisymmetric %s: using lazy CWSourceMatrix (%d src x %d steps, %.1f MB saved)\n', ...
                    parameters.simulation.code_type, num_points, nTime, num_points * nTime * 8 / 1e6);
            else
                source.p = zeros(num_points, nTime);
                for el_i = 1:tr.(tr.type).elem_n
                    weights_el  = squeeze(grid_weights_cropped(el_i, :, :));
                    weights_vec = weights_el(source_idx);
                    source.p    = source.p + bsxfun(@times, weights_vec, cw_signal(el_i, :));
                end
            end

            source_labels = sum(grid_weights_cropped > 0.5, 1) > 0;
            source_labels = squeeze(source_labels);

            tr_arr(1) = tr;

        else

            % ---- Non-axisymmetric path: supports N transducers ----
            % Elements from all transducers are added to one karray; a global
            % element index (offset per transducer) maps rows of the combined
            % CW signal matrix to karray element slots.

            elem_counts    = zeros(1, transducer_N);
            cw_signals_all = cell(1, transducer_N);

            for it = 1:transducer_N
                tr           = tr_arr(it);
                trans_pos_it = trans_pos(it,:);
                focus_pos_it = focus_pos(it,:);

                switch tr.type
                    case 'annular'
                        is_flat = isinf(tr.annular.curv_radius_mm);
                        if is_flat && numel(parameters.grid.dims) == 3
                            % Flat disc (plane piston), normal along trans_pos -> focus_pos
                            if tr.annular.elem_n ~= 1 || tr.annular.elem_id_mm(1) ~= 0
                                error('PRESTUS:sourceCreate:flatRings3D', ...
                                    ['Transducer %d: a flat (curv_radius_mm = Inf) annular source in 3D must be a ' ...
                                     'single disc (elem_n = 1, elem_id_mm = 0); kWaveArray has no flat ring element.'], it);
                            end
                            pos_vec   = [kgrid.x_vec(trans_pos_it(1)), kgrid.y_vec(trans_pos_it(2)), kgrid.z_vec(trans_pos_it(3))];
                            focus_vec = [kgrid.x_vec(focus_pos_it(1)), kgrid.y_vec(focus_pos_it(2)), kgrid.z_vec(focus_pos_it(3))];
                            karray.addDiscElement(pos_vec, tr.annular.elem_od_mm(1) * 1e-3, focus_vec);
                        elseif is_flat
                            % Flat element in 2D: a line through trans_pos normal to the beam
                            pos_vec   = [kgrid.x_vec(trans_pos_it(1)), kgrid.y_vec(trans_pos_it(2))];
                            focus_vec = [kgrid.x_vec(focus_pos_it(1)), kgrid.y_vec(focus_pos_it(2))];
                            beam = (focus_vec - pos_vec) / norm(focus_vec - pos_vec);
                            perp = [-beam(2), beam(1)];
                            for el_i = 1:tr.annular.elem_n
                                ri = tr.annular.elem_id_mm(el_i) / 2 * 1e-3;
                                ro = tr.annular.elem_od_mm(el_i) / 2 * 1e-3;
                                if ri == 0
                                    karray.addLineElement(pos_vec - ro * perp, pos_vec + ro * perp);
                                else
                                    error('PRESTUS:sourceCreate:flatRings2D', ...
                                        'Transducer %d: flat rings are supported in axisymmetric mode only.', it);
                                end
                            end
                        elseif numel(parameters.grid.dims) == 3
                            pos_vec   = [kgrid.x_vec(trans_pos_it(1)), kgrid.y_vec(trans_pos_it(2)), kgrid.z_vec(trans_pos_it(3))];
                            focus_vec = [kgrid.x_vec(focus_pos_it(1)), kgrid.y_vec(focus_pos_it(2)), kgrid.z_vec(focus_pos_it(3))];
                            karray.addAnnularArray(pos_vec, ...
                                                   tr.annular.curv_radius_mm * 1e-3, ...
                                                   [tr.annular.elem_id_mm; tr.annular.elem_od_mm] * 1e-3, ...
                                                   focus_vec);
                        elseif numel(parameters.grid.dims) == 2
                            pos_vec   = [kgrid.x_vec(trans_pos_it(1)), kgrid.y_vec(trans_pos_it(2))];
                            focus_vec = [kgrid.x_vec(focus_pos_it(1)), kgrid.y_vec(focus_pos_it(2))];
                            karray.addArcElement(pos_vec, ...
                                                 tr.annular.curv_radius_mm * 1e-3, ...
                                                 tr.annular.elem_od_mm * 1e-3, ...
                                                 focus_vec);
                        end

                    case 'matrix'
                        trans_pos_m = [kgrid.x_vec(trans_pos_it(1));
                                       kgrid.y_vec(trans_pos_it(2));
                                       kgrid.z_vec(trans_pos_it(3))];
                        focus_pos_m = [kgrid.x_vec(focus_pos_it(1));
                                       kgrid.y_vec(focus_pos_it(2));
                                       kgrid.z_vec(focus_pos_it(3))];

                        switch tr.matrix.matrix_shape.type
                            case 'define_here'
                                [elem_pos_m, tr] = convert_to_element_pos(parameters, tr, trans_pos_m, focus_pos_m);
                            case 'extract_from_file'
                                elem_pos_m = extract_element_pos(parameters, tr, trans_pos_m);
                            otherwise
                                error('Matrix shape %s is unknown or not implemented.', tr.matrix.matrix_shape.type)
                        end

                        if parameters.simulation.debug == 1
                            elem_pos_mm = elem_pos_m' * 1e3;
                            h = figure;
                            scatter3(elem_pos_mm(:,1), elem_pos_mm(:,2), elem_pos_mm(:,3), 60, 'filled');
                            axis equal; xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
                            view([0 90]); title('Transducer Element Distribution'); grid on;
                            src_trans_suffix = trans_suffix(it, transducer_N);
                            fig_filename = fullfile(parameters.io.dir_debug_source, ...
                                sprintf('sub-%03d_%s_transducer_element_distribution%s%s.fig', ...
                                parameters.subject_id, parameters.simulation.medium, src_trans_suffix, parameters.io.output_affix));
                            png_filename = fullfile(parameters.io.dir_debug_source, ...
                                sprintf('sub-%03d_%s_transducer_element_distribution%s%s.png', ...
                                parameters.subject_id, parameters.simulation.medium, src_trans_suffix, parameters.io.output_affix));
                            saveas(h, fig_filename, 'fig'); saveas(h, png_filename, 'png'); close(h);
                        end

                        [karray, tr] = create_matrix_karray(kgrid, karray, parameters, tr, elem_pos_m, trans_pos_it, focus_pos_it);

                    otherwise
                        error('Array shape %s is unknown or not implemented.', tr.type)
                end

                elem_counts(it)    = tr.(tr.type).elem_n;
                cw_signals_all{it} = createCWSignals(kgrid.t_array, tr.freq_hz, ...
                                                      tr.(tr.type).elem_p_amp, tr.(tr.type).elem_phase_rad);
                tr_arr(it) = tr;
            end

            total_elements = sum(elem_counts);
            offsets        = [0, cumsum(elem_counts(1:end-1))];

            % Two-pass weight computation: never store more than one element's
            % weight vector at a time (avoids total_elements x Ngrid peak allocation).
            Ngrid = kgrid.Nx * kgrid.Ny * max(kgrid.Nz, 1);
            Nt    = size(cw_signals_all{1}, 2);

            % Pass 1: build union binary mask using getElementBinaryMask, which
            % returns a logical array and skips the sinc BLI accumulation
            % entirely (MaskOnly=true path in offGridPoints).
            acc_mask = false(Ngrid, 1);
            for ind = 1:total_elements
                fprintf('Computing weights for element %i (pass 1)...', ind);
                w = karray.getElementBinaryMask(kgrid, ind);
                acc_mask = acc_mask | w(:);
                fprintf(' done\n');
            end

            binary_mask = reshape(acc_mask, kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
            mask_ind    = find(acc_mask);
            num_src     = numel(mask_ind);

            % Pass 2: build sparse W_src (num_src x total_elements) via COO triplets.
            % Pre-allocate conservatively; each element typically covers O(1000) mask pts.
            coo_i = zeros(num_src * total_elements, 1, 'double');  % over-alloc, trimmed below
            coo_j = zeros(num_src * total_elements, 1, 'double');
            coo_v = zeros(num_src * total_elements, 1, 'double');
            ptr   = 0;
            for ind = 1:total_elements
                fprintf('Computing weights for element %i (pass 2)...', ind);
                w    = karray.getElementGridWeights(kgrid, ind);
                col  = w(mask_ind);
                clear w;   % release 604 MB before next iteration
                nz   = find(col);
                nnz_ = numel(nz);
                if nnz_ > 0
                    coo_i(ptr+1:ptr+nnz_) = nz;
                    coo_j(ptr+1:ptr+nnz_) = ind;
                    coo_v(ptr+1:ptr+nnz_) = col(nz);
                    ptr = ptr + nnz_;
                end
                fprintf(' done\n');
            end
            W_src = sparse(coo_i(1:ptr), coo_j(1:ptr), coo_v(1:ptr), num_src, total_elements);
            clear coo_i coo_j coo_v;

            % For matlab_cpu and matlab_gpu: use lazy rank-2 CWSourceMatrix
            % instead of the full num_src x Nt dense matrix. Every CW voxel
            % signal is:
            %   p(i,t) = Re(C_i)*cos(omega*t) - Im(C_i)*sin(omega*t)
            % where C_i = W_src(i,:) * C_elem — a single complex scalar per
            % voxel. This saves O(num_src * Nt) RAM (typically tens–hundreds
            % of MB). k-Wave's dataCast passes through harmlessly:
            %   single(source.p)   → CWSourceMatrix.single() → obj unchanged
            %   gpuArray(source.p) → CWSourceMatrix.gpuArray() → Re_C/Im_C
            %                        uploaded to GPU; subsref then returns
            %                        gpuArray slices for the GPU time loop.
            % NOT valid for cpp backends (source.p written to HDF5 → needs
            % a dense matrix), or when transducers have different frequencies.
            % jwave_wrapper writes a CWSourceMatrix in factored form (code_type
            % does not apply to jwave), so jwave always takes the lazy path.
            is_jwave    = isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'jwave');
            use_cw_lazy = (ismember(parameters.simulation.code_type, {'matlab_cpu', 'matlab_gpu'}) || is_jwave) && ...
                          numel(unique([tr_arr.freq_hz])) == 1;

            if use_cw_lazy
                % Build per-element complex amplitude: A_j * exp(i*phi_j)
                C_elem = zeros(total_elements, 1);
                omega_cw = 2 * pi * tr_arr(1).freq_hz;
                for it = 1:transducer_N
                    r = offsets(it) + (1:elem_counts(it));
                    C_elem(r) = tr_arr(it).(tr_arr(it).type).elem_p_amp(:) .* ...
                                 exp(1i * tr_arr(it).(tr_arr(it).type).elem_phase_rad(:));
                end
                % Per-voxel complex amplitude: sparse matrix multiply (negligible cost).
                % full(): with a single element C_elem is a scalar and the product
                % stays sparse, which CWSourceMatrix cannot cast to single.
                C_src = full(W_src * C_elem);   % num_src x 1 complex
                distributed_source_signal = CWSourceMatrix( ...
                    real(C_src), imag(C_src), omega_cw, kgrid.dt, Nt);
                fprintf('[source_create] %s: using lazy CWSourceMatrix (%d src x %d steps, %.1f MB saved)\n', ...
                    parameters.simulation.code_type, num_src, Nt, num_src * Nt * 4 / 1e6);
            else
                % Stack all element CW signals: total_elements x Nt
                all_cw = zeros(total_elements, Nt, parameters.simulation.precision);
                for it = 1:transducer_N
                    r = offsets(it) + (1:elem_counts(it));
                    all_cw(r, :) = cw_signals_all{it};
                end
                % Sparse CPU multiply. Cast immediately to simulation precision:
                % W_src (double sparse) * all_cw (single) would upcast to double otherwise.
                distributed_source_signal = cast(W_src * double(all_cw), parameters.simulation.precision);
            end

            % source_labels: sum of element indices weighted by binary occupancy.
            labels_at_src = (W_src ~= 0) * (1:total_elements)';  % num_src x 1
            source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
            source_labels(mask_ind) = labels_at_src;

            source.p_mask = binary_mask;
            source.p      = distributed_source_signal;

            if parameters.simulation.debug == 1
                show_binary_mask_transducer(karray, kgrid, parameters);
            end
        end
    end

    fprintf('[source_create] ---- source summary ----\n');
    fprintf('[source_create]   transducers       : %d\n', transducer_N);
    for it = 1:transducer_N
        tr_s = tr_arr(it);
        elem_n_s = tr_s.(tr_s.type).elem_n;
        amps_s = tr_s.(tr_s.type).elem_p_amp;
        if numel(amps_s) == 1
            fprintf('[source_create]   transducer %d elements : %d | p_amp : %.2f Pa\n', it, elem_n_s, amps_s);
        else
            fprintf('[source_create]   transducer %d elements : %d | p_amp : %.2f – %.2f Pa\n', it, elem_n_s, min(amps_s), max(amps_s));
        end
    end
    fprintf('[source_create]   active voxels     : %d\n', nnz(source.p_mask));
    freqs_s = [tr_arr.freq_hz];
    if numel(unique(freqs_s)) == 1
        fprintf('[source_create]   frequency         : %.0f Hz (%.3f kHz)\n', freqs_s(1), freqs_s(1)/1e3);
    else
        for it = 1:transducer_N
            fprintf('[source_create]   freq transducer %d  : %.0f Hz (%.3f kHz)\n', it, freqs_s(it), freqs_s(it)/1e3);
        end
    end
    if exist('axisymmetric', 'var') && axisymmetric
        fprintf('[source_create]   mode              : axisymmetric (mirrored-grid approximation)\n');
    end
    fprintf('[source_create] ------------------------------------\n');
end
