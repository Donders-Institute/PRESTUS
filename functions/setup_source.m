function [source, source_labels, transducer_pars] = setup_source(parameters, kgrid, trans_pos, focus_pos)

% SETUP_SOURCE Creates a source for k-Wave simulations based on transducer parameters.
%
% This function generates a source structure for k-Wave simulations using the provided
% transducer parameters, computational grid, and transducer/focus positions. It supports
% multiple transducers for custom element geometries (use_kWaveArray == 0), and a single
% transducer for kWaveArray-based setups (use_kWaveArray ~= 0).
%
% Input:
%   parameters      - Struct containing simulation parameters (must contain parameters.transducers)
%   kgrid           - k-Wave grid struct (e.g., kWaveGrid)
%   trans_pos       - [nT x n_sim_dims] array of transducer positions in grid coordinates,
%                     or single [1 x n_sim_dims] / [n_sim_dims x 1] position broadcast to all nT
%   focus_pos       - same shape rules as trans_pos
%
% Output:
%   source          - Struct with source.p_mask and source.p (size depends on branch)
%   source_labels   - Grid of integer labels identifying active source regions
%   transducer_pars - Struct array of transducer parameters with extra grid-based fields
%
% Notes:
%   - For multiple transducers, all source_freq_hz are currently assumed equal; if they differ,
%     the frequency of the first transducer is used.
%   - Multiple transducers with kWaveArray (use_kWaveArray ~= 0) are not implemented.

    nT = numel(parameters.transducers);

    if parameters.use_kWaveArray ~= 0 && nT > 1
    error(['Multiple transducers with kWaveArray (use_kWaveArray ~= 0) ' ...
           'are not implemented yet. Use use_kWaveArray == 0 or a single transducer.']);
    end

    % Validate input dimensions and broadcast positions

    if isequal(size(trans_pos), [parameters.n_sim_dims 1]); trans_pos = trans_pos'; end
    if isequal(size(focus_pos), [parameters.n_sim_dims 1]); focus_pos = focus_pos'; end
    
    if size(trans_pos,1) == 1; trans_pos = repmat(trans_pos, nT, 1); end
    if size(focus_pos,1) == 1; focus_pos = repmat(focus_pos, nT, 1); end
    
    if size(trans_pos,1) ~= nT || size(focus_pos,1) ~= nT
        error('trans_pos and focus_pos must have one row per transducer.')
    end
    if size(trans_pos,2) ~= parameters.n_sim_dims || size(focus_pos,2) ~= parameters.n_sim_dims
        error('trans_pos and focus_pos columns must match parameters.n_sim_dims.')
    end

    % Convert element diameters from mm to grid points
    transducer_pars = parameters.transducers;
    grid_step_mm = parameters.grid_step_mm;

    for it = 1:nT
        tp = transducer_pars(it);
    
        tp.Elements_OD = 2 * floor(tp.Elements_OD_mm / grid_step_mm / 2) + 1;
        tp.Elements_ID = 2 * floor(tp.Elements_ID_mm / grid_step_mm / 2) + 1;
        tp.Elements_ID(tp.Elements_ID_mm == 0) = 0;
    
        tp.radius_grid = round(tp.curv_radius_mm / grid_step_mm);
    
        transducer_pars(it) = tp;
    end

    %% Setup source signal
    
    freqs  = [transducer_pars.source_freq_hz];
    amps   = [transducer_pars.source_amp];
    phases = [transducer_pars.source_phase_rad];
    
    if numel(unique(freqs)) > 1
        warning('Multiple source frequencies not yet supported. Using %i Hz.', freqs(1));
    end
    
    cw_signal = createCWSignals(kgrid.t_array, freqs(1), amps, phases);   % nT x Nt
    
    source = struct();

    %% Custom element geometry (non-kWaveArray setup)
    if parameters.use_kWaveArray == 0
        grid_dims = parameters.grid_dims;
        transducer_mask = zeros(grid_dims);
        source_labels = zeros(grid_dims);
    
        % one label (1..nT) and one row in cw_signal per transducer
        for it = 1:nT
            tp = transducer_pars(it);
    
            if parameters.n_sim_dims == 3
                bowl = makeBowl(grid_dims, trans_pos(it, :), tp.radius_grid, ...
                                tp.Elements_OD, focus_pos(it, :));
            else
                bowl = makeArc(grid_dims, trans_pos(it, :), tp.radius_grid, ...
                               tp.Elements_OD, focus_pos(it, :));
            end
    
            % inner bowl if needed
            if tp.Elements_ID > 0
                if length(grid_dims) == 3
                    bowl = bowl - makeBowl(grid_dims, trans_pos(it, :), ...
                                           tp.radius_grid, tp.Elements_ID, focus_pos(it, :));
                else
                    bowl = bowl - makeArc(grid_dims, trans_pos(it, :), ...
                                          tp.radius_grid, tp.Elements_ID, focus_pos(it, :));
                end
            end
    
            transducer_mask = transducer_mask + bowl;
            source_labels = source_labels + it * bowl;   % it = transducer index
        end
    
        % assign pressure signals: one row per transducer index
        p_mask_source_p = reshape(source_labels(source_labels > 0), [], 1);
        Nt = size(cw_signal, 2);
        source.p = zeros(length(p_mask_source_p), Nt);
    
        for ii = 1:length(p_mask_source_p)
            source.p(ii, :) = cw_signal(p_mask_source_p(ii), :);
        end
    
        source.p_mask = transducer_mask;

    %% kWaveArray-based setup (for advanced simulations)
    else
        disp('Setting up kWaveArray (might take a bit of time)');

        tp = transducer_pars(1);  % use first (and only) transducer here

        % Determine if axisymmetric mode should be enabled
        if parameters.n_sim_dims == 2 && ...
                isfield(parameters, 'axisymmetric') && ...
                parameters.axisymmetric == 1
            axisymmetric = true;
            disp("Using axisymmetric setup for 2D input...")
        else
            axisymmetric = false;
        end
        
        % Create kWaveArray object (always without axisymmetry)
        karray = kWaveArray('Axisymmetric', false, ...
                           'BLITolerance', 0.1, ...
                           'UpsamplingRate', 10, ...
                           'BLIType', 'sinc');
        
        % Set focus position and transducer position vectors in physical coordinates
        if parameters.n_sim_dims == 3
            % 3D annular array
            pos_vec = [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))];
            focus_vec = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))];
            
            karray.addAnnularArray(pos_vec, ...
                                   tp.curv_radius_mm * 1e-3, ...
                                   [tp.Elements_ID_mm; tp.Elements_OD_mm] * 1e-3, ...
                                   focus_vec);
        elseif parameters.n_sim_dims == 2 && axisymmetric == false

            % 2D arc-shaped element
            pos_vec = [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2))];
            focus_vec = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2))];
            
            karray.addArcElement(pos_vec, ...
                                 tp.curv_radius_mm * 1e-3, ...
                                 tp.Elements_OD_mm * 1e-3, ...
                                 focus_vec);
        end
        
        if axisymmetric
            % Create a mirrored grid with odd number of points
            kgrid_mirrored = kWaveGrid(kgrid.Nx, kgrid.dx, 2*kgrid.Ny - 1, kgrid.dy);
            karray_full = kWaveArray('Axisymmetric', false, 'BLITolerance', 0.01, 'UpsamplingRate', 100);
            
            % Base position in physical units
            x_offset = (trans_pos(1)-1)*(1/parameters.grid_step_mm);
            position_base = [kgrid.x_vec(1) + x_offset*kgrid.dx, 0+eps];
            focus_pos_full = [0, 0+eps];
            
            % Add each element separately with lateral displacement
            for el_i = 1:tp.n_elements
                % Element-specific parameters
                el_OD_m = tp.Elements_OD_mm(el_i) * 1e-3; % Outer diameter converted from millimeters to meters
                
                % Compute lateral y-position shift to arrange elements side by side,
                % centered around the base position. The formula centers the elements so that
                % element indices run from negative to positive offsets.
                % The spacing is set as the element diameter.
                y_shift = (el_i - (tp.n_elements+1)/2) * el_OD_m;
                
                % Calculate the element's exact position by shifting the base position in y.
                % position_base is the central position in x; y is offset by y_shift.
                element_pos = position_base + [0, y_shift];
                
                % Add the arc element to the kWaveArray with:
                % - element_pos: position of the element center on the rear surface [x, y] in meters
                % - tp.curv_radius_mm * 1e-3: radius of curvature in meters
                % - el_OD_m: aperture diameter in meters
                % - focus_pos_full: a point on the beam axis defining its angular direction [x, y] (meters)
                karray_full.addArcElement(element_pos, tp.curv_radius_mm * 1e-3, el_OD_m, focus_pos_full);
            end
    
            % Initialize mask and weights
            source_mask_full = false(kgrid_mirrored.Nx, kgrid_mirrored.Ny);
            grid_weights_3d = zeros(tp.n_elements, kgrid_mirrored.Nx, kgrid_mirrored.Ny);
            
            % Compute weights and aggregate mask
            for el_i = 1:tp.n_elements
                % Retrieve element-specific grid weights
                grid_weights_3d(el_i, :, :) = karray_full.getElementGridWeights(kgrid_mirrored, el_i);
                % Update the mask
                source_mask_full = source_mask_full | (squeeze(grid_weights_3d(el_i, :, :)) > 0);
            end
            
            % Crop to half space for axisymmetry
            source.p_mask = source_mask_full(:, kgrid.Ny:end);
            grid_weights_cropped = grid_weights_3d(:, :, kgrid.Ny:end);
            
            % Number of source points
            source_idx = find(source.p_mask);
            num_points = numel(source_idx);
            nTime = size(cw_signal, 2);
            
            % Initialize pressure matrix
            source.p = zeros(num_points, nTime);
            
            % Sum weighted signals over all elements
            for el_i = 1:tp.n_elements
                weights_el = squeeze(grid_weights_cropped(el_i, :, :));
                weights_vec = weights_el(source_idx);
                source.p = source.p + bsxfun(@times, weights_vec, cw_signal(el_i, :));
            end
            
            % Optional: create source_labels (thresholded weights)
            source_labels = sum(grid_weights_cropped > 0, 1) > 0;
            source_labels = squeeze(source_labels);
            
        else
            % --- Non-axisymmetric mode: compute weights and distribute signals manually ---
            
            % Initialize 4D array for element weights: [elements, Nx, Ny, Nz]
            grid_weights_4d = zeros(tp.n_elements, kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));
            
            for ind = 1:tp.n_elements
                fprintf('Computing weights for element %i...', ind);
                grid_weights_4d(ind,:,:,:) = karray.getElementGridWeights(kgrid, ind);
                fprintf(' done\n');
            end
            
            % Sum weights over elements and create binary mask of active source points
            binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;
            
            Nt = size(cw_signal, 2);
            mask_ind = find(binary_mask);
            num_source_points = sum(binary_mask(:));
            
            distributed_source_signal = zeros(num_source_points, Nt);
            source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz,1));
            
            if canUseGPU
                cw_signal = gpuArray(cw_signal);
                distributed_source_signal = gpuArray(distributed_source_signal);
            end
            
            % Distribute signals weighted by element weights
            for ind = 1:tp.n_elements
                source_weights = squeeze(grid_weights_4d(ind,:,:,:));
                el_binary_mask = source_weights ~= 0;
                
                if canUseGPU
                    source_weights = gpuArray(source_weights);
                end
                
                element_mask_ind = find(el_binary_mask);
                local_ind = ismember(mask_ind, element_mask_ind);
                
                distributed_source_signal(local_ind,:) = ...
                    distributed_source_signal(local_ind,:) + ...
                    source_weights(element_mask_ind) * cw_signal(ind,:);
                
                source_labels = source_labels + ind * el_binary_mask;
            end
            
            % Assign outputs
            source.p_mask = binary_mask;
            source.p = distributed_source_signal;
        end


    end

end
