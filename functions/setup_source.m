function [source, source_labels, transducer_pars] = setup_source(parameters, kgrid, trans_pos, focus_pos)

% Function that creates a source for kwave based on the parameters structure, kwave grid, the transducer coordinates
% in the computational grid, the geometric focus coordinates. Also returns
% a 3d array with different numbers indicating different transducer
% elements and the updated transducer parameters.
 

if ~isequal(size(focus_pos), size(trans_pos))
    error('Transducer and focus positions should be arrays of equal size')
end

if isequal(size(focus_pos), [parameters.n_sim_dims 1])
    focus_pos = focus_pos';
    trans_pos = trans_pos';
elseif ~isequal(size(focus_pos), [1 parameters.n_sim_dims])
	error('Transducer and focus positions should have the size [N 1] or [1 N] where N is equal to the number of simulation dimensions (parameters.n_sim_dims).')
end

% Convert from [mm] to [grid points]
% and round the diameter to the nearest odd integer
transducer_pars = parameters.transducer;
grid_step_mm = parameters.grid_step_mm;

transducer_pars.Elements_OD = 2*floor(transducer_pars.Elements_OD_mm / grid_step_mm / 2) + 1; % [grid points]
transducer_pars.Elements_ID = 2*floor(transducer_pars.Elements_ID_mm / grid_step_mm / 2) + 1; % [grid points]

transducer_pars.Elements_ID(transducer_pars.Elements_ID_mm==0) = 0;

% Radius of curvature (ROC) indirect from the focal length or direct
% transducer_pars.radius_mm = sqrt((transducer_pars.Elements_OD_mm(2)/2)^2 + (focal_dist_mm)^2);   % [mm] % indirect
% transducer_pars.focal_dist_mm = focal_dist_mm;
% transducer_pars.focal_dist_based_on_focal_pos_mm = norm(trans_pos - focus_pos)*grid_step_mm;

% Convert from [mm] to [grid points]
transducer_pars.radius_grid = round(transducer_pars.curv_radius_mm / grid_step_mm);  % [grid points]

% Setup source
% define the input signal for each element
cw_signal = createCWSignals(kgrid.t_array, parameters.transducer.source_freq_hz, parameters.transducer.source_amp, parameters.transducer.source_phase_rad);

source = struct();

if parameters.use_kWaveArray == 0
    grid_dims = parameters.grid_dims;
    transducer_mask = zeros(grid_dims);
    source_labels = zeros(grid_dims);

    % create element bowls one by one
    % for the elements after the first, the inner bowl is subtracted
    for el_i = 1:transducer_pars.n_elements
        if parameters.n_sim_dims == 3
            bowl = makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_OD(el_i), focus_pos);
        else
            bowl = makeArc(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_OD(el_i), focus_pos);
        end
        if transducer_pars.Elements_ID(el_i) > 0
            if length(grid_dims)==3
                bowl = bowl - makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_ID(el_i), focus_pos);
            else
                bowl = bowl - makeArc(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_ID(el_i), focus_pos);
            end
        end

        % Define the binary source mask
        transducer_mask = transducer_mask + bowl;

        % Label the elements for setting the source signal and to identify the elements in a figure
        source_labels = source_labels + el_i*bowl;
    end
    % ensure each source point has an assigned time series of the source signal
    p_mask_source_p = source_labels;
    p_mask_source_p(p_mask_source_p(:,:,:) == 0) = [];
    p_mask_source_p = reshape(p_mask_source_p,[],1);

    source.p = zeros(length(p_mask_source_p),length(cw_signal));
    
    for ii = 1 : length(p_mask_source_p)
        source.p(ii, :) = cw_signal(p_mask_source_p(ii), :);
    end
    source.p_mask = transducer_mask;

else
    disp('Setting up kWaveArray (might take a bit of time)');
    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.1, 'UpsamplingRate', 10, 'BLIType', 'sinc');

    % add bowl shaped element
    karray.addAnnularArray([kgrid.x_vec(trans_pos(1)) kgrid.y_vec(trans_pos(2)) kgrid.z_vec(trans_pos(3))], transducer_pars.curv_radius_mm*1e-3, [transducer_pars.Elements_ID_mm; transducer_pars.Elements_OD_mm]*1e-3, [kgrid.x_vec(focus_pos(1)) kgrid.y_vec(focus_pos(2)) kgrid.z_vec(focus_pos(3))])
	
    %karray.addBowlElement([kgrid.x_vec(trans_pos(1)) kgrid.y_vec(trans_pos(2)) kgrid.z_vec(trans_pos(3))], transducer_pars.curv_radius_mm*1e-3, transducer_pars.Elements_OD_mm*1e-3, [kgrid.x_vec(focus_pos(1)) kgrid.y_vec(focus_pos(2)) kgrid.z_vec(focus_pos(3))])
    grid_weights_4d = zeros(transducer_pars.n_elements, kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
    for ind = 1:transducer_pars.n_elements
        fprintf('Computing weights for element %i...', ind); 
        grid_weights_4d(ind,:,:,:) = karray.getElementGridWeights(kgrid, ind);                
        fprintf(' done\n');
    end
    % get the binary mask
    binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;
    % number of time points in the signal
    Nt = size(cw_signal, 2);

    mask_ind = find(binary_mask);
    num_source_points = sum(binary_mask(:));
     
    % bits of code from kWaveArray.m
    % initialise the source signal
    distributed_source_signal = zeros(num_source_points, Nt);
    
    source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
    if canUseGPU
        cw_signal = gpuArray(cw_signal);
        distributed_source_signal = gpuArray(distributed_source_signal);
    end

    % loop through the elements
    for ind = 1:transducer_pars.n_elements

        % get the offgrid source weights
        source_weights = squeeze(grid_weights_4d(ind,:,:,:));
        el_binary_mask = source_weights ~= 0;
        if canUseGPU
            source_weights = gpuArray(source_weights);
        end
        % get indices of the non-zero points 
        element_mask_ind = find(el_binary_mask);
        % convert these to indices in the distributed source
        local_ind = ismember(mask_ind, element_mask_ind);

        % add to distributed source
        distributed_source_signal(local_ind, :) = ...
            distributed_source_signal(local_ind, :) ...
            + source_weights(element_mask_ind) * cw_signal(ind, :);

        source_labels = source_labels+ind*el_binary_mask;

    end
            

    % assign binary mask
    source.p_mask = binary_mask;

    % assign source signals
    source.p = distributed_source_signal;
end
end
