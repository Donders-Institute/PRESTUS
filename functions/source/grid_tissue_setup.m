function [parameters, medium_masks, segmentation, skull_edge, planimg] = grid_tissue_setup(parameters)
% GRID_TISSUE_SETUP sets up the computational grid and tissue masks
% for ultrasound neuromodulation simulations, either by preprocessing a 
% subject-specific head model (T1-weighted MRI) or loading phantom/alternative grids.

if contains(parameters.simulation_medium, {'skull'; 'layered'})
    % preprocessing 
    % [multi-transducer] based on first transducer–focus pair in T1 space
    % Note: do not pass converted positions directly, handle in later function
    [medium_masks, segmentation, skull_edge, ~, ~, ...
        planimg.t1_image_orig, planimg.t1_header, planimg.transf, planimg.inv_transf] = ...
        preproc_head(parameters);
    
    if isempty(medium_masks)
        filename_output_table = '';
        return;
    end
    parameters.grid_dims  = size(medium_masks);
    parameters.n_sim_dims = numel(parameters.grid_dims);

else
    % In case simulations are not run in a skull of layered tissue, 
    % alternative grid dimensions are set up
    assert(isfield(parameters, 'default_grid_dims'), ...
        'The parameters structure should have the field grid_dims for the grid dimensions')
    parameters.grid_dims = parameters.default_grid_dims;
    if any(parameters.grid_dims==1)||length(parameters.grid_dims)==2
        parameters.grid_dims = squeeze(parameters.grid_dims);
        parameters.n_sim_dims = length(parameters.grid_dims);
        disp('One of the simulation grid dimensions is of length 1. Assuming you want 2D simulations ... dropping this dimension')
    end

    % read in phantoms directly as medium masks 
    if strcmp(parameters.simulation_medium, 'phantom')
        segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', parameters.subject_id));
        filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');
        segmented_img = niftiread(filename_segmented);
        if size(segmented_img) == parameters.grid_dims
            disp('Check passed: phantom dimensions fit requested grid...');
        elseif size(segmented_img') == parameters.grid_dims
            segmented_img = segmented_img';
            disp('Check passed: phantom dimensions fit requested grid...');
        else
            warning('WARNING: phantom dimensions DO NOT fit requested grid...')
        end
        % create medium mask according to indices in parameters.layer_labels (see preproc_smooth_and_crop.m)
        [medium_masks] = preproc_medium_mask(segmented_img, parameters, 1);
        segmentation = segmented_img; clear segmented_img;
    else
        medium_masks = [];
        segmentation = zeros(parameters.grid_dims);
    end

    % specify that no transformation was applied
    planimg.t1_image_orig = [];
    planimg.t1_header = [];
    planimg.transf = zeros(parameters.grid_dims);
    planimg.inv_transf = zeros(parameters.grid_dims);

    % create empty outputs for other variables
    skull_edge = zeros(parameters.grid_dims);

end