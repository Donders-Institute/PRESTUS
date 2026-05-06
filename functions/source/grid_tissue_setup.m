function [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = grid_tissue_setup(parameters)
% GRID_TISSUE_SETUP  Set up computational grid and tissue masks for TUS simulation
%
% Dispatches to one of three setup paths depending on parameters.simulation.medium:
%   'layered'  - preprocesses a SimNIBS T1-weighted head segmentation
%   'phantom'  - loads a pre-built NIfTI phantom from the SimNIBS m2m folder
%   other      - initialises a homogeneous water grid of default_dims
%
% Use as:
%   [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = ...
%       grid_tissue_setup(parameters)
%
% Input:
%   parameters - (struct) PRESTUS config; must contain:
%                  simulation.medium (char)
%                  grid.default_dims [1x2] or [1x3] voxels
%                  path.seg (char) root path to SimNIBS m2m folders
%                  subject_id (int)
%
% Output:
%   parameters   - updated: grid.dims set from loaded/created mask
%   medium_masks - [Nx x Ny (x Nz)] uint8, medium layer label map
%   segmentation - [Nx x Ny (x Nz)] uint8, tissue label map
%   bone_mask    - [Nx x Ny (x Nz)] logical, binary skull mask (always present)
%   pseudoCT     - [Nx x Ny (x Nz)] numeric, Hounsfield-unit skull image
%                  ([] when parameters.pct.enabled ~= 1 or medium is not layered)
%   planimg      - (struct) planning image data (t1_image_orig, t1_header,
%                    transf, inv_transf); fields are empty for non-layered media
%
% See also: GRID_AXISYMMETRY, PREPROC_HEAD, PREPROC_MEDIUM_MASK

arguments
    parameters (1,1) struct
end

if contains(parameters.simulation.medium, {'layered'})
    
    % Set up a grid containing a layered medium
    % (1) Align planning image to transducer and rescale to requested grid resolution
    %       [multi-transducer] based on first transducer–focus pair in T1 space
    % (2) Convert segmentation to requested layers, crop grid to head +
    % transducer + PML
    %
    % Note: The new grid positions are not yet encoded in the parameters.
    % The original T1 positions are later transformed using the transf
    % matrix.

    % Preprocess layered head
    [medium_masks, segmentation, bone_mask, pseudoCT, ~, ~, planimg.t1_image_orig, ...
        planimg.t1_header, planimg.transf, planimg.inv_transf] = ...
        preproc_head(parameters);

    % update present layers (following medium mask creation)
    [parameters] = check_layers(parameters, segmentation);

    if isempty(medium_masks)
        return;
    end
    
    parameters.grid.dims  = size(medium_masks);
    planimg.origin_ras_mm = [];   % defined by T1 header for layered grids

else
    % In case simulations are not run in a skull of layered tissue, 
    % alternative grid dimensions are set up
    if strcmp(parameters.simulation.medium, 'phantom')
        % read in phantoms directly as medium masks 
        segmentation_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
        filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');
        segmented_img = squeeze(niftiread(filename_segmented));
        if size(segmented_img) == parameters.grid.default_dims
            parameters.grid.dims = parameters.grid.default_dims;
            disp('Check passed: phantom dimensions fit requested grid...');
        elseif length(size(segmented_img))==2 && all(size(segmented_img') == parameters.grid.default_dims)
            parameters.grid.dims = parameters.grid.default_dims;
            segmented_img = segmented_img';
            disp('Check passed: phantom dimensions fit requested grid after rotating phantom...');
        else
            parameters.grid.dims = size(segmented_img);
            disp('Setting grid according to phantom dimensions...')
            sprintf('%dD grid dimensions: [%d, %d, %d]', numel(parameters.grid.dims), parameters.grid.dims);
        end
        % create medium mask according to indices in parameters.layers (see preproc_smooth_and_crop.m)
        [medium_masks] = preproc_medium_mask(segmented_img, parameters);
        segmentation = segmented_img; clear segmented_img;
        mask = tissuemask_binary(parameters, medium_masks);
        bone_mask = mask.skull;
        pseudoCT  = [];
        % update present layers (following medium mask creation)
        [parameters] = check_layers(parameters, segmentation);
    else % e.g., water
        % set up default grid dimensions
        assert(isfield(parameters, 'grid') && isfield(parameters.grid, 'default_dims'), ...
            'parameters.grid.default_dims must be set for water simulations')
        parameters.grid.dims = squeeze(parameters.grid.default_dims);
        sprintf('Using default %dD grid dimensions: [%d, %d, %d]', numel(parameters.grid.dims), parameters.grid.dims);
        % set up empty medium masks and segmentations
        medium_masks = ones(parameters.grid.dims); % single water layer
        segmentation = zeros(parameters.grid.dims);
        bone_mask = zeros(parameters.grid.dims);
        pseudoCT  = [];
        % update present layer (only water)
        [parameters] = check_layers(parameters, segmentation);
    end

    % no planning image transform for non-layered grids
    planimg.t1_image_orig = [];
    planimg.t1_header     = [];
    planimg.transf        = [];
    planimg.inv_transf    = [];
    % optional RAS world-space anchor (set via parameters.grid.origin_ras_mm)
    if isfield(parameters, 'grid') && isfield(parameters.grid, 'origin_ras_mm')
        planimg.origin_ras_mm = parameters.grid.origin_ras_mm(:)';
    else
        planimg.origin_ras_mm = [];
    end

end