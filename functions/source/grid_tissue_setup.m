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
%                    transf, inv_transf); fields are empty for non-layered media.
%                    Phantoms add phantom_header, phantom_transposed (2D phantom
%                    transposed to default_dims) and phantom_slice (struct dim,
%                    index of the plane cut from a 3D phantom; [] otherwise)
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
    % parameters.transducer is updated with positions in the (RAS+) T1 grid
    % that planimg.transf starts from.
    [medium_masks, segmentation, bone_mask, pseudoCT, ~, ~, planimg.t1_image_orig, ...
        planimg.t1_header, planimg.transf, planimg.inv_transf, parameters.transducer] = ...
        preproc_head(parameters);

    % Planned positions for TARGETING_ANALYSIS must refer to the same (RAS+)
    % T1 grid as planimg.t1_header, which the T1w output NIfTIs carry.
    if isfield(parameters, 'state')
        for fn = {'trans_pos', 'focus_pos'}
            f = ['planned_' strrep(fn{1}, '_pos', '') '_t1_vox'];
            if isfield(parameters.state, f) && ~isempty(parameters.state.(f))
                parameters.state.(f) = cell2mat(arrayfun(@(t) double(t.(fn{1})(:)'), ...
                    parameters.transducer(:), 'UniformOutput', false));
            end
        end
    end

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
        nii_info      = niftiinfo(filename_segmented);
        segmented_img = squeeze(niftiread(filename_segmented));
        planimg.phantom_header = nii_info;
        % What is done to the phantom array below, so output NIfTIs can be
        % returned to the file's voxel order (PHANTOM_VOXEL_ORDER)
        planimg.phantom_transposed = false;
        planimg.phantom_slice      = [];

        % Read voxel size from NIfTI header; the grid has one resolution, so
        % the phantom voxels must be isotropic in its spatial dims
        spatial = find(nii_info.ImageSize(1:min(3, end)) > 1);
        pix_mm  = nii_info.PixelDimensions(spatial);
        if any(abs(pix_mm - pix_mm(1)) > 1e-6 * max(1, pix_mm(1)))
            error('grid_tissue_setup:anisotropicPhantom', ...
                ['Phantom NIfTI voxels are not isotropic (%s mm): %s\n' ...
                 'Resample the phantom to isotropic voxels.'], num2str(pix_mm), filename_segmented);
        end
        nii_res_mm = pix_mm(1);
        if abs(nii_res_mm - parameters.grid.resolution_mm) > 1e-6
            warning(['Phantom NIfTI voxel size (%.4f mm) differs from parameters.grid.resolution_mm (%.4f mm). ' ...
                     'Overriding parameters.grid.resolution_mm with NIfTI value.'], ...
                    nii_res_mm, parameters.grid.resolution_mm);
            parameters.grid.resolution_mm = nii_res_mm;
        end

        if isequal(size(segmented_img), parameters.grid.default_dims)
            parameters.grid.dims = parameters.grid.default_dims;
            fprintf('[grid_tissue_setup] %s\n', 'Check passed: phantom dimensions fit requested grid...');
        elseif length(size(segmented_img))==2 && isequal(size(segmented_img'), parameters.grid.default_dims)
            parameters.grid.dims = parameters.grid.default_dims;
            segmented_img = segmented_img';
            planimg.phantom_transposed = true;
            fprintf('[grid_tissue_setup] %s\n', 'Check passed: phantom dimensions fit requested grid after rotating phantom...');
        elseif isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric && ...
                ndims(segmented_img) == 3 && numel(parameters.grid.default_dims) == 2
            % 2D axisymmetric sim with a 3D phantom: extract the dim-2 plane
            % through the transducer
            [slice_y, parameters] = phantom_axis_slice(parameters, size(segmented_img));
            segmented_img = squeeze(segmented_img(:, slice_y, :));
            planimg.phantom_slice = struct('dim', 2, 'index', slice_y);
            if isequal(size(segmented_img), parameters.grid.default_dims)
                parameters.grid.dims = parameters.grid.default_dims;
            else
                parameters.grid.dims = size(segmented_img);
            end
            fprintf('[grid_tissue_setup] Axisymmetric 2D: extracted slice (y=%d) from 3D phantom.\n', slice_y);
            fprintf('[grid_tissue_setup] %dD grid dimensions: [%s]\n', numel(parameters.grid.dims), num2str(parameters.grid.dims));
        else
            parameters.grid.dims = size(segmented_img);
            fprintf('[grid_tissue_setup] %s\n', 'Setting grid according to phantom dimensions...')
            fprintf('[grid_tissue_setup] %dD grid dimensions: [%s]\n', numel(parameters.grid.dims), num2str(parameters.grid.dims));
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
        if ischar(parameters.grid.default_dims) || isstring(parameters.grid.default_dims)
            assert(strcmpi(parameters.grid.default_dims, 'auto'), 'grid_tissue_setup:defaultDims', ...
                'grid.default_dims must be numeric or ''auto''.');
            % Annular transducers are placed on the auto grid (GRID_TRANSDUCER_LOCATION
            % clears their positions in water); other types keep trans_pos/focus_pos,
            % which are indices of a grid of known size: keep the grid they were
            % set for
            tr1 = parameters.transducer(1);
            pos_set = any(cellfun(@(f) isfield(tr1, f) && ~isempty(tr1.(f)), ...
                {'trans_pos', 'focus_pos', 'trans_pos_mm', 'focus_pos_mm'}));
            if ~strcmp(tr1.type, 'annular') && pos_set
                legacy = [144, 144, 400];
                if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric; legacy = [144, 400]; end
                warning('grid_tissue_setup:autoDimsPositions', ['grid.default_dims ''auto'' with a %s transducer ' ...
                    'placed by trans_pos/focus_pos (grid indices): using [%s]. Set grid.default_dims explicitly.'], ...
                    tr1.type, num2str(legacy));
                parameters.grid.default_dims = legacy;
            end
        end
        if ischar(parameters.grid.default_dims) || isstring(parameters.grid.default_dims)
            % Sized to the transducer aperture and focal depth (WATER_GRID_DIMS)
            [parameters.grid.default_dims, gi] = water_grid_dims(parameters);
            fprintf('[grid_tissue_setup] default_dims auto: aperture %.1f mm + margins → %.1f mm lateral, focus %.1f mm (bowl) + beyond → %.1f mm axial\n', ...
                gi.aperture_mm, gi.lateral_mm, gi.focal_distance_bowl_mm, gi.axial_mm);
        end
        parameters.grid.dims = squeeze(parameters.grid.default_dims);
        fprintf('[grid_tissue_setup] Using default %dD grid dimensions: [%s]\n', numel(parameters.grid.dims), num2str(parameters.grid.dims));
        % set up empty medium masks and segmentations
        medium_masks = ones(parameters.grid.dims); % single water layer
        segmentation = zeros(parameters.grid.dims);
        bone_mask = zeros(parameters.grid.dims);
        pseudoCT  = [];
        % update present layer (only water)
        [parameters] = check_layers(parameters, segmentation);
    end

    % no planning image transform for non-layered grids
    planimg.t1_image_orig  = [];
    planimg.t1_header      = [];
    planimg.transf         = [];
    planimg.inv_transf     = [];
    if ~isfield(planimg, 'phantom_header')
        planimg.phantom_header = [];
    end
    % optional RAS world-space anchor (set via parameters.grid.origin_ras_mm)
    if isfield(parameters, 'grid') && isfield(parameters.grid, 'origin_ras_mm')
        planimg.origin_ras_mm = parameters.grid.origin_ras_mm(:)';
    else
        planimg.origin_ras_mm = [];
    end

end
end

function [k, parameters] = phantom_axis_slice(parameters, sz)
% Dim-2 index of the plane cut from a 3D phantom for an axisymmetric run: the
% transducer's when transducer(1) gives a 3D trans_pos (or trans_pos_mm,
% converted as in GRID_TRANSDUCER_LOCATION), which is then reduced to the
% plane; else the centre slice.
    tr  = parameters.transducer(1);
    res = parameters.grid.resolution_mm;
    pos = struct();
    for f = {'trans_pos', 'focus_pos'}
        mm = [f{1} '_mm'];
        if isfield(tr, mm) && numel(tr.(mm)) == 3
            pos.(f{1}) = grid_mm_to_voxel(tr.(mm), res);
        elseif isfield(tr, f{1}) && numel(tr.(f{1})) == 3
            pos.(f{1}) = tr.(f{1})(:)';
        end
    end
    if ~isfield(pos, 'trans_pos')
        k = round(sz(2) / 2);
        if mod(sz(2), 2) == 0
            warning('grid_tissue_setup:ambiguousPhantomCentre', ...
                ['The 3D phantom has an even number of slices (%d) in dim 2, so its centre lies ' ...
                 'between two slices; using slice %d. Give transducer.trans_pos in 3D to choose it.'], ...
                sz(2), k);
        end
        return
    end
    k = pos.trans_pos(2);
    if k < 1 || k > sz(2)
        error('grid_tissue_setup:phantomSliceOutside', ...
            'transducer.trans_pos dim-2 index %d lies outside the 3D phantom (1..%d).', k, sz(2));
    end
    if isfield(pos, 'focus_pos') && pos.focus_pos(2) ~= k
        error('grid_tissue_setup:phantomSliceAxis', ...
            ['Axisymmetric run on a 3D phantom: transducer and focus must lie in the same dim-2 ' ...
             'slice (trans_pos %d, focus_pos %d).'], k, pos.focus_pos(2));
    end
    % positions in the extracted [dim1 x dim3] plane
    for f = {'trans_pos', 'focus_pos', 'trans_pos_mm', 'focus_pos_mm'}
        if isfield(tr, f{1}) && numel(tr.(f{1})) == 3
            parameters.transducer(1).(f{1}) = tr.(f{1})([1 3]);
        end
    end
end
