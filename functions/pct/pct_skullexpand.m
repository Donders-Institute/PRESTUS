function pct_skullexpand(seg_path, path_pct)
% PCT_SKULLEXPAND  Load SimNIBS NIfTIs, run skull rubber wrap, and save output
%
% Loads tissue labels and skull mask from seg_path, constructs the
% parameters struct for skull_rubber_wrap, runs the balloon inflation, and
% saves the updated skull mask back to pct_path.
%
% Use as:
%   pct_skullexpand(seg_path, path_pct)
%
% Input:
%   seg_path - path containing final_tissues.nii.gz (SimNIBS tissue labels)
%   path_pct - pCT processing folder containing skull_mask.nii.gz
%
% See also: SKULL_RUBBER_WRAP, PCT_SKULLMAPPING

    arguments
        seg_path (1,:) char
        path_pct (1,:) char
    end

    %% --- Paths ---
    final_tissues_path      = fullfile(seg_path, 'final_tissues.nii.gz');
    skull_path = fullfile(path_pct, 'skull_mask.nii.gz');

    assert(exist(final_tissues_path,'file')==2,      "Missing: %s", final_tissues_path);
    assert(exist(skull_path,'file')==2,"Missing: %s", skull_path);

    %% --- Load NIfTI volumes + headers ---
    info_tissues = niftiinfo(final_tissues_path);                 % header/meta [web:149]
    final_tissues = niftiread(info_tissues);                      % data [web:148]

    info_skull = niftiinfo(skull_path);             % header/meta [web:149]
    final_tissues_skull = niftiread(info_skull);                  % data [web:148]

    % Sanity: same grid
    assert(isequal(size(final_tissues), size(final_tissues_skull)), ...
        'final_tissues and final_tissues_skull have different sizes.');

    %% --- Build masks/inputs for skull_rubber_wrap ---
    % BW: binary skull mask
    BW = final_tissues_skull ~= 0;   % treat any nonzero as skull

    % segmented_img: SimNIBS tissue labeling volume
    segmented_img = final_tissues;

    % medium_masks: your pipeline assumes identical masks/seg image
    medium_masks = segmented_img;

    %% --- Parameters ---
    parameters = struct();
    parameters.simulation.debug = 1;

    parameters.path.seg   = seg_path;
    parameters.debug_path = path_pct;

    % If you use this elsewhere in your pipeline, keep it consistent:
    parameters.path.t1_pattern = 'T1.nii.gz';

    parameters.headmodel.skull_wrap_radius = 10;
    parameters.headmodel.skull_wrap_visualize = 0;

    parameters.io.debug_dir_preproc = path_pct;       % where skull_rubber_wrap_visualize writes images
    parameters.io.output_affix = '';

    % Grid step (mm) from header voxel size if available
    % niftiinfo.PixelDimensions is [dx dy dz] in mm for most NIfTIs.
    if isfield(info_tissues,'PixelDimensions') && numel(info_tissues.PixelDimensions) >= 3
        parameters.grid.resolution_mm = double(info_tissues.PixelDimensions(1)); % First dimension, assuming isometric
    else
        parameters.grid.resolution_mm = 1; % fallback
    end

    % SimNIBS tissue label conventions (charm)
    parameters.layers.brain = [1, 2];     % GM/WM
    parameters.layers.skin  = [5];        % skin

    %% --- Run skull wrap (must exist on path) ---
    skull_rubber_wrap(parameters, BW, medium_masks, segmented_img);

    %% Replace skull mask
    disp('Updating skull mask...')
    cd(path_pct);
    movefile('balloon_mask_final.nii.gz', skull_path);

end
