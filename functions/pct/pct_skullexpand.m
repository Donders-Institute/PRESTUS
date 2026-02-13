function pct_skullexpand(seg_path, path_pct)
%PCT_SKULLEXPAND Load SimNIBS NIfTIs + headers, run skull rubber wrap, save output.
%
% Expected inputs under seg_path:
%   - final_tissues.nii.gz            (SimNIBS tissue labels)
%   - final_tissues_skull.nii.gz      (binary skull mask; or skull label image)
%
% Output written under pct_path (and/or parameters.debug_dir).

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
    parameters.debug = 1;

    parameters.seg_path   = seg_path;
    parameters.debug_path = path_pct;

    % If you use this elsewhere in your pipeline, keep it consistent:
    parameters.t1_path_template = 'T1.nii.gz';

    parameters.wrapradius = 10;
    parameters.skullwrap_visualize = 0;

    parameters.debug_dir = path_pct;                 % where skull_rubber_wrap_visualize writes images
    parameters.results_filename_affix = '';

    % Grid step (mm) from header voxel size if available
    % niftiinfo.PixelDimensions is [dx dy dz] in mm for most NIfTIs.
    if isfield(info_tissues,'PixelDimensions') && numel(info_tissues.PixelDimensions) >= 3
        parameters.grid_step_mm = double(info_tissues.PixelDimensions(1)); % First dimension, assuming isometric
    else
        parameters.grid_step_mm = 1; % fallback
    end

    % SimNIBS tissue label conventions: set these to what YOUR data uses
    parameters.layers.brain = [1, 2];     % GM/WM (example)
    parameters.layers.skin  = [5];        % skin (example)
    parameters.seg_labels.csf = [3];      % CSF (example)

    %% --- Run skull wrap (must exist on path) ---
    skull_rubber_wrap(parameters, BW, medium_masks, segmented_img);

    %% Replace skull mask
    disp('Updating skull mask...')
    cd(path_pct);
    movefile('balloon_mask_final.nii.gz', skull_path);

end
