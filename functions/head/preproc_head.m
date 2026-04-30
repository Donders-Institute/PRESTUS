function [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, ...
    focus_pos_final, t1_image, t1_header, final_transformation_matrix, inv_final_transformation_matrix] = preproc_head(parameters)

% PREPROC_HEAD  Preprocess structural MRI for k-Wave simulation
%
% Combines SimNIBS segmentation, focal-axis alignment, and grid cropping
% to prepare tissue masks, bone images, and transducer/focus positions for
% k-Wave simulation. Produces diagnostic figures and saves intermediate
% outputs.
%
% Use as:
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, ...
%    focus_pos_final, t1_image, t1_header, final_transformation_matrix, ...
%    inv_final_transformation_matrix] = preproc_head(parameters)
%
% Input:
%   parameters - (1,1) simulation configuration struct
%
% Output:
%   medium_masks                    - [Nx x Ny x Nz] tissue medium label array
%   segmentation_crop               - [Nx x Ny x Nz] cropped tissue segmentation
%   bone_mask_crop                  - [Nx x Ny x Nz] cropped binary skull mask (always present)
%   pseudoCT_crop                   - [Nx x Ny x Nz] cropped Hounsfield-unit skull image
%                                     ([] when parameters.pct.enabled ~= 1)
%   trans_pos_final                 - [1x3] transducer position in final grid
%   focus_pos_final                 - [1x3] focus position in final grid
%   t1_image                        - [Nx x Ny x Nz] T1-weighted image
%   t1_header                       - NIfTI header struct
%   final_transformation_matrix     - [4x4] combined alignment + crop transform
%   inv_final_transformation_matrix - [4x4] inverse of final_transformation_matrix
%
% See also: PREPROC_ALIGN_TO_FOCAL_AXIS, HEAD_SMOOTH_AND_CROP, PREPROC_CROP_GRID

    arguments
        parameters (1,1) struct
    end

    %% CHECK INPUTS

    disp('Checking inputs for head preprocessing ...');
    
    % Define path to T1 image (user-specified)
    filename_t1 = fullfile(parameters.path.anat, sprintf(parameters.path.t1_pattern, parameters.subject_id));

    % Define path to segmentation results
    segmentation_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
    filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');

    % Define path to T1 image (simnibs; aligned with segmentation space)
    filename_t1_simnibs = fullfile(segmentation_folder, 'T1.nii.gz');

    % Validate existence of files
    files_to_check = {filename_t1, filename_segmented, filename_t1_simnibs};
    check_availability(files_to_check)

    %% Load planning, mask, and segmentation images + header information
    % Default: use T1 planning image from the segmentation directory
    % Experimental: use user-specified T1 planning image when position is based on localite

    disp('Loading images...');

    if isfield(parameters.placement,'localite') && isfield(parameters.placement.localite,'enabled') && parameters.placement.localite.enabled
        t1_image = niftiread(filename_t1);
        t1_header = niftiinfo(filename_t1);
    else
        t1_image = niftiread(filename_t1_simnibs);
        t1_header = niftiinfo(filename_t1_simnibs);
    end

    if parameters.pct.enabled == 1
        % Load pseudoCT from the resolved pCT directory (parameters.io.dir_pct).
        % path_log_setup sets this to the SimNIBS m2m folder by default, or to
        % {path.pct}/sub-NNN/ when a dedicated pCT base directory is configured.
        pct_dir = parameters.io.dir_pct;
        filename_pseudoCT = fullfile(pct_dir, 'pseudoCT.nii.gz');
        pseudoCT_image = niftiread(filename_pseudoCT);
        pseudoCT_header = niftiinfo(filename_pseudoCT);

        % Load pseudoCT tissues mask
        filename_tissues_mask = fullfile(pct_dir, 'tissues_mask.nii.gz');
        tissues_mask_image = niftiread(filename_tissues_mask);
        tissues_mask_header = niftiinfo(filename_tissues_mask);
    else
        % Load traditional SimNIBS segmentation
        tissues_mask_image = niftiread(filename_segmented);
        tissues_mask_header = niftiinfo(filename_segmented);
    end

    %% Determine transducer position in T1 grid
    % based on configuration file [default] or localite neuronavigation data.
    % Note: localite coordinates refer to the planning T1 header, which may
    % differ from the SimNIBS segmentation T1.
    % [Multi-transducer] preprocessing is based on the first transducer only.

    if isfield(parameters.placement,'localite') && isfield(parameters.placement.localite,'enabled') && parameters.placement.localite.enabled

        lc = parameters.placement.localite;
        use_file = isfield(lc, 'file') && ~isempty(lc.file);
        use_neuronav = ~use_file && isfield(parameters, 'path') && ...
                       isfield(parameters.path, 'localite') && ~isempty(parameters.path.localite);

        if use_file
            % --- Simple mode: single XML file specified directly ---
            check_availability({lc.file})
            [trans_pos_grid, focus_pos_grid, ~, ~] = ...
                position_transducer_localite(lc.file, t1_header, parameters);

        elseif use_neuronav
            % --- Full neuronav mode: automatic file selection + series averaging ---
            % Required: parameters.path.localite (data folder)
            % Optional: placement.localite.session   (e.g., 'ses-01' or 1; default: 1)
            %           placement.localite.markertype (default: 'TriggerMarkers')
            %           placement.localite.position   (series index to use; default: 1)
            pn.data_postlocalite = parameters.path.localite;
            sub_id  = sprintf('sub-%03d', parameters.subject_id);
            ses_id  = 1;
            if isfield(lc, 'session') && ~isempty(lc.session), ses_id = lc.session; end
            markertype = 'TriggerMarkers';
            if isfield(lc, 'markertype') && ~isempty(lc.markertype), markertype = lc.markertype; end
            position = 1;
            if isfield(lc, 'position') && ~isempty(lc.position), position = lc.position; end

            localite = neuronav_select_localite(pn, sub_id, ses_id, markertype);
            if isempty(localite)
                error('PRESTUS:localite:noFile', ...
                    'neuronav_select_localite found no valid %s file for %s session %s in:\n  %s', ...
                    markertype, sub_id, num2str(ses_id), parameters.path.localite);
            end

            stats = neuronav_compute_series_statistics(localite, 1, [], markertype);
            if isempty(stats) || position > numel(stats)
                error('PRESTUS:localite:noSeries', ...
                    'Requested position index %d not found in Localite data (%d series detected).', ...
                    position, numel(stats));
            end

            coord_matrix = reshape(squeeze(stats{position}.matrix4d_mean), [4, 4])';
            [trans_pos_grid, focus_pos_grid] = ...
                localite_matrix_to_positions(coord_matrix, t1_header, parameters);

        else
            error('PRESTUS:localite:noSource', ...
                ['placement.localite.enabled = 1 but no Localite source is configured. ' ...
                 'Set placement.localite.file (single XML) or path.localite (folder for ' ...
                 'automatic file selection).']);
        end
    else
        if numel(parameters.transducer) > 1
            warn('Multiple transducers defined; alignment, cropping, and intermediate debug plots will be based only on the first transducer.');
        end
        trans_pos_grid = parameters.transducer(1).trans_pos;
        focus_pos_grid = parameters.transducer(1).focus_pos;
    end

    % ensure that dimensions of transducer in grid are sensible
    % transducer position should be [1 x 2/3]
    if size(trans_pos_grid, 1) > size(trans_pos_grid, 2)
        warn('Transducer may be incorrectly specified. Switching dimensions ...')
        trans_pos_grid = trans_pos_grid';
        focus_pos_grid = focus_pos_grid';
    end

    %% Rotate and scale images to match the stimulation trajectory

    disp('Rotating images to focal axis and rescaling to grid resolution ...')

    % If the segmentation process was not successful, it will stop preprocessing
    assert(exist(filename_segmented,'file') > 0, ...
        'Head segmentation is not completed (%s does not exist), see logs in the batch_logs folder and in %s folder',...
            filename_segmented, segmentation_folder)

    % Defines output file location and name.
    % Use io.preproc_affix when set (uncertainty pipeline stages 2–4 point
    % this at stage 1's cached files with affix ''), otherwise io.output_affix.
    if isfield(parameters.io, 'preproc_affix')
        preproc_affix = parameters.io.preproc_affix;
    else
        preproc_affix = parameters.io.output_affix;
    end
    filename_reoriented_scaled_data = fullfile(parameters.io.dir_cache, ...
        sprintf('sub-%03d_%s_rotated_scaled%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, preproc_affix));

    if confirm_overwriting(filename_reoriented_scaled_data, parameters)

        log_timer('start','preproc_rotscale', parameters.io.dir_output);

        %% [Planning image] (rotation matrix will be established for planning image)

        % Note: This primarily serves plotting purposes.
        if t1_header.ImageSize(3) > 1
            % rescale the image to the desired grid resolution
            scale_factor_t1 = t1_header.PixelDimensions(1)/parameters.grid.resolution_mm;
            [t1_img_rr, trans_pos_rescaled, focus_pos_rescaled, ...
            scale_rotate_recenter_matrix, rotation_matrix, ~, ~, t1_rr_img_montage] = ...
                preproc_align_to_focal_axis(...
                t1_image, ...
                t1_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                scale_factor_t1, ...
                parameters);

            % [DEBUG] visualize original and rotated planning image
            if parameters.simulation.debug == 1
                h = figure;
                imshow(t1_rr_img_montage)
                title('Original (left) and rotated (right) planning image');
                output_plot_filename = fullfile(parameters.io.dir_debug_preproc, ...
                    sprintf('sub-%03d_t1_rotated_scaled%s.png',...
                    parameters.subject_id, parameters.io.output_affix));
                saveas(h, output_plot_filename, 'png');
                close(h);
            end; clear t1_rr_img_montage;

        else
            warn('Currently no support for reorienting 2D planning image...');
        end
        
        %% [Tissue segmentation]

        scale_factor_seg = tissues_mask_header.PixelDimensions(1)/parameters.grid.resolution_mm;
        [segmented_img_rr, ~, ~, ...
            ~, ~, ~, ~, segm_img_montage] = ...
            preproc_align_to_focal_axis(...
            tissues_mask_image, ...
            tissues_mask_header, ...
            trans_pos_grid, ...
            focus_pos_grid, ...
            scale_factor_seg, ...
            parameters);
        
        % [DEBUG] visualize original and rotated segmentation image
        if parameters.simulation.debug == 1
            h = figure;
            imshow(segm_img_montage)
            title('Original (left) and rotated (right) tissue segmentation');
            output_plot_filename = fullfile(parameters.io.dir_debug_preproc, ...
                sprintf('sub-%03d_segmented_rotated_scaled%s.png',...
                parameters.subject_id, parameters.io.output_affix));
            saveas(h, output_plot_filename, 'png');
            close(h);
        end;  clear segm_img_montage;
        
        %% [bone mask and pseudoCT]

        % Binary skull mask is always derived from the segmentation labels.
        bone_img = ismember(tissues_mask_image, charm_seg_labels().bonemask);
        [bone_mask_img_rr, ~, ~, ~, ~, ~, ~, bone_mask_montage] = ...
            preproc_align_to_focal_axis(...
            bone_img, ...
            tissues_mask_header, ...
            trans_pos_grid, ...
            focus_pos_grid, ...
            tissues_mask_header.PixelDimensions(1)/parameters.grid.resolution_mm, ...
            parameters);
        clear bone_img;

        % pseudoCT (HU values) is only available in the pCT pipeline variant.
        if parameters.pct.enabled == 1
            [pseudoCT_img_rr, ~, ~, ~, ~, ~, ~, pseudoCT_montage] = ...
                preproc_align_to_focal_axis(...
                pseudoCT_image, ...
                pseudoCT_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                pseudoCT_header.PixelDimensions(1)/parameters.grid.resolution_mm, ...
                parameters);
        else
            pseudoCT_img_rr = [];
        end

        % [DEBUG] visualize original and rotated bone mask
        if parameters.simulation.debug == 1
            h = figure;
            imshow(bone_mask_montage)
            title('Original (left) and rotated (right) bone mask');
            output_plot_filename = fullfile(parameters.io.dir_debug_preproc, ...
                sprintf('sub-%03d_rotated_scaled_orig%s.png',...
                parameters.subject_id, parameters.io.output_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
            if parameters.pct.enabled == 1
                h = figure;
                imshow(pseudoCT_montage)
                title('Original (left) and rotated (right) pseudoCT');
                output_plot_filename = fullfile(parameters.io.dir_debug_preproc, ...
                    sprintf('sub-%03d_rotated_scaled_pseudoCT%s.png',...
                    parameters.subject_id, parameters.io.output_affix));
                saveas(h, output_plot_filename, 'png')
                close(h);
            end
        end; clear bone_mask_montage pseudoCT_montage;

        assert(isequal(size(trans_pos_rescaled(1:2)),size(focus_pos_rescaled(1:2))),...
            "After reorientation, the first two coordinates of the focus and the transducer should be the same")

        % Save rotated and rescaled data.
        % The transformation matrices are always saved — they are small and
        % needed for any downstream use of this checkpoint file.
        % The large image volumes (t1_img_rr, segmented_img_rr, bone_img_rr)
        % are only needed as inputs to the smoothing/cropping step (file 2).
        % Once file 2 exists they are never used again, so we omit them from
        % the checkpoint to save disk space.
        filename_cropped_smoothed_skull_data_check = fullfile(parameters.io.dir_cache, ...
            sprintf('sub-%03d_%s_cropped_smoothed%s.mat',...
            parameters.subject_id, parameters.simulation.medium, preproc_affix));
        if ~isfile(filename_cropped_smoothed_skull_data_check)
            % File 2 does not yet exist — include volumes so the next block
            % can load them if this session is interrupted and restarted.
            save(filename_reoriented_scaled_data, ...
                't1_img_rr', ...
                'segmented_img_rr', ...
                'bone_mask_img_rr', ...
                'pseudoCT_img_rr', ...
                'trans_pos_rescaled', ...
                'focus_pos_rescaled', ...
                'scale_rotate_recenter_matrix', ...
                'rotation_matrix');
        else
            % File 2 already exists — volumes are no longer needed; save
            % only the small transformation matrices.
            save(filename_reoriented_scaled_data, ...
                'trans_pos_rescaled', ...
                'focus_pos_rescaled', ...
                'scale_rotate_recenter_matrix', ...
                'rotation_matrix');
        end
        log_timer('stop','preproc_rotscale');
    else
        disp('Skipping reorientation; the file already exists. Loading it instead.');
        % Only load the large volumes if file 2 does not yet exist —
        % they are only needed as inputs to the smoothing/cropping step.
        filename_cropped_smoothed_skull_data_check = fullfile(parameters.io.dir_cache, ...
            sprintf('sub-%03d_%s_cropped_smoothed%s.mat',...
            parameters.subject_id, parameters.simulation.medium, preproc_affix));
        cropped_exists = isfile(filename_cropped_smoothed_skull_data_check);
        if ~cropped_exists
            load(filename_reoriented_scaled_data);   % includes volumes + matrices
            % Backward compat: old caches stored bone_img_rr as the dual-purpose variable.
            if exist('bone_img_rr', 'var') && ~exist('bone_mask_img_rr', 'var')
                if parameters.pct.enabled == 1
                    error('prestus:preproc:oldCacheFormat', ...
                        ['Cache %s uses the old bone_img_rr format (pseudoCT and bone mask ' ...
                         'were combined). Delete it to regenerate with separated variables.'], ...
                        filename_reoriented_scaled_data);
                else
                    bone_mask_img_rr = bone_img_rr;
                    pseudoCT_img_rr  = [];
                    clear bone_img_rr;
                end
            end
        else
            load(filename_reoriented_scaled_data, ...
                'trans_pos_rescaled', 'focus_pos_rescaled', ...
                'scale_rotate_recenter_matrix', 'rotation_matrix');
        end
    end

    % volumes_available tracks whether t1_img_rr / segmented_img_rr /
    % bone_img_rr are in the workspace (they are omitted on load when
    % cropped_smoothed.mat already exists).
    volumes_available = exist('t1_img_rr', 'var');

    %% [DEBUG] Plot the skin & skull from the segmented image and an overlay for comparison

    if parameters.simulation.debug == 1 && volumes_available
        % Create a T1 slice for comparison to SimNIBS segmented data
        t1_slice = repmat(mat2gray(squeeze(t1_img_rr(:,trans_pos_rescaled(2),:))), [1 1 3]);
        % Create slices of segmented SimNIBS data
        % Segmentation labels can be different for the pseudoCT mask. See create_pseudoCT.sh.
        skullidx = getidx(parameters.layers,'skull'); % Includes all "skull_" segmentations
        skinidx = getidx(parameters.layers,'skin');
        % Unsmoothed skull & skin masks
        skull_mask_unsmoothed = ismember(segmented_img_rr,skullidx);
        skin_mask_unsmoothed = ismember(segmented_img_rr,skinidx);
        skin_slice = squeeze(skin_mask_unsmoothed(:,trans_pos_rescaled(2),:));
        skull_slice = squeeze(skull_mask_unsmoothed(:,trans_pos_rescaled(2),:));
        skin_skull_img = cat(3, zeros(size(skull_slice)), skin_slice, skull_slice); % Skull is blue; skin is green.

        h = figure;
        montage(cat(4,t1_slice*255 ,skin_skull_img*255 ,...
            imfuse(mat2gray(t1_slice), skin_skull_img,'blend')) ,'size',[1 NaN]);
        title('T1 and SimNIBS skin (green) and skull (blue) masks');
        output_plot_filename = fullfile(parameters.io.dir_debug_preproc,...
            sprintf('sub-%03d_t1_skin_skull%s.png',parameters.subject_id, parameters.io.output_affix));
        saveas(h ,output_plot_filename ,'png')
        close(h);

        clear t1_slice skullidx skinidx skull_mask_unsmoothed skin_mask_unsmoothed skin_slice skull_slice skin_skull_img;
    end
    
    %% SEG-> MEDIUM, SMOOTH & CROP LAYERED
    % see preproc_smooth_and_crop
    
    disp('Creating layered medium by smoothing and cropping the head segmentation...')

    % Defines output file location and name
    filename_cropped_smoothed_skull_data = fullfile(parameters.io.dir_cache, ...
        sprintf('sub-%03d_%s_cropped_smoothed%s.mat',...
        parameters.subject_id, parameters.simulation.medium, preproc_affix));
    
   if confirm_overwriting(filename_cropped_smoothed_skull_data, parameters)
        log_timer('start','preproc_skullsmooth', parameters.io.dir_output);
        % postprocess skull segmentation
        [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, focus_pos_final, crop_translation_matrix] = ...
            head_smooth_and_crop(...
            parameters, segmented_img_rr, bone_mask_img_rr, trans_pos_rescaled, focus_pos_rescaled, pseudoCT_img_rr);

        % Combine transformations from focal axis alignment with crop
        final_transformation_matrix = scale_rotate_recenter_matrix*crop_translation_matrix';
        inv_final_transformation_matrix = maketform('affine', inv(final_transformation_matrix')');

        % Diagnostic NIfTIs: medium_masks, segmentation, and skull_mask
        % back-projected to T1 space. Written to debug_dir/preproc/ only
        % when simulation.debug == 1 — not intended as pipeline outputs.
        if parameters.simulation.debug == 1
            % save medium mask
            orig_hdr = t1_header; % header is based on original T1w (always present)
            orig_hdr.Datatype = 'single';
            segmented_file = fullfile(parameters.io.dir_debug_preproc,...
                sprintf('sub-%03d_medium_masks_final', parameters.subject_id));
            plotdata = single(tformarray(uint8(medium_masks), inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
            if ~isfile(segmented_file)
                try
                    niftiwrite(plotdata, segmented_file, orig_hdr, 'Compressed',true);
                catch ME
                    warn(ME)
                end
            end
            clear segmented_file plotdata orig_hdr;
    
            % save segmentation skull mask/pseudoCT
            orig_hdr = t1_header; % header is based on original T1w (always present)
            orig_hdr.Datatype = 'single';
            segmentation_file = fullfile(parameters.io.dir_debug_preproc,...
                sprintf('sub-%03d_segmentation_final', parameters.subject_id));
            plotdata = single(tformarray(uint8(segmentation_crop), inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
            if ~isfile(segmentation_file)
                try
                    niftiwrite(plotdata, segmentation_file, orig_hdr, 'Compressed',true);
                catch ME
                    warn(ME)
                end
            end
            clear segmentation_file plotdata orig_hdr;

            % save binary skull mask
            orig_hdr = t1_header; % header is based on original T1w (always present)
            orig_hdr.Datatype = 'double';
            skull_mask_file = fullfile(parameters.io.dir_debug_preproc,...
                sprintf('sub-%03d_skull_mask_final', parameters.subject_id));
            plotdata = double(tformarray(bone_mask_crop, inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0));
            if ~isfile(skull_mask_file)
                try
                    niftiwrite(plotdata, skull_mask_file, orig_hdr, 'Compressed',true);
                catch ME
                    warn(ME)
                end
            end
            clear skull_mask_file plotdata;

            % save pseudoCT (if available)
            if ~isempty(pseudoCT_crop)
                orig_hdr.Datatype = 'single';
                pseudoCT_file = fullfile(parameters.io.dir_debug_preproc,...
                    sprintf('sub-%03d_pseudoCT_final', parameters.subject_id));
                plotdata = single(tformarray(pseudoCT_crop, inv_final_transformation_matrix, ...
                    makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0));
                if ~isfile(pseudoCT_file)
                    try
                        niftiwrite(plotdata, pseudoCT_file, orig_hdr, 'Compressed',true);
                    catch ME
                        warn(ME)
                    end
                end
                clear pseudoCT_file plotdata;
            end
            clear orig_hdr;
        end

        % save cropped and smoothed skull data
        save(filename_cropped_smoothed_skull_data, ...
            'medium_masks', ...
            'segmentation_crop', ...
            'bone_mask_crop', ...
            'pseudoCT_crop', ...
            'trans_pos_final', ...
            'focus_pos_final', ...
            'crop_translation_matrix',...
            'final_transformation_matrix',...
            'inv_final_transformation_matrix')
        log_timer('stop','preproc_skullsmooth');
    else
        disp('Skipping head smoothing and cropping, the file already exists, loading it instead.')
        load(filename_cropped_smoothed_skull_data);
        % Backward compat: old caches used bone_crop for the dual-purpose variable.
        if exist('bone_crop', 'var') && ~exist('bone_mask_crop', 'var')
            if parameters.pct.enabled == 1
                error('prestus:preproc:oldCacheFormat', ...
                    ['Cache %s uses the old bone_crop format (pseudoCT and bone mask were ' ...
                     'combined). Delete it to regenerate with separated variables.'], ...
                    filename_cropped_smoothed_skull_data);
            else
                bone_mask_crop = bone_crop;
                pseudoCT_crop  = [];
                clear bone_crop;
            end
        end
    end
    
    %% Sanity-check transformations

    % Check that transformations are correct by inverting locations and comparing them to the original 
    % If the transformation cannot be correctly inverted, this will be displayed
    % If the inverse transformation is off by > 1 voxel, the script exits

    backtransf_coordinates = round(tformfwd([trans_pos_final; focus_pos_final], inv_final_transformation_matrix));
    diff_voxels = abs(backtransf_coordinates - [trans_pos_grid; focus_pos_grid]);
    
    if any(diff_voxels(:) > 1)
        warn('Backtransformed positions exceed 1-voxel tolerance.');
        disp('Original coordinates:'); disp([trans_pos_grid; focus_pos_grid]);
        disp('Final local:'); disp([trans_pos_final; focus_pos_final]);
        disp('Backtransformed:'); disp(backtransf_coordinates);
        disp('Voxel differences:'); disp(diff_voxels);
        error('Coordinate mismatch >1 voxel/dim—check cropping/padding.');
    elseif any(diff_voxels(:) > 0)
        warn('Minor backtransform offset (≤1 voxel): %s', mat2str(diff_voxels));
    end

    %% Plot placement of up to 2 Transducers

    % Update the parameters (for plottign below only)  
    parameters.grid.dims = size(medium_masks);

    max_plots = min(2, numel(parameters.transducer));
    if numel(parameters.transducer) > max_plots
        warn('More than two transducers defined; only the first 2 will be shown in debug plots');
    end
    for ti = 1:max_plots            
        if ti == 1
            % canonical pair already used to define the transform
            tpos_t1 = trans_pos_grid;
            fpos_t1 = focus_pos_grid;
            tpos_sim = trans_pos_final;
            fpos_sim = focus_pos_final;
        else % ti == 2
            tr = parameters.transducer(ti);
            
            % T1-grid positions
            tpos_t1 = tr.trans_pos(:).';
            fpos_t1 = tr.focus_pos(:).';
    
            % map to simulation grid using the same global transform
            pts_sim = round(tformfwd([tpos_t1; fpos_t1], maketform('affine', final_transformation_matrix))); % 2×3
            tpos_sim = pts_sim(1,:);
            fpos_sim = pts_sim(2,:);
        end

        if parameters.simulation.debug == 1
            % [DEBUG] Plot brain segmentation
            [seg_with_trans_img, ~] = plot_t1_with_transducer(...
                medium_masks, parameters.grid.resolution_mm, tpos_sim, fpos_sim, parameters);

            h = figure;
            imshow(seg_with_trans_img);
            title('Segmentation with transducer');
            output_plot_filename = fullfile(parameters.io.dir_debug_preproc, ...
                sprintf('sub-%03d_%s_segmented_brain_final_T%02d%s.png',...
                parameters.subject_id, parameters.simulation.medium, ti, parameters.io.output_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
        end

        % Plot positioning of transducer on segmentation
        if numel(parameters.transducer) > 1
            tnn_suffix = sprintf('_T%02d', ti);
        else
            tnn_suffix = '';
        end
        output_plot_filename = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_positioning%s%s.png', ...
            parameters.subject_id, parameters.io.output_affix, tnn_suffix));

        % Subplot 1: Original segmentation with initial transducer and focus positions
        % Subplot 2: Original segmentation with slice cap applied
        % Subplot 3: Final segmentation with updated transducer and focus positions
        show_positioning_plots(...
            tissues_mask_image, ...
            tissues_mask_header.PixelDimensions(1), ...
            tpos_t1, ...
            fpos_t1, ...
            medium_masks, ...
            tpos_sim, ...
            fpos_sim, ...
            parameters, ...
            output_plot_filename);
    end

end