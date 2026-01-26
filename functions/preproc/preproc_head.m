function [medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, focus_pos_final, ...
    t1_image, t1_header, final_transformation_matrix, inv_final_transformation_matrix] = preproc_head(parameters)

% PREPROC_HEAD Preprocessing of structural data for simulations.
%
% This function combines MATLAB and SimNIBS to segment, realign, and crop structural data 
% for use in simulations. It prepares the data for k-wave simulations by segmenting tissues 
% (e.g., skin, bone, neural tissue), aligning the transducer's axis, and creating figures 
% to visualize simulation results.
%
% Input:
%   parameters  - Struct containing simulation parameters (e.g., paths, transducer settings).
%   .subject_id  - Integer specifying the subject ID.
%
% Output:
%   medium_masks                - Processed masks for different tissue types.
%   segmented_image_cropped     - Cropped segmented image after preprocessing.
%   skull_edge                  - Edge information for the skull segmentation.
%   trans_pos_final             - Final transducer position after preprocessing.
%   focus_pos_final             - Final focus position after preprocessing.
%   t1_image                    - Original T1-weighted image.
%   t1_header                   - Header information for the T1-weighted image.
%   final_transformation_matrix - Transformation matrix combining alignment and cropping steps.
%   inv_final_transformation_matrix - Inverse transformation matrix.

    %% CHECK INPUTS

    disp('Checking inputs for head preprocessing ...');
    
    % Define path to T1 image (user-specified)
    filename_t1 = fullfile(parameters.data_path, sprintf(parameters.t1_path_template, parameters.subject_id));

    % Define path to segmentation results
    segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', parameters.subject_id));
    if strcmp(parameters.segmentation_software, 'charm')
        filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');
    else 
        filename_segmented = fullfile(segmentation_folder, sprintf('sub-%03d_final_contr.nii.gz', parameters.subject_id));
    end

    % Define path to T1 image (simnibs; aligned with segmentation space)
    filename_t1_simnibs = fullfile(segmentation_folder, 'T1.nii.gz');

    % Validate existence of files
    files_to_check = {filename_t1, filename_segmented, filename_t1_simnibs};
    check_availability(files_to_check)

    %% Load planning, mask, and segmentation images + header information
    % Default: use T1 planning image from the segmentation directory
    % Experimental: use user-specified T1 planning image when position is based on localite

    disp('Loading images...');

    if isfield(parameters,'transducer_from_localite') && parameters.transducer_from_localite
        t1_image = niftiread(filename_t1);
        t1_header = niftiinfo(filename_t1);
    else
        t1_image = niftiread(filename_t1_simnibs);
        t1_header = niftiinfo(filename_t1_simnibs);
    end

    if parameters.usepseudoCT == 1
        % Load pseudoCT
        filename_pseudoCT = fullfile(segmentation_folder,'pseudoCT.nii.gz');
        pseudoCT_image = niftiread(filename_pseudoCT);
        pseudoCT_header = niftiinfo(filename_pseudoCT);
        
        % Load pseudoCT tissues mask
        filename_tissues_mask = fullfile(segmentation_folder,'tissues_mask.nii.gz');
        tissues_mask_image = niftiread(filename_tissues_mask);
        tissues_mask_header = niftiinfo(filename_tissues_mask);
    else
        % Load traditional SimNIBS segmentation
        tissues_mask_image = niftiread(filename_segmented);
        tissues_mask_header = niftiinfo(filename_segmented);
    end

    %% Determine transducer position in T1 grid
    % based on configuration file [default] or localite [experimental, undocumented]
    % Note: localite coordinates may refer to different header than e.g., simnibs segmentation
    % [Multi-transducer] the preprocessing will be based on the first transducer

    if isfield(parameters,'transducer_from_localite') && parameters.transducer_from_localite
        % Validate existence of localite file
        check_availability({localite_file})
        % Determine transducer position [experimental]
        [trans_pos_grid, focus_pos_grid, ~, ~] = ...
            position_transducer_localite(localite_file, t1_header, parameters);
    else
        if numel(parameters.transducer) > 1
            warning('Multiple transducers defined; alignment, cropping, and intermediate debug plots will be based only on the first transducer.');
        end
        trans_pos_grid = parameters.transducer(1).trans_pos;
        focus_pos_grid = parameters.transducer(1).focus_pos;
    end

    % ensure that dimensions of transducer in grid are sensible
    % transducer position should be [1 x 2/3]
    if size(trans_pos_grid, 1) > size(trans_pos_grid, 2)
        warning('Transducer may be incorrectly specified. Switching dimensions ...')
        trans_pos_grid = trans_pos_grid';
        focus_pos_grid = focus_pos_grid';
    end

    %% [Optional] Create unprocessed T1 slice oriented along the transducer's axis

    if parameters.debug==1
        h = figure;
        imshowpair(plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), ...
            trans_pos_grid, focus_pos_grid, parameters), ...
            plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), ...
            trans_pos_grid, focus_pos_grid, parameters, 'slice_dim', 1), 'montage');
        title('T1 with transducer');
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_t1_with_transducer_og_%s.png', ...
            parameters.subject_id, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);
    end

    %% Rotate and scale images to match the stimulation trajectory

    disp('Rotating images to focal axis and rescaling to grid resolution ...')

    % If the headreco process was not successful, it will stop preprocessing
    assert(exist(filename_segmented,'file') > 0, ...
        'Head segmentation is not completed (%s does not exist), see logs in the batch_logs folder and in %s folder',...
            filename_segmented, segmentation_folder)

    % Defines output file location and name
    filename_reoriented_scaled_data = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_after_rotating_and_scaling%s.mat', ...
        parameters.subject_id, parameters.results_filename_affix));

    if confirm_overwriting(filename_reoriented_scaled_data, parameters)

        log_timer('start','preproc_rotscale', parameters.output_dir);

        %% [Planning image] (rotation matrix will be established for planning image)

        % Note: This primarily serves plotting purposes.
        if t1_header.ImageSize(3) > 1
            % rescale the image to the desired grid resolution
            scale_factor_t1 = t1_header.PixelDimensions(1)/parameters.grid_step_mm;
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
            if parameters.debug == 1
                h = figure;
                imshow(t1_rr_img_montage)
                title('Original (left) and rotated (right) planning image');
                output_plot_filename = fullfile(parameters.debug_dir, ...
                    sprintf('sub-%03d_t1_after_rotating_and_scaling%s.png', ...
                    parameters.subject_id, parameters.results_filename_affix));
                saveas(h, output_plot_filename, 'png');
                close(h);
            end; clear t1_rr_img_montage;

        else
            warning('Currently no support for reorienting 2D planning image...');
        end
        
        %% [Tissue segmentation]

        scale_factor_seg = tissues_mask_header.PixelDimensions(1)/parameters.grid_step_mm;
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
        if parameters.debug == 1
            h = figure;
            imshow(segm_img_montage)
            title('Original (left) and rotated (right) tissue segmentation');
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_segmented_after_rotating_and_scaling%s.png', ...
                parameters.subject_id, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png');
            close(h);
        end;  clear segm_img_montage;
        
        %% [bone mask/pCT]

        if parameters.usepseudoCT == 1
            [bone_img_rr, ~, ~, ~, ~, ~, ~, bone_img_montage] = ...
                preproc_align_to_focal_axis(...
                pseudoCT_image, ...
                pseudoCT_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                pseudoCT_header.PixelDimensions(1)/parameters.grid_step_mm, ...
                parameters);
        else
            if strcmp(parameters.segmentation_software, 'charm') 
                % create filled bone mask as charm doesn't make it itself
                if isfield(parameters, 'seg_labels') && any(strcmp(fieldnames(parameters.seg_labels), 'bonemask'))
                    bone_img = ismember(tissues_mask_image,getidx(parameters.seg_labels,'bonemask'));
                else
                    bone_img = tissues_mask_image>0&(tissues_mask_image<=4|tissues_mask_image>=7);
                    warning("Using hardcoded labels for bonemask...");
                end
            else % load bone mask created by simnibs
                filename_bone_headreco = fullfile(segmentation_folder, 'bone.nii.gz');
                bone_img = niftiread(filename_bone_headreco);
            end
            [bone_img_rr, ~, ~, ~, ~, ~, ~, bone_img_montage] = ...
                preproc_align_to_focal_axis(...
                bone_img, ...
                tissues_mask_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                tissues_mask_header.PixelDimensions(1)/parameters.grid_step_mm, ...
                parameters);
        end

        % [DEBUG] visualize original and rotated bone mask image
        if parameters.debug == 1
            h = figure;
            imshow(bone_img_montage)
            title('Original (left) and rotated (right) original bone mask');
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_after_rotating_and_scaling_orig%s.png', ...
                parameters.subject_id, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h); 
        end; clear bone_img_montage;

        assert(isequal(size(trans_pos_rescaled(1:2)),size(focus_pos_rescaled(1:2))),...
            "After reorientation, the first two coordinates of the focus and the transducer should be the same")

        % Save rotated and rescaled data
        save(filename_reoriented_scaled_data, ...
            't1_img_rr', ...
            'segmented_img_rr', ...
            'bone_img_rr', ...
            'trans_pos_rescaled', ...
            'focus_pos_rescaled', ...
            'scale_rotate_recenter_matrix', ...
            'rotation_matrix');
        log_timer('stop','preproc_rotscale');
    else 
        disp('Skipping reorientation; the file already exists. Loading it instead.');
        load(filename_reoriented_scaled_data);
    end
    
    %% [DEBUG] Plot the skin & skull from the segmented image and an overlay for comparison

    if parameters.debug == 1
        % Create a T1 slice for comparison to SimNIBS segmented data
        t1_slice = repmat(mat2gray(squeeze(t1_img_rr(:,trans_pos_rescaled(2),:))), [1 1 3]);
        % Create slices of segmented SimNIBS data
        % Segmentation labels can be different for the pseudoCT mask. See create_pseudoCT.sh.
        skullidx = getidx(parameters.layer_labels,'skull'); % Includes all "skull_" segmentations
        skinidx = getidx(parameters.layer_labels,'skin');
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
        output_plot_filename = fullfile(parameters.debug_dir,...
            sprintf('sub-%03d_t1_skin_skull%s.png',parameters.subject_id, parameters.results_filename_affix));
        saveas(h ,output_plot_filename ,'png')
        close(h);

        clear t1_slice skullidx skinidx skull_mask_unsmoothed skin_mask_unsmoothed skin_slice skull_slice skin_skull_img;
    end
    
    %% SEG-> MEDIUM, SMOOTH & CROP LAYERED
    % see preproc_smooth_and_crop
    
    disp('Creating layered medium by smoothing and cropping the head segmentation...')

    % Defines output file location and name
    filename_cropped_smoothed_skull_data = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_after_cropping_and_smoothing%s.mat', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
   if confirm_overwriting(filename_cropped_smoothed_skull_data, parameters)
        log_timer('start','preproc_skullsmooth', parameters.output_dir);
        % postprocess skull segmentation
        [medium_masks, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, crop_translation_matrix] = ...
            head_smooth_and_crop(...
            parameters, segmented_img_rr, bone_img_rr, trans_pos_rescaled, focus_pos_rescaled);

        % Combine transformations from focal axis alignment with crop
        final_transformation_matrix = scale_rotate_recenter_matrix*crop_translation_matrix';
        inv_final_transformation_matrix = maketform('affine', inv(final_transformation_matrix')');

        % [DEBUG] save transformed medium mask and skull to debug dir
        if parameters.debug == 1
            % save medium mask
            orig_hdr = t1_header; % header is based on original T1w (always present)
            orig_hdr.Datatype = 'single';
            segmented_file = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_medium_masks_final', parameters.subject_id));
            plotdata = single(tformarray(uint8(medium_masks), inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
            if ~isfile(segmented_file)
                niftiwrite(plotdata, segmented_file, orig_hdr, 'Compressed',true);
            end
            clear segmented_file plotdata orig_hdr;
    
            % save skull mask/pseudoCT
            orig_hdr = t1_header; % header is based on original T1w (always present)
            orig_hdr.Datatype = 'single';
            skull_mask_file = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_skull_final', parameters.subject_id));
            plotdata = single(tformarray(uint8(segmented_image_cropped), inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
            if ~isfile(skull_mask_file)
                niftiwrite(plotdata, skull_mask_file, orig_hdr, 'Compressed',true);
            end
            clear skull_mask_file plotdata orig_hdr;
        end

        % save cropped and smoothed skull data
        save(filename_cropped_smoothed_skull_data, ...
            'medium_masks', ...
            'skull_edge', ...
            'segmented_image_cropped', ...
            'trans_pos_final', ...
            'focus_pos_final', ...
            'crop_translation_matrix',...
            'final_transformation_matrix',...
            'inv_final_transformation_matrix')
        log_timer('stop','preproc_skullsmooth');
    else 
        disp('Skipping head smoothing and cropping, the file already exists, loading it instead.')
        load(filename_cropped_smoothed_skull_data);
    end    
    
    %% Sanity-check transformations

    % Check that transformations are correct by inverting locations and comparing them to the original 
    % If the transformation cannot be correctly inverted, this will be displayed
    % If the inverse transformation is off by > 1 voxel, the script exits

    backtransf_coordinates = round(tformfwd([trans_pos_final; focus_pos_final], inv_final_transformation_matrix));
    diff_voxels = abs(backtransf_coordinates - [trans_pos_grid; focus_pos_grid]);
    
    if any(diff_voxels(:) > 1)
        warning('Backtransformed positions exceed 1-voxel tolerance.');
        disp('Original coordinates:'); disp([trans_pos_grid; focus_pos_grid]);
        disp('Final local:'); disp([trans_pos_final; focus_pos_final]);
        disp('Backtransformed:'); disp(backtransf_coordinates);
        disp('Voxel differences:'); disp(diff_voxels);
        error('Coordinate mismatch >1 voxel/dim—check cropping/padding.');
    elseif any(diff_voxels(:) > 0)
        warning('Minor backtransform offset (≤1 voxel): %s', mat2str(diff_voxels));
    end

    %% Plot placement of up to 2 Transducers

    % Update the parameters (for plottign below only)  
    parameters.grid_dims = size(medium_masks);

    max_plots = min(2, numel(parameters.transducer));
    if numel(parameters.transducer) > max_plots
        warning('More than two transducers defined; only the first 2 will be shown in debug plots');
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

        if parameters.debug == 1
            % [DEBUG] Plot brain segmentation
            [seg_with_trans_img, ~] = plot_t1_with_transducer(...
                medium_masks, parameters.grid_step_mm, tpos_sim, fpos_sim, parameters);
            
            h = figure;
            imshow(seg_with_trans_img);
            title('Segmentation with transducer');
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_segmented_brain_final_T%02d%s.png', ...
                parameters.subject_id, parameters.simulation_medium, ti, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
        end
        
        % Plot positioning of transducer on segmentation
        output_plot_filename = fullfile(parameters.output_dir, ...
            sprintf('sub-%03d_positioning_T%02d%s.png', ...
            parameters.subject_id, ti, parameters.results_filename_affix));

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