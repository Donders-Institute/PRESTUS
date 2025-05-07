function [medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, focus_pos_final, ...
    t1_image, t1_header, final_transformation_matrix, inv_final_transformation_matrix] = preprocess_brain(parameters, subject_id)

% PREPROCESS_BRAIN Preprocessing of structural data for simulations.
%
% This function combines MATLAB and SimNIBS to segment, realign, and crop structural data 
% for use in simulations. It prepares the data for k-wave simulations by segmenting tissues 
% (e.g., skin, bone, neural tissue), aligning the transducer's axis, and creating figures 
% to visualize simulation results.
%
% Input:
%   parameters  - Struct containing simulation parameters (e.g., paths, transducer settings).
%   subject_id  - Integer specifying the subject ID.
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

    %% CHECK INPUTS AND TRANSLATE PATTERNS
    disp('Checking inputs...');
    
    % Define paths to T1 and T2 images
    filename_t1 = fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id));
    filename_t2 = fullfile(parameters.data_path, sprintf(parameters.t2_path_template, subject_id));
    
    % Check if Localite transducer location data is enabled
    files_to_check = ["filename_t1", "filename_t2"];
    if isfield(parameters,'transducer_from_localite') && parameters.transducer_from_localite
        localite_file = fullfile(parameters.data_path, sprintf(parameters.localite_instr_file_template, subject_id));
        files_to_check = [files_to_check, "localite_file"];
    end

    % Validate existence of required files
    for filename_var = files_to_check
        eval(sprintf('filename = %s;',  filename_var ));
        if contains(filename, '*') % Handle wildcards in file paths
            matching_files = dir(filename);
            if length(matching_files) > 1
                error('More than 1 file matches the template %s', filename);
            elseif isempty(matching_files)
                error('No files match the template %s', filename);
            else 
                filename = fullfile(matching_files.folder , matching_files.name);
                eval(sprintf('%s = filename;',  filename_var ));
            end
        end

        if ~isfile(filename)
            error('File does not exist: \r\n%s', filename);
        end
    end

    %% LOAD AND PROCESS T1 IMAGE (OR CT IMAGE OR CONTRAST IMAGE)
    disp('Loading T1...');
    
    % Load T1-weighted MRI image and header information
    t1_image = niftiread(filename_t1);
    t1_header = niftiinfo(filename_t1);

%     % 2D T1 image: copy matrix to a third dimension
%     if length(size(t1_image)) == 2
%         t1_image = cat(3, t1_image, t1_image);
%         plot_options.slice_dim = 3;
%     else
%         plot_options.slice_dim = 2;
%     end

    % Determine transducer position based on Localite or configuration file
    if isfield(parameters,'transducer_from_localite') && parameters.transducer_from_localite
        [trans_pos_grid, focus_pos_grid, ~, ~] = ...
            position_transducer_localite(localite_file, t1_header, parameters);
    else
        trans_pos_grid = parameters.transducer.pos_t1_grid';
        focus_pos_grid = parameters.focus_pos_t1_grid';
    end
    
    if size(trans_pos_grid, 1) > size(trans_pos_grid, 2)
        trans_pos_grid = trans_pos_grid';
        focus_pos_grid = focus_pos_grid';
    end

    % Create unprocessed T1 slice oriented along the transducer's axis
%     h = figure;
%     imshowpair(plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), trans_pos_grid, focus_pos_grid, parameters, 'slice_dim', plot_options.slice_dim), ...
%                plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), trans_pos_grid, focus_pos_grid, parameters, 'slice_dim', 1), 'montage');
%     title('T1 with transducer');
%     output_plot_filename = fullfile(parameters.debug_dir, ...
%         sprintf('sub-%03d_t1_with_transducer_before_smoothing_and_cropping%s.png', ...
%         subject_id, parameters.results_filename_affix));
%     saveas(h, output_plot_filename, 'png');
%     close(h);

    %% SEGMENTATION USING SIMNIBS
    disp('Starting segmentation...');
    
    % Define output folder for segmentation results
    segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', subject_id));

    if strcmp(parameters.segmentation_software, 'charm')
        filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');
    else 
        filename_segmented = fullfile(segmentation_folder, sprintf('sub-%03d_final_contr.nii.gz', subject_id));
    end

    % Run segmentation if necessary
    if confirm_overwriting(filename_segmented, parameters) && ...
       (~isfield(parameters,'overwrite_simnibs') || parameters.overwrite_simnibs || ~exist(filename_segmented,'file'))
        if parameters.usepseudoCT == 1
            warning("SimNIBS integration not supported when requesting pseudoCT. Please ensure SimNIBS has been run.");
        end

        if parameters.interactive == 0 || confirmation_dlg('This will run SEGMENTATION WITH SIMNIBS that takes a long time. Are you sure?', 'Yes', 'No')
            run_segmentation(parameters.data_path, subject_id, filename_t1, filename_t2, parameters);
            return;
        end
    else
        disp('Skipping segmentation; loading existing file instead.');
    end
    
    %% Rotate to match the stimulation trajectory
    disp('Rotating to match the focus axis...')

    % If the headreco process was not successful, it will stop preprocessing
    assert(exist(filename_segmented,'file') > 0, ...
        'Head segmentation is not completed (%s does not exist), see logs in the batch_logs folder and in %s folder',...
            filename_segmented, segmentation_folder)

    % Defines output file location and name
    filename_reoriented_scaled_data = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_after_rotating_and_scaling%s.mat', ...
        subject_id, parameters.results_filename_affix));

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

    if confirm_overwriting(filename_reoriented_scaled_data, parameters)

        % Rotate and scale the original 3D T1 to line up with the transducer's axis
        if t1_header.ImageSize(3) > 1
            [t1_img_rr, ~, ~, ~, ~, ~, ~, t1_rr_img_montage] = ...
                align_to_focus_axis_and_scale(...
                t1_image, ...
                t1_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                t1_header.PixelDimensions(1)/parameters.grid_step_mm, ...
                parameters);
        else
            [t1_img_rr, ~, ~, ~, ~, t1_rr_img_montage] = ...
                align_to_focus_2d(...
                t1_image, ...
                t1_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                t1_header.PixelDimensions(1)/parameters.grid_step_mm, ...
                parameters);
        end
        
        h = figure;
        imshow(t1_rr_img_montage)
        title('Rotated (left) and original (right) original T1');
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_t1_after_rotating_and_scaling%s.png', ...
            subject_id, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png');
        close(h); clear t1_rr_img_montage;
        
        % Rotate and scale the segmented T1 (tissue mask) to line up with the transducer's axis
        [segmented_img_rr, trans_pos_upsampled_grid, focus_pos_upsampled_grid, ...
            scale_rotate_recenter_matrix, rotation_matrix, ~, ~, segm_img_montage] = ...
            align_to_focus_axis_and_scale(...
            tissues_mask_image, ...
            tissues_mask_header, ...
            trans_pos_grid, ...
            focus_pos_grid, ...
            tissues_mask_header.PixelDimensions(1)/parameters.grid_step_mm, ...
            parameters);
        
        h = figure;
        imshow(segm_img_montage)
        title('Rotated (left) and original (right) segmented T1');
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_segmented_after_rotating_and_scaling%s.png', ...
            subject_id, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png');
        close(h); clear segm_img_montage;
        
        % Rotate and scale the bone mask or pseudoCT to line up with the transducer's axis
        if parameters.usepseudoCT == 1
            [bone_img_rr, ~, ~, ~, ~, ~, ~, bone_img_montage] = ...
                align_to_focus_axis_and_scale(...
                pseudoCT_image, ...
                pseudoCT_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                pseudoCT_header.PixelDimensions(1)/parameters.grid_step_mm, ...
                parameters);
        else
            if strcmp(parameters.segmentation_software, 'charm') % create filled bone mask as charm doesn't make it itself
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
                align_to_focus_axis_and_scale(...
                bone_img, ...
                tissues_mask_header, ...
                trans_pos_grid, ...
                focus_pos_grid, ...
                tissues_mask_header.PixelDimensions(1)/parameters.grid_step_mm, ...
                parameters);
        end

        h = figure;
        imshow(bone_img_montage)
        title('Rotated (left) and original (right) original bone mask');
        output_plot_filename = fullfile(parameters.debug_dir, sprintf('sub-%03d_after_rotating_and_scaling_orig%s.png', ...
            subject_id, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h); clear bone_img_montage;

        assert(isequal(size(trans_pos_upsampled_grid(1:2)),size(focus_pos_upsampled_grid(1:2))),...
            "After reorientation, the first two coordinates of the focus and the transducer should be the same")

        % Saves the output according to the naming convention set in the
        % beginning of this section
        save(filename_reoriented_scaled_data, 'segmented_img_rr', 'trans_pos_upsampled_grid', 'bone_img_rr', 'focus_pos_upsampled_grid', 'scale_rotate_recenter_matrix', 'rotation_matrix', 't1_img_rr');
    else 
        disp('Skipping; the file already exists. Loading it instead.');
        load(filename_reoriented_scaled_data);
    end
    
    %% Plot the skin & skull from the segmented image
    
    % From SimNIBS charm segmentations
    % Segmentation labels can be different for the pseudoCT mask. See create_pseudoCT.sh.
    skullidx = getidx(parameters.layer_labels,'skull'); % Includes all "skull_" segmentations
    skinidx = getidx(parameters.layer_labels,'skin');

    % Unsmoothed skull & skin masks
    skull_mask_unsmoothed = ismember(segmented_img_rr,skullidx);
    skin_mask_unsmoothed = ismember(segmented_img_rr,skinidx);

    skin_slice = squeeze(skin_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:));
    skull_slice = squeeze(skull_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:));

    % Create a T1 slice for comparison to SimNIBS segmented data
    t1_slice = repmat(mat2gray(squeeze(t1_img_rr(:,trans_pos_upsampled_grid(2),:))), [1 1 3]);
    % Create a slice of segmented SimNIBS data
    skin_skull_img = cat(3, zeros(size(skull_slice)), skin_slice, skull_slice); % Skull is blue; skin is green.

    % Plot the different slices and an overlay for comparison
    h = figure;
    montage(cat(4,t1_slice*255 ,skin_skull_img*255 ,imfuse(mat2gray(t1_slice), skin_skull_img,'blend')) ,'size',[1 NaN]);
    title('T1 and SimNIBS skin (green) and skull (blue) masks');
    output_plot_filename = fullfile(parameters.debug_dir,sprintf('sub-%03d_t1_skin_skull%s.png',subject_id,...
                                parameters.results_filename_affix));
    saveas(h ,output_plot_filename ,'png')
    close(h);
    
    %% SMOOTH & CROP SKULL
    disp('Smoothing and cropping the skull...')

    % Defines output file location and name
    filename_cropped_smoothed_skull_data = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_after_cropping_and_smoothing%s.mat', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    % Uses one of three functions to crop and smooth the skull
    % See each respected function for more documentation
    if confirm_overwriting(filename_cropped_smoothed_skull_data, parameters)
        % postprocess skull segmentation (varies between skull, layered, and pseudoCT calls)
        [medium_masks, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, ~, ~, new_grid_size, crop_translation_matrix] = ...
            smooth_and_crop(...
            segmented_img_rr, bone_img_rr, parameters.grid_step_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters);

        % Combines the matrix that defined the alignment with the transducer
        % axis with the new matrix that defines the cropping of the skull
        final_transformation_matrix = scale_rotate_recenter_matrix*crop_translation_matrix';
        inv_final_transformation_matrix = maketform('affine', inv(final_transformation_matrix')');

        % save medium mask
        orig_hdr = t1_header; % header is based on original T1w (always present)
        orig_hdr.Datatype = 'single';
        segmented_file = fullfile(parameters.debug_dir, sprintf('sub-%03d_medium_masks_final', parameters.subject_id));
        plotdata = single(tformarray(uint8(medium_masks), inv_final_transformation_matrix, ...
            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
        if ~isfile(segmented_file)
            niftiwrite(plotdata, segmented_file, orig_hdr, 'Compressed',true);
        end
        clear segmented_file plotdata orig_hdr;

        % save skull mask/pseudoCT
        orig_hdr = t1_header; % header is based on original T1w (always present)
        orig_hdr.Datatype = 'single';
        skull_mask_file = fullfile(parameters.debug_dir, sprintf('sub-%03d_skull_final', parameters.subject_id));
        plotdata = single(tformarray(uint8(segmented_image_cropped), inv_final_transformation_matrix, ...
            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
        if ~isfile(skull_mask_file)
            niftiwrite(plotdata, skull_mask_file, orig_hdr, 'Compressed',true);
        end
        clear skull_mask_file plotdata orig_hdr;

        % Saves the output according to the naming convention set in the
        % beginning of this section
        save(filename_cropped_smoothed_skull_data, 'medium_masks', 'skull_edge', 'segmented_image_cropped', 'trans_pos_final', 'focus_pos_final', 'new_grid_size', 'crop_translation_matrix','final_transformation_matrix','inv_final_transformation_matrix')
    else 
        disp('Skipping, the file already exists, loading it instead.')
        load(filename_cropped_smoothed_skull_data);
    end    
    parameters.grid_dims = new_grid_size;

    % Creates and saves a figure with the segmented brain
    [seg_with_trans_img, transducer_pars] = plot_t1_with_transducer(...
        medium_masks, parameters.grid_step_mm, trans_pos_final, focus_pos_final, parameters);
    
    h = figure;
    imshow(seg_with_trans_img);
    title('Segmentation with transducer');
    output_plot_filename = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_segmented_brain_final%s.png', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'png')
    close(h);

    % Check that the transformations are correct by inverting them and
    % comparing to the original 
    if ~exist('inv_final_transformation_matrix','var')
        final_transformation_matrix = scale_rotate_recenter_matrix*crop_translation_matrix';
        inv_final_transformation_matrix = maketform('affine', inv(final_transformation_matrix')');
    end

    % If the transformation cannot be correctly inverted, this will be displayed
    backtransf_coordinates = round(tformfwd([trans_pos_final; focus_pos_final], inv_final_transformation_matrix));
    if ~all(all(backtransf_coordinates ==[trans_pos_grid; focus_pos_grid]))
        disp('Backtransformed focus and transducer parameters differ from the original ones. Something went wrong (but note that small rounding errors could be possible.')
        disp('Original coordinates')
        disp([trans_pos_final, focus_pos_final]')
        disp('Backtransformed coordinates')
        disp(backtransf_coordinates)
        exit()
    end
    
    % Plot positioning of transducer on segmentation
    output_plot_filename = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_positioning%s.png', ...
        subject_id, parameters.results_filename_affix));
    show_positioning_plots(...
        tissues_mask_image,...
        tissues_mask_header.PixelDimensions(1), ...
        trans_pos_grid, ...
        focus_pos_grid, ...
        medium_masks, ...
        trans_pos_final, ...
        focus_pos_final, ...
        parameters, ...
        output_plot_filename)

end