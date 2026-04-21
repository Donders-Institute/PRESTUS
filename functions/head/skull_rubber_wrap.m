function [SKULL_BALLON] = skull_rubber_wrap(parameters, BW, medium_masks, segmented_img)
% SKULL_RUBBER_WRAP  Fill local holes in the skull layer using a balloon inflation approach
%
% Iteratively inflates a balloon from the skin surface inward to fill gaps
% in the binary skull mask, producing a watertight skull layer for
% simulation. Saves the result as a NIfTI file.
%
% Use as:
%   SKULL_BALLON = skull_rubber_wrap(parameters, BW, medium_masks, segmented_img)
%
% Input:
%   parameters    - (1,1) simulation configuration struct
%   BW            - [Nx x Ny x Nz] binary skull mask
%   medium_masks  - [Nx x Ny x Nz] medium label array
%   segmented_img - [Nx x Ny x Nz] SimNIBS tissue label volume
%
% Output:
%   SKULL_BALLON - [Nx x Ny x Nz] filled binary skull mask
%
% See also: SKULL_FILL_HOLES, PCT_SKULLEXPAND

    disp('Using rubber wrap to fill potential local holes in skull layer...')

    % ==============================================================================
    % Nifti info (to save new niftis)
    % ==============================================================================
    
    % Load the T1w input image as a reference
    info = niftiinfo(fullfile(parameters.path.seg, parameters.path.t1_pattern));

    % ==============================================================================
    % CLEAN SKULL MASK
    % ==============================================================================
    
    % Keep the largest connected skull component
    
    sz = size(BW);
    CC = bwconncomp(BW, 26);
    if CC.NumObjects == 0, error("Skull mask is empty."); end
    sizes = cellfun(@numel, CC.PixelIdxList);
    [~, iBig] = max(sizes);
    BW1 = false(sz);
    BW1(CC.PixelIdxList{iBig}) = true;

    % ==============================================================================
    % RUBBER WRAP AROUND SKULL = BALLOON MASK
    % ==============================================================================

    %% --- Rubber wrap parameters ---

    if isfield(parameters.headmodel, 'skull_wrap_radius') && ~isempty(parameters.headmodel.skull_wrap_radius)
        wrapRadius = parameters.headmodel.skull_wrap_radius;
    else
        wrapRadius = 10; % skull rubber wrap radius [grid voxels];
    end

    se = strel("sphere", wrapRadius);

    % "Shrink-wrap" / morphological hull:
    % Bridge concavities that are smaller than wrapRadius (-> rubber won't go in).
    BWwrap = imclose(BW1, se);

    % Ensure that the wrap is a continuous, filled region
    BWwrap = imfill(BWwrap, "holes");

    % Optional: remove any stray bits and keep largest component again
    CCw = bwconncomp(BWwrap, 26);
    sizesw = cellfun(@numel, CCw.PixelIdxList);
    [~, iw] = max(sizesw);
    BWwrap2 = false(sz);
    BWwrap2(CCw.PixelIdxList{iw}) = true;

    % Choose balloon mask, ensure logical / uint8
    BW = uint8(BWwrap2); clear BWwrap2;

    % [debug] save balloon wrap mask
    if parameters.simulation.debug == 1

        outNii = fullfile(parameters.io.debug_dir_preproc,...
            sprintf('balloon_mask%s.nii', parameters.io.output_affix));

        infoOut = info;
        infoOut.ImageSize = size(BW);
        infoOut.PixelDimensions = repmat(parameters.grid.resolution_mm,1,numel(size(BW)));
        infoOut.Datatype = 'uint8';
        infoOut.BitsPerPixel = 8;
        infoOut.Filename = outNii;
        infoOut.Filemoddate = char(datetime("now"));

        niftiwrite(BW, outNii, infoOut, 'Compressed', false);

        gzip(outNii);
        delete(outNii);   % remove uncompressed .nii
    end

    %% Manage layered head tissues

    % Note: the final_tissue outputs are more fine-grained than the layered masks inside of PRESTUS.
    % Currently, we use the latter.

    T = medium_masks;

    if ~isequal(size(BW), size(T))
        error("Size mismatch: balloon %s vs tissues %s", mat2str(size(BW)), mat2str(size(T)));
    end

    tissues_available = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium_properties);

    if ismember(tissues_available, 'brain')
        BRAIN_LABEL   = find(strcmp(medium_labels, 'brain'));
        BRAIN = (T == BRAIN_LABEL);
    else
        BRAIN = zeros(size(T));
    end
    
    if ismember(tissues_available, 'brain')
        SKIN_LABEL = find(strcmp(medium_labels, 'skin'));
        SKIN = (T == SKIN_LABEL);
    else
        SKIN = zeros(size(T));
    end

    %% "Touch" masks: balloon voxel touches BRAIN if it lies within 1-voxel dilation of BRAIN

    touchRadius = 1;   % voxels: 1 = immediate neighbors (26-neighborhood)

    se = strel("sphere", touchRadius);

    GM_touch_region   = imdilate(BRAIN,   se);
    SKIN_touch_region = imdilate(SKIN, se);

    % Keep balloon voxels that touch BOTH
    B_touch_both =  BW & GM_touch_region & SKIN_touch_region;

    % Remove overlap with BRAIN from the remaining mask
    B_out = B_touch_both & ~BRAIN;

    %% Include CSF voxels at the SKIN interface (strict 6-neigh), near B_out

    % Here, we refer to the more detailed segmentation
    CSF_LABEL = charm_seg_labels().csf;
    CSF  = (ismember(segmented_img, CSF_LABEL));

    % --- Strict 6-neighborhood kernel (faces only) ---
    K6 = zeros(3,3,3,'logical');
    K6(2,2,2) = true;
    K6(1,2,2) = true; K6(3,2,2) = true;
    K6(2,1,2) = true; K6(2,3,2) = true;
    K6(2,2,1) = true; K6(2,2,3) = true;

    % SKIN dilated by 6-neighborhood => voxels face-adjacent to SKIN
    SKIN_touch_region_6 = convn(single(SKIN), single(K6), 'same') > 0;

    % CSF voxels that touch SKIN (faces only)
    CSF_touch_SKIN_6 = CSF & SKIN_touch_region_6;

    % --- Restrict to CSF voxels near existing balloon mask ---
    nearRadius = 2;                 % 1–2 suggested
    seNear = strel("sphere", nearRadius);
    nearBalloon = imdilate(B_out, seNear);

    CSF_touch_SKIN_6_near = CSF_touch_SKIN_6 & nearBalloon;

    % fprintf("CSF voxels touching SKIN (6-neigh): %d\n", nnz(CSF_touch_SKIN_6));
    % fprintf("...of those, near existing B_out (r=%d): %d\n", nearRadius, nnz(CSF_touch_SKIN_6_near));

    % --- Add MEN interface voxels to balloon mask ---
    B_out2 = B_out | CSF_touch_SKIN_6_near;

    % Keep your original constraint: no BRAIN overlap
    B_out2 = B_out2 & ~BRAIN;

    B_out = BW1 | B_out2;

    %% Fill along Z where there is skull in front AND behind (same Y,X)

    B0 = (B_out ~= 0);   % canonical logical

    BWm1 = circshift(B0, [0 0  1]);
    BWp1 = circshift(B0, [0 0 -1]);

    BWm1(:,:,1)   = false;
    BWp1(:,:,end) = false;

    %% Collect final outputs

    SKULL = BW1;
    BALLOON = B0 & ~SKULL;      % balloon-approach voxels
    ZADDED = ~B0 & BWm1 & BWp1;
    SKULL_BALLON = B0 | ZADDED;   % final filled ballon mask (logical)

    %% [DEBUG] Save as .nii.gz

    if parameters.simulation.debug == 1

        outNii = fullfile(parameters.io.debug_dir_preproc,...
            sprintf('balloon_mask_final%s.nii', parameters.io.output_affix));

        infoOut = info;
        infoOut.ImageSize = size(SKULL_BALLON);
        infoOut.PixelDimensions = repmat(parameters.grid.resolution_mm,1,numel(size(SKULL_BALLON)));
        infoOut.Datatype = 'uint8';
        infoOut.BitsPerPixel = 8;
        infoOut.Filename = outNii;
        infoOut.Filemoddate = char(datetime("now"));

        niftiwrite(uint8(SKULL_BALLON), outNii, infoOut, 'Compressed', false);
        gzip(outNii);
        delete(outNii);

        % fprintf("Saved: %s.gz\n", outNii);
    end

    %% [DEBUG] Visualize the expanded skull
    % Note: this seems computationally heavy on the HPC
    % There are like;ly ways to optimize, but it will be inactive for now.
    
    if parameters.headmodel.skull_wrap_visualize == 1
        skull_rubber_wrap_visualize(parameters, SKULL, ZADDED, BALLOON);
    end

end