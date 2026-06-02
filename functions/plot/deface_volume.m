function vol_out = deface_volume(vol, seg, voxel_size_mm, margin_mm)
% DEFACE_VOLUME  Remove face voxels in a volume using the SimNIBS segmentation
%
% Identifies the face as scalp/skin/muscle voxels anterior to a cut plane
% that is MARGIN_MM posterior to the posterior edge of the eye labels.
% The default margin (10 mm) places the cut plane just behind the eye socket,
% removing the nose/mouth region while preserving the ears and most of the
% lateral face. All skull, brain, and CSF tissue is always preserved
% regardless of position.
%
% For floating-point volumes (T1 intensity) face voxels are set to NaN so
% that plot functions can render them as transparent rather than black.
% For integer volumes (tissue label maps) face voxels are set to 0 (background
% label) so that isosurface extraction simply omits them.
%
% Use as:
%   vol_defaced = deface_volume(vol, seg)
%   vol_defaced = deface_volume(vol, seg, voxel_size_mm)
%   vol_defaced = deface_volume(vol, seg, voxel_size_mm, margin_mm)
%
% Input:
%   vol            - [Nx x Ny x Nz] image volume to deface (any numeric type)
%   seg            - [Nx x Ny x Nz] integer SimNIBS charm tissue label volume,
%                    same size as VOL
%   voxel_size_mm  - isotropic voxel size in mm (default: 1.0)
%   margin_mm      - posterior margin past the back of the eye socket (default: 10)
%
% Output:
%   vol_out - copy of VOL with face voxels removed (NaN for float, 0 for integer)
%
% See also: CHARM_SEG_LABELS, PLOT_PLACEMENT_T1_OVERLAY, SHOW_3D_HEAD

    arguments
        vol            (:,:,:)
        seg            (:,:,:)
        voxel_size_mm  (1,1) double = 1.0
        margin_mm      (1,1) double = 10.0
    end

    if ~isequal(size(vol), size(seg))
        error('deface_volume: VOL and SEG must be the same size.');
    end

    seg_labels = charm_seg_labels();

    % Locate eye voxels to anchor the face plane
    eye_mask = (seg == seg_labels.eye);
    if ~any(eye_mask(:))
        warning('deface_volume: no eye voxels found in segmentation; defacing skipped.');
        vol_out = vol;
        return;
    end

    % Determine which Y direction is anterior by comparing the eye centroid
    % to the brain centroid. Eyes are always anterior to the brain, so this
    % is orientation-agnostic and works for both RAS and LAS/other NIfTIs.
    [~, y_eye, ~] = ind2sub(size(seg), find(eye_mask));
    brain_mask = ismember(seg, [seg_labels.wm, seg_labels.gm]);
    [~, y_brain, ~] = ind2sub(size(seg), find(brain_mask));
    anterior_is_high_y = mean(y_eye) > mean(y_brain);

    % Cut plane: posterior edge of eye socket shifted margin_mm into the head,
    % placing the cut at roughly the ear-canal level.
    margin_vox = round(margin_mm / voxel_size_mm);
    if anterior_is_high_y
        % Face at high Y — mask from cut plane upward
        y_cutoff = max(min(y_eye) - margin_vox, 1);
        face_mask = false(size(seg));
        face_mask(:, y_cutoff:end, :) = true;
    else
        % Face at low Y — mask from cut plane downward
        y_cutoff = min(max(y_eye) + margin_vox, size(seg, 2));
        face_mask = false(size(seg));
        face_mask(:, 1:y_cutoff, :) = true;
    end

    % Tissues that should never be removed regardless of position
    inside_head = ismember(seg, [seg_labels.bonemask, seg_labels.csf]);
    face_mask = face_mask & ~inside_head;

    vol_out = vol;
    if isinteger(vol)
        % Integer label maps: 0 is the background label — isosurface ignores it
        vol_out(face_mask) = 0;
    else
        % Float intensity images: NaN lets plot functions render the region
        % as transparent rather than black
        vol_out = double(vol_out);
        vol_out(face_mask) = NaN;
    end
end
