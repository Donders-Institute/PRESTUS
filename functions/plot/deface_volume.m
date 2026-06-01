function vol_out = deface_volume(vol, seg, voxel_size_mm, margin_mm)
% DEFACE_VOLUME  Remove face voxels in a volume using the SimNIBS segmentation
%
% Identifies the face as scalp/skin/muscle voxels anterior to a cut plane
% that is MARGIN_MM posterior to the posterior edge of the eye labels.
% The default margin (30 mm) places the cut plane at roughly the ear-canal
% level, matching the aggressiveness of standard defacing tools. All skull,
% brain, and CSF tissue is always preserved regardless of position.
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
%   margin_mm      - posterior margin past the back of the eye socket (default: 30)
%
% Output:
%   vol_out - copy of VOL with face voxels removed (NaN for float, 0 for integer)
%
% See also: CHARM_SEG_LABELS, PLOT_PLACEMENT_T1_OVERLAY, SHOW_3D_HEAD

    arguments
        vol            (:,:,:)
        seg            (:,:,:)
        voxel_size_mm  (1,1) double = 1.0
        margin_mm      (1,1) double = 30.0
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

    % Cut plane: posterior edge of eye socket + margin_mm, clamped to volume.
    % This places the cut at roughly the ear-canal level for a standard head.
    [~, y_eye, ~] = ind2sub(size(seg), find(eye_mask));
    margin_vox = round(margin_mm / voxel_size_mm);
    y_cutoff   = min(max(y_eye) + margin_vox, size(seg, 2));

    % Tissues that should never be removed regardless of position
    inside_head = ismember(seg, [seg_labels.bonemask, seg_labels.csf]);

    face_mask = false(size(seg));
    face_mask(:, 1:y_cutoff, :) = true;
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
