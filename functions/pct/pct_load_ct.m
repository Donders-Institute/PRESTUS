function [ct_image, ct_header, mask_image, mask_header, ct_file] = pct_load_ct(parameters)
% PCT_LOAD_CT  The CT image the skull mapping reads, with its tissue labels
%
% The skull mapping (pct.mapping_*) reads one CT image in Hounsfield units:
%   - a measured CT (path.ct_pattern set): the CT is resampled (pct.interpolation) onto
%     the SimNIBS segmentation grid through both images' world (RAS mm)
%     affines, so it must already be co-registered to the subject's T1w.
%     Voxels outside the CT are air (-1000 HU). Tissue labels come from the
%     segmentation (PCT_TISSUES_MASK). Nothing is written to disk.
%   - otherwise the pseudo-CT generated from the UTE image:
%     pseudoCT.nii.gz and tissues_mask.nii.gz in the pCT folder
%     (FIND_PSEUDOCT_DIR).
% Images are returned as stored on disk; callers apply ENSURE_RAS_PLUS.
%
% Use as:
%   [ct_image, ct_header, mask_image, mask_header, ct_file] = pct_load_ct(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (path.ct_pattern, path.anat,
%                path.seg, path.pct, subject_id)
%
% Output:
%   ct_image    - [Nx x Ny x Nz] HU image on the segmentation grid
%   ct_header   - its niftiinfo header (the segmentation's grid)
%   mask_image  - tissue labels on the same grid (bone = 4)
%   mask_header - header of mask_image
%   ct_file     - the CT file that was read
%
% See also: PCT_CT_FILE, PCT_TISSUES_MASK, FIND_PSEUDOCT_DIR, PREPROC_HEAD

arguments
    parameters (1,1) struct
end

    ct_file = pct_ct_file(parameters);
    if isempty(ct_file)
        pct_dir = find_pseudoct_dir(parameters);
        if isempty(pct_dir)
            error('PRESTUS:pct:notFound', ...
                ['No CT for the skull mapping: path.ct_pattern is empty and no pseudoCT.nii.gz ' ...
                 'was found in:\n  %s'], get_pct_dir(parameters));
        end
        ct_file     = fullfile(pct_dir, 'pseudoCT.nii.gz');
        ct_image    = niftiread(ct_file);
        ct_header   = niftiinfo(ct_file);
        mask_file   = fullfile(pct_dir, 'tissues_mask.nii.gz');
        mask_image  = niftiread(mask_file);
        mask_header = niftiinfo(mask_file);
        return
    end

    seg_file = fullfile(char(parameters.path.seg), sprintf('m2m_sub-%03d', parameters.subject_id), ...
        'final_tissues.nii.gz');
    if ~isfile(seg_file)
        error('PRESTUS:pct:noSegmentation', 'A measured CT needs the segmentation: %s not found.', seg_file);
    end
    mask_header = niftiinfo(seg_file);
    mask_image  = pct_tissues_mask(niftiread(seg_file));

    ct_hdr = niftiinfo(ct_file);
    % niftiread returns stored values; apply the NIfTI scaling (slope 0 = none) to get HU
    ct_in = double(niftiread(ct_file));
    if ~isempty(ct_hdr.MultiplicativeScaling) && ct_hdr.MultiplicativeScaling ~= 0
        ct_in = ct_in * ct_hdr.MultiplicativeScaling + ct_hdr.AdditiveOffset;
    end

    % 1-based CT voxel -> world -> 1-based segmentation voxel, row-vector
    % convention ([vox0 1] * T = [ras 1]) as AFFINE_RESAMPLE_3D expects
    to0 = eye(4); to0(4, 1:3) = -1;
    to1 = eye(4); to1(4, 1:3) =  1;
    M = to0 * ct_hdr.Transform.T / mask_header.Transform.T * to1;
    seg_size = mask_header.ImageSize(1:3);
    ct_image = single(affine_resample_3d(ct_in, M, seg_size, ...
        char(get_field(parameters, {'pct', 'interpolation'}, 'linear')), -1000));

    ct_header = mask_header;
    ct_header.Datatype = 'single';
    ct_header.BitsPerPixel = 32;
    ct_header.MultiplicativeScaling = 1;
    ct_header.AdditiveOffset = 0;

    report_overlap(ct_image, mask_image, ct_file);
end

function report_overlap(ct_image, mask_image, ct_file)
% A CT that is not co-registered to the T1 puts air or soft tissue in the skull
    bone = mask_image == 4;
    if ~any(bone(:))
        return
    end
    hu = sort(double(ct_image(bone)));
    q = @(p) hu(max(1, round(p * numel(hu))));     % no Statistics Toolbox (prctile)
    fprintf('[pct_load_ct] Measured CT %s: skull voxels median %.0f HU (5-95%%: %.0f to %.0f HU)\n', ...
        ct_file, median(hu), q(0.05), q(0.95));
    if median(hu) < 300
        warning('PRESTUS:pct:ctMisregistered', ...
            ['The skull voxels of the segmentation have a median of %.0f HU in the measured CT %s. ' ...
             'Bone is usually several hundred HU or more: check that the CT is co-registered ' ...
             'to the T1w (same world coordinates) and in Hounsfield units.'], median(hu), ct_file);
    end
end
