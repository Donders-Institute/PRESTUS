function tissues_mask = pct_tissues_mask(segmentation)
% PCT_TISSUES_MASK  Tissue labels for a CT-mapped skull, from the SimNIBS segmentation
%
% The MATLAB counterpart of the tissues_mask.nii.gz step of
% pct_create_pseudoCT.sh, used for a measured CT (path.ct_pattern), where
% that script does not run:
%   1. skull     = charm labels 7 (compact) and 8 (spongy bone)
%   2. interface = voxels outside the head within one voxel of the skull
%                  (skull-air partial volume)
%   3. the skull and interface, dilated by one voxel (3x3x3 box, as
%      fslmaths -dilF), are labelled 4 (bone); other voxels keep their
%      charm label
%
% Use as:
%   tissues_mask = pct_tissues_mask(segmentation)
%
% Input:
%   segmentation - [Nx x Ny x Nz] charm final_tissues labels
%
% Output:
%   tissues_mask - labels of the same size and class, bone = 4
%
% See also: PCT_LOAD_CT, PCT_CREATE_PCT_RUN

arguments
    segmentation {mustBeNumeric}
end

    box = ones(3, 3, 3, 'single');
    dilate = @(m) convn(single(m), box, 'same') > 0;

    skull     = segmentation == 7 | segmentation == 8;
    interface = dilate(skull) & segmentation == 0;
    bone      = dilate(skull | interface);

    tissues_mask = segmentation;
    tissues_mask(bone) = 4;
end
