function [native_target_mm, native_trans_mm, target_pos, trans_pos] = ...
    neuronav_convert_MNI_to_native(sub_id, parameters, pn, trans_mm, targ_mm)

% NEURONAV_CONVERT_MNI_TO_NATIVE - Transform coordinates from MNI space to native subject space
%
% This function converts stimulation coordinates given in MNI space (millimeters)
% back into the native anatomical space of a subject, expressed both as continuous
% RAS millimeter coordinates and as voxel indices within the subject’s segmentation image.
%
% INPUTS:
%   sub_id      - Subject identifier string, e.g., 'sub-010'
%   parameters  - Struct containing subject-specific processing parameters and settings
%   pn          - Struct with standardized file paths and toolbox locations
%   trans_mm    - Nx3 array of transducer coordinates in MNI space (millimeters, RAS)
%   targ_mm     - Nx3 array of target coordinates in MNI space (millimeters, RAS)
%
% OUTPUTS:
%   native_target_mm - Nx3 array of target coordinates transformed to native subject RAS space (mm)
%   native_trans_mm  - Nx3 array of transducer coordinates transformed to native subject RAS space (mm)
%   target_pos       - Nx3 array of voxel indices corresponding to targets within the subject's segmentation image
%   trans_pos        - Nx3 array of voxel indices corresponding to transducers within the subject's segmentation image
%
% PROCESS:
%   1. Load segmentation NIfTI header info for the subject.
%   2. Use inverse nonlinear warp fields (SimNIBS mni2subject_coords_LDfix) to map MNI mm coords into native mm coords.
%   3. Convert native RAS mm coordinates to voxel indices in segmentation space via 'ras_to_grid'.
%
% This function is critical for relating group-level or atlas-derived MNI coordinates
% to the individual subject anatomy for precise spatial localization and accurate simulation.

    Npos = size(trans_mm,1); % Number of positions

    t1seg_info = niftiinfo(fullfile(pn.data_seg, sprintf('m2m_%s', sub_id), 'final_tissues.nii.gz'));
    m2m_folder = fullfile(pn.data_seg, sprintf('m2m_%s', sub_id));
    transformation_type = 'nonl';

    for i = 1:Npos
        % [1] inverse nonlinear deformation to warp MNI mm → native mm
        native_target_mm(i,:) = mni2subject_coords_LDfix(targ_mm(i,1:3), m2m_folder, parameters, transformation_type);
        native_trans_mm(i,:) = mni2subject_coords_LDfix(trans_mm(i,1:3), m2m_folder, parameters, transformation_type);
        % [2] from native mm convert to voxel location (in segmentation image)
        target_pos(i,:) = ras_to_grid(native_target_mm(i,:)', t1seg_info)';
        trans_pos(i,:) = ras_to_grid(native_trans_mm(i,:)', t1seg_info)';
    end

end