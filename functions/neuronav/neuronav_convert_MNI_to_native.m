function [native_target_mm, native_trans_mm, target_pos, trans_pos] = ...
    neuronav_convert_MNI_to_native(sub_id, parameters, pn, trans_mm, targ_mm)

% NEURONAV_CONVERT_MNI_TO_NATIVE  Transform coordinates from MNI space to native subject space
%
% Converts stimulation coordinates from MNI space back to native anatomical
% space via SimNIBS nonlinear warp fields, returning both RAS mm
% coordinates and voxel indices in the subject’s segmentation image.
%
% Use as:
%   [native_target_mm, native_trans_mm, target_pos, trans_pos] = ...
%       neuronav_convert_MNI_to_native(sub_id, parameters, pn, trans_mm, targ_mm)
%
% Input:
%   sub_id     - subject identifier string (e.g. ‘sub-010’)
%   parameters - (1,1) simulation parameters struct
%   pn         - (1,1) path names struct with data_seg field
%   trans_mm   - [Nx3] transducer coordinates in MNI space [mm RAS]
%   targ_mm    - [Nx3] target coordinates in MNI space [mm RAS]
%
% Output:
%   native_target_mm - [Nx3] target coordinates in native subject RAS space [mm]
%   native_trans_mm  - [Nx3] transducer coordinates in native subject RAS space [mm]
%   target_pos       - [Nx3] target voxel indices in segmentation image
%   trans_pos        - [Nx3] transducer voxel indices in segmentation image
%
% See also: NEURONAV_CONVERT_NATIVE_TO_MNI, LOCALITE_MATRIX_TO_POSITIONS

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