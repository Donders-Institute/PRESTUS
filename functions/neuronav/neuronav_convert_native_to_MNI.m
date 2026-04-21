function [trans_ras_seg, ...
    target_ras_seg, ...
    trans_mni_pos, ...
    target_mni_pos, ...
    trans_mni_ras, ...
    targ_mni_ras] = ...
    neuronav_convert_native_to_MNI...
    (sub_id, parameters, pn, trans_ras, target_ras, side_order)

% NEURONAV_CONVERT_NATIVE_TO_MNI  Convert native RAS coordinates to MNI space
%
% Transforms transducer and target coordinates from native subject RAS
% space to standardised MNI space via SimNIBS nonlinear warp fields.
% Clamps coordinates to valid overlapping bounds and optionally adjusts
% for stimulation laterality.
%
% Use as:
%   [trans_ras_seg, target_ras_seg, trans_mni_pos, target_mni_pos, ...
%    trans_mni_ras, targ_mni_ras] = neuronav_convert_native_to_MNI( ...
%       sub_id, parameters, pn, trans_ras, target_ras)
%   [...] = neuronav_convert_native_to_MNI(..., side_order)
%
% Input:
%   sub_id     - subject identifier string (e.g. 'sub-010')
%   parameters - (1,1) simulation parameters struct
%   pn         - (1,1) path names struct
%   trans_ras  - [Nx3] transducer RAS coordinates in native space [mm]
%   target_ras - [Nx3] target RAS coordinates in native space [mm]
%   side_order - cell array of laterality strings {'left','right',...} (optional)
%
% Output:
%   trans_ras_seg  - [Nx3] transducer coordinates in segmentation RAS space [mm]
%   target_ras_seg - [Nx3] target coordinates in segmentation RAS space [mm]
%   trans_mni_pos  - [Nx3] transducer voxel indices in MNI space
%   target_mni_pos - [Nx3] target voxel indices in MNI space
%   trans_mni_ras  - [Nx3] transducer RAS coordinates in MNI space [mm]
%   targ_mni_ras   - [Nx3] target RAS coordinates in MNI space [mm]
%
% See also: NEURONAV_CONVERT_MNI_TO_NATIVE, NEURONAV_GET_GROUP_MEAN_MNI

    % Number of positions
    Npos = size(trans_ras,1);

    % Convert RAS locations to MNI space (RAS mm) & (voxel)

    % Note: the planning and segmentation images have the same size,
    % but with different headers / location systems

    t1plan_info = niftiinfo(fullfile(pn.data_prelocalite, sprintf('%s_T1*.nii*',sub_id)));
    t1seg_info = niftiinfo(fullfile(pn.data_seg, sprintf('m2m_%s', sub_id), 'final_tissues.nii.gz'));
    mni_info = niftiinfo(fullfile(pn.data_seg, sprintf('m2m_%s', sub_id), 'toMNI', 'final_tissues_MNI.nii.gz'));
    
    warp_file = fullfile(pn.data_seg, sprintf('m2m_%s', sub_id), 'toMNI', 'Conform2MNI_nonl.nii.gz');
    warp_field = niftiread(warp_file);     % [X Y Z 3] nonlinear field: native→MNI displacement (mm)
    warp_info  = niftiinfo(warp_file);     % Contains the affine
    
    plan_affine = t1plan_info.Transform.T;
    seg_affine = t1seg_info.Transform.T;
    warp_affine = warp_info.Transform.T;
    mni_affine = mni_info.Transform.T;

    m2m_folder = fullfile(pn.data_seg, sprintf('m2m_%s', sub_id));
    transformation_type = 'nonl';

    % Intermission: ensure that segmentation coordinates map onto MNI
    % Use the following logic: identify segmentation space mm bounds
    % based on maximum MNI voxels; clamp segmentation coordinates by
    % those bounds.

    % --- 1. Compute min/max voxel for the target MNI image ---
    min_voxel_mni = [1 1 1 1]';                   % 1-based indexing (MATLAB convention)
    max_voxel_mni = [mni_info.ImageSize(1:3),1]';      % Make sure to make it a column 4-vector
    
    % --- 2. Convert voxel bounds to RAS (mm) bounds in MNI space ---
    min_ras_mni = mni_affine' * min_voxel_mni;
    max_ras_mni = mni_affine' * max_voxel_mni;
    
    % --- 3. Convert those RAS bounds into segmentation RAS space ---
    
    min_ras_seg = mni2subject_coords_LDfix(min_ras_mni(1:3)', m2m_folder, parameters, transformation_type);
    max_ras_seg = mni2subject_coords_LDfix(max_ras_mni(1:3)', m2m_folder, parameters, transformation_type);
    
    % Take the (x,y,z) coordinates:
    min_ras_seg_coords = min(min_ras_seg(1:3), max_ras_seg(1:3));  % Account for potential axis flip
    max_ras_seg_coords = max(min_ras_seg(1:3), max_ras_seg(1:3));
    
    % --- 4. Regularize/clamp all segmentation-space coordinates to that intersection ---
    clip_to_bounds = @(coord, min_val, max_val) ...
        max(min(coord, max_val), min_val);
    
    % This ensures trans_ras_seg are always valid for subsequent mapping to MNI and stay inside the overlap volume.

    for i = 1:max(Npos)
        if any(isnan(trans_ras(i,:))) || any(isnan(target_ras(i,:)))
            warning('Side %d: native coordinate missing, skipping warp.', i);
            continue;
        end

        % [1] transform from planning voxel -> RAS -> segmentation RAS

        vox_plan = [target_ras(i,1:3) 1] * inv(plan_affine);            % Planning voxel (not rounded) = target_pos(i,:)
        R = seg_affine(1:3,1:3);
        T = seg_affine(4,1:3)';
        ras_seg = (vox_plan(1:3) * R) + T';
        target_ras_seg(i,:) = ras_seg(1:3);                             % RAS mm in segmentation space

        vox_plan = [trans_ras(i,1:3) 1] * inv(plan_affine);
        R = seg_affine(1:3,1:3);
        T = seg_affine(4,1:3)';
        ras_seg = (vox_plan(1:3) * R) + T';
        trans_ras_seg(i,:) = ras_seg(1:3);

        % clamp transducer values to min/max MNI space
        orig_coord = trans_ras_seg(i,1:3);
        regularized_coord = clip_to_bounds(orig_coord, min_ras_seg_coords, max_ras_seg_coords);
        trans_ras_seg(i,:) = regularized_coord;

        % [2] transform from warp RAS to MNI RAS (use dedicated SimNIBS function)
        % Note: the FOV of the MNI image seems very limited for the
        % transducer location. So we'll adjust it laterally until it
        % fits.

        targ_mni_mm = subject2mni_coords_LDfix(target_ras_seg(i,1:3), m2m_folder, parameters, transformation_type);
        trans_mni_mm = [0 0 0];
        trans_mni_mm = subject2mni_coords_LDfix(trans_ras_seg(i,1:3), m2m_folder, parameters, transformation_type);
        while all(trans_mni_mm==0)
            warning('MNI space too contrained for original segmentation-space locations ... adjusting')
            trans_mni_mm = subject2mni_coords_LDfix(trans_ras_seg(i,1:3), m2m_folder, parameters, transformation_type);
            if trans_ras_seg(i,1)>0
                trans_ras_seg(i,1) = trans_ras_seg(i,1)-1; % change the right-left location to fall within the frame
            elseif trans_ras_seg(i,1)<0
                trans_ras_seg(i,1) = trans_ras_seg(i,1)+1;
            end
        end
        % alternative: try to transform manually with matrices
%             targ_mni_mm  = neuronav_apply_deformation(target_ras_seg(1:3), warp_field, warp_affine);
%             trans_mni_mm = neuronav_apply_deformation(trans_ras_seg(1:3), warp_field, warp_affine);
        targ_mni_ras(i,:)  = round(targ_mni_mm);
        trans_mni_ras(i,:) = round(trans_mni_mm);

        % [3] transform MNI from RAS to voxel grid
        target_mni_pos(i,:) = ras_to_grid(targ_mni_mm', mni_info)';
        trans_mni_pos(i,:) = ras_to_grid(trans_mni_mm', mni_info)';

        % [4] shift MNI voxel right-left to min/max of MNI space
        % this tries to deal with the relatively cropped out-of-skull window
        % depends on optional input side_order

        % individual image is RAS (standard): smaller values are left
        % important: MNI is LAS (smaller values are right)
        if exist('side_order')
            if strcmp(side_order{i}, 'left')
                trans_mni_pos(i,1) = max_voxel_mni(1);
            elseif strcmp(side_order{i}, 'right')
                trans_mni_pos(i,1) = min_voxel_mni(1);
            end
        end
        
        % [5] Convert clamped voxel positions back to RAS mm
        tmp  = (mni_affine' * [trans_mni_pos(i,:) 1]')';
        trans_mni_ras(i,:)  = tmp(1:3); clear tmp;
    end
end