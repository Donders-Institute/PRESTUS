function [trans_ras_seg, ...
    target_ras_seg, ...
    trans_mni_pos, ...
    target_mni_pos, ...
    trans_mni_ras, ...
    targ_mni_ras] = ...
    neuronav_convert_native_to_MNI...
    (sub_id, parameters, pn, trans_ras, target_ras, side_order)

% NEURONAV_CONVERT_NATIVE_TO_MNI - Convert native RAS coordinates to MNI space
%
% This function transforms transducer and target coordinates from native subject 
% space (RAS mm) into segmentation space and then into standardized MNI space,
% producing both voxel indices and RAS coordinates in MNI space. It also clamps 
% coordinates to the valid overlapping spatial bounds between segmentation and MNI spaces 
% to ensure anatomical plausibility. Optionally, the output voxel positions are adjusted 
% based on the recorded stimulation hemisphere ('left' or 'right').
%
% INPUTS:
%   sub_id       - Subject ID string, e.g., 'sub-010'
%   parameters   - Struct containing subject-specific parameters and settings
%   pn           - Struct with paths to necessary data and toolboxes
%   trans_ras    - Nx3 array of transducer RAS coordinates in native space (in mm)
%   target_ras   - Nx3 array of target RAS coordinates in native space (in mm)
%   side_order   - Cell array of strings {'left','right',...} indicating laterality per coordinate (optional)
%
% OUTPUTS:
%   trans_ras_seg   - Nx3 array of transducer coordinates converted to segmentation RAS space (mm)
%   target_ras_seg  - Nx3 array of target coordinates converted to segmentation RAS space (mm)
%   trans_mni_pos   - Nx3 array of transducer voxel indices in MNI space (integer)
%   target_mni_pos  - Nx3 array of target voxel indices in MNI space (integer)
%   trans_mni_ras   - Nx3 array of transducer continuous RAS coordinates in MNI space (mm)
%   targ_mni_ras    - Nx3 array of target continuous RAS coordinates in MNI space (mm)
%
% STEPS:
%   1. Load image headers for planning, segmentation, and MNI reference spaces.
%   2. Load nonlinear warp field transforming native segmentation to MNI.
%   3. Compute valid coordinate boundaries by mapping MNI voxel bounds back to segmentation space.
%   4. Transform input native RAS points (transducer and target) to segmentation RAS coordinates.
%   5. Clamp segmentation coordinates within computed spatial overlap bounds.
%   6. Use SimNIBS helper functions to map segmentation RAS to MNI RAS space.
%   7. Convert MNI RAS coordinates to voxel indices using image affine info.
%   8. Adjust voxel indices laterally based on hemisphere side to compensate for limited MNI field-of-view.
%   9. Convert adjusted voxel indices back to MNI RAS mm coordinates.
%
% NOTES:
%   - Input and output coordinates use the RAS convention (Right-Anterior-Superior)
%   - Voxel indices are 1-based MATLAB indexing (x, y, z)
%   - side_order input influences lateral voxel adjustments to correct for partial MNI coverage
%   - Coordinates reported in multiple spaces enable precise spatial localization and group comparisons
%
% DEPENDENCIES:
%   - niftiinfo, niftiread (MATLAB NIfTI toolbox)
%   - mni2subject_coords_LDfix.m (SimNIBS coordinate transformation)
%   - ras_to_grid.m (convert RAS mm to voxel indices)

    % Number of positions
    Npos = size(trans_ras,1);

    % Convert RAS locations to MNI space (RAS mm) & (voxel)

    % Note: the planning and segmentation images have the same size,
    % but with different headers / location systems

    t1plan_info = niftiinfo(fullfile(pn.data_prelocalite, sprintf('%s_T1_forneuronav.nii.gz',sub_id)));
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