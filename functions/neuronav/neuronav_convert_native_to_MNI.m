function [trans_ras_seg, ...
    target_ras_seg, ...
    trans_mni_pos, ...
    target_mni_pos, ...
    trans_mni_ras, ...
    targ_mni_ras] = ...
    neuronav_convert_native_to_MNI...
    (sub_id, parameters, trans_ras, target_ras, side_order)
% NEURONAV_CONVERT_NATIVE_TO_MNI  Convert native RAS coordinates to MNI space
%
% Transforms transducer and target coordinates from native subject RAS
% space to standardised MNI space by sampling SimNIBS's nonlinear
% coordinate field directly (see NEURONAV_APPLY_DEFORMATION), redirecting
% positions outside its defined domain to the nearest valid voxel, and
% optionally adjusts for stimulation laterality.
%
% Use as:
%   [trans_ras_seg, target_ras_seg, trans_mni_pos, target_mni_pos, ...
%    trans_mni_ras, targ_mni_ras] = neuronav_convert_native_to_MNI( ...
%       sub_id, parameters, trans_ras, target_ras)
%   [...] = neuronav_convert_native_to_MNI(..., side_order)
%
% Input:
%   sub_id     - subject identifier string (e.g. 'sub-010')
%   parameters - PRESTUS config struct: path.localite_t1 (or localite_raw)
%                for the planning T1 (NEURONAV_FIND_PLANNING_T1), and path.seg
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
% See also: NEURONAV_APPLY_DEFORMATION, NEURONAV_BUILD_DEFORMATION_CACHE,
%           NEURONAV_CONVERT_MNI_TO_NATIVE, NEURONAV_GET_GROUP_MEAN_MNI

    % Number of positions
    Npos = size(trans_ras,1);

    % Convert RAS locations to MNI space (RAS mm) & (voxel)

    % Note: the planning and segmentation images have the same size,
    % but with different headers / location systems

    t1plan_info = niftiinfo(neuronav_find_planning_t1(sub_id, parameters));
    t1seg_info = niftiinfo(fullfile(parameters.path.seg, sprintf('m2m_%s', sub_id), 'final_tissues.nii.gz'));
    mni_info = niftiinfo(fullfile(parameters.path.seg, sprintf('m2m_%s', sub_id), 'toMNI', 'final_tissues_MNI.nii.gz'));

    warp_file = fullfile(parameters.path.seg, sprintf('m2m_%s', sub_id), 'toMNI', 'Conform2MNI_nonl.nii.gz');
    warp_field = niftiread(warp_file);     % [X Y Z 3] nonlinear field: native voxel -> MNI RAS coordinate (mm)
    warp_info  = niftiinfo(warp_file);     % Contains the affine

    plan_affine = t1plan_info.Transform.T;
    seg_affine = t1seg_info.Transform.T;
    warp_affine = warp_info.Transform.T;

    % --- Min/max voxel for the target MNI image (used for the side_order clamp below) ---
    min_voxel_mni = [1 1 1 1]';                   % 1-based indexing (MATLAB convention)
    max_voxel_mni = [mni_info.ImageSize(1:3),1]';      % Make sure to make it a column 4-vector

    % Nearest-valid-voxel lookup for the native->MNI coordinate field, so
    % positions falling outside its defined domain resolve directly
    % instead of via a per-point CLI retry.
    nearest_idx = neuronav_build_deformation_cache(warp_field);

    for i = 1:max(Npos)
        if any(isnan(trans_ras(i,:))) || any(isnan(target_ras(i,:)))
            warn('Side %d: native coordinate missing, skipping warp.', i);
            continue;
        end

        % [1] transform from planning voxel -> RAS -> segmentation RAS
        %
        % The planning -> voxel step is computed directly (not via
        % transform_coordinates('ras_plus','grid',...)) because that helper
        % rounds to an integer voxel index, while this conversion must keep
        % the continuous (sub-voxel) position — the planning and
        % segmentation images share the same voxel grid but have different
        % headers, so reinterpreting the *unrounded* voxel position under
        % the segmentation affine is what carries the position across.
        % The voxel -> segmentation-RAS step (forward, no rounding involved)
        % reuses transform_coordinates per doc_coordinate_systems.md's rule
        % against re-implementing affine application elsewhere.

        % inv(plan_affine) gives the 0-based voxel; 'grid' is 1-based (+1)
        vox_plan = [target_ras(i,1:3) 1] * inv(plan_affine);            % Planning voxel (not rounded)
        target_ras_seg(i,:) = transform_coordinates(parameters, vox_plan(1:3) + 1, 'grid', 'ras_plus', t1seg_info);

        vox_plan = [trans_ras(i,1:3) 1] * inv(plan_affine);
        trans_ras_seg(i,:) = transform_coordinates(parameters, vox_plan(1:3) + 1, 'grid', 'ras_plus', t1seg_info);

        % [2] transform from segmentation RAS to MNI RAS by sampling the
        % SimNIBS coordinate field directly (see NEURONAV_APPLY_DEFORMATION).
        targ_mni_mm  = neuronav_apply_deformation(target_ras_seg(i,1:3), warp_field, warp_affine, nearest_idx);
        trans_mni_mm = neuronav_apply_deformation(trans_ras_seg(i,1:3), warp_field, warp_affine, nearest_idx);
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
        trans_mni_ras(i,:) = grid_to_ras(trans_mni_pos(i,:), mni_info);
    end
end