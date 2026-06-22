function [parameters, kwave_medium, medium_masks, segmentation, trans_pos, focus_pos, ...
          sensor_data, planimg, highlighted_pos] = ...
    prepare_thermal_domain(parameters, kwave_medium, medium_masks, segmentation, ...
                           trans_pos, focus_pos, has_acoustic_fov, thermal_use_acoustic_fov, ...
                           fov_full_state, fov_offset_ac, grid_dims_full, ...
                           sensor_data, planimg, highlighted_pos)
% PREPARE_THERMAL_DOMAIN  Restore full-grid state or adjust planimg for thermal FOV.
%
% When an acoustic FOV crop was active, this function decides how to hand off
% to thermal simulation:
%
%   Default (thermal_use_acoustic_fov == false):
%     Back-projects sensor_data pressure fields onto the full acoustic grid,
%     then restores kwave_medium, medium_masks, segmentation, trans_pos,
%     focus_pos, and parameters.grid.dims from the pre-crop snapshot.
%
%   Thermal-in-FOV (thermal_use_acoustic_fov == true):
%     Adjusts planimg.transf and planimg.inv_transf so that thermal_grid_setup
%     and nifti_thermal produce correctly georeferenced NIfTI outputs without
%     expanding to the full grid.
%
% When has_acoustic_fov is false the function is a no-op.
%
% Use as:
%   [parameters, kwave_medium, medium_masks, segmentation, trans_pos, focus_pos, ...
%    sensor_data, planimg, highlighted_pos] = ...
%       prepare_thermal_domain(parameters, kwave_medium, medium_masks, segmentation, ...
%                              trans_pos, focus_pos, has_acoustic_fov, thermal_use_acoustic_fov, ...
%                              fov_full_state, fov_offset_ac, grid_dims_full, ...
%                              sensor_data, planimg, highlighted_pos)

if ~has_acoustic_fov
    return
end

if ~thermal_use_acoustic_fov
    % Back-project FOV-sized sensor_data fields to full acoustic grid.
    p_fields = {'p_max_all', 'p_max_async', 'p_min_async', 'p_min_all'};
    for fi = 1:numel(p_fields)
        f = p_fields{fi};
        if ~isempty(sensor_data) && isfield(sensor_data, f)
            sensor_data.(f) = fov_backproject_volume( ...
                single(sensor_data.(f)), fov_offset_ac, grid_dims_full);
        end
    end

    % Restore full-grid medium and geometry from pre-crop snapshot.
    kwave_medium  = fov_full_state.kwave_medium;
    medium_masks  = fov_full_state.medium_masks;
    segmentation  = fov_full_state.segmentation;
    trans_pos     = fov_full_state.trans_pos;
    focus_pos     = fov_full_state.focus_pos;
    parameters.grid.dims = fov_full_state.grid_dims;
    parameters.transducer(1).trans_pos = fov_full_state.trans_pos;
    if isfield(parameters.transducer(1), 'focus_pos')
        parameters.transducer(1).focus_pos = fov_full_state.focus_pos;
    end
    if ~isempty(highlighted_pos)
        highlighted_pos = highlighted_pos + fov_offset_ac - 1;
    end

else
    % Stay in acoustic FOV coords: adjust planimg affines for thermal NIfTI export.
    %   planimg.transf:    T1w vox → full-ac vox  (row-vector convention)
    %   We want:           T1w vox → fov-ac vox
    %   fov_vox = ac_vox - (fov_offset_ac - 1)
    M_ac2fov        = eye(4);
    M_ac2fov(4,1:3) = -(fov_offset_ac - 1);
    planimg.transf  = planimg.transf * M_ac2fov;
    if isfield(planimg, 'inv_transf')
        %   inv_transf: full-ac vox → T1w vox
        %   We want:    fov-ac vox  → T1w vox
        M_fov2ac           = eye(4);
        M_fov2ac(4,1:3)    = (fov_offset_ac - 1);
        planimg.inv_transf = M_fov2ac * planimg.inv_transf;
    end
    fprintf('[pipeline] Thermal FOV = acoustic FOV (offset [%s], dims [%s]).\n', ...
        num2str(fov_offset_ac), num2str(parameters.grid.dims));
end

end
