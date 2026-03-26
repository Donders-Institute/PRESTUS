function parameters = focal_distance_calculation(parameters)
% Compute expected focal distances for (multi-)transducer setup
%
% DESCRIPTION:
%   Ensures all transducers in parameters.transducer have a valid expected_focal_distance_mm field.
%   Supports legacy single-value propagation, transducer-specific values, and automatic derivation 
%   from T1-weighted MRI grid positions when needed. Essential preprocessing for acoustic simulations 
%   requiring precise geometric focus-targeting.
%
% INPUTS:
%   parameters  - (struct) Simulation configuration with:
%                   .transducer          - (1×N struct array) Transducer configurations where each element may contain:
%                     .expected_focal_distance_mm  - (scalar) [mm] Expected geometric focal distance (optional)
%                     .trans_pos        - (1×3 vector) Transducer position [voxels] in T1 grid (optional)
%                     .focus_pos  - (1×3 vector) Target focus position [voxels] in T1 grid (optional)
%                   .expected_focal_distance_mm  - (scalar) [mm] Legacy global focal distance (optional)
%                   .data_path           - (char) Path to subject data
%                   .t1_path_template    - (char) sprintf template for T1w filename e.g. 'sub-%03d_T1w.nii.gz'
%
% OUTPUTS:
%   parameters  - (struct) Updated with expected_focal_distance_mm populated for all transducers
%
% ALGORITHM:
%   1. PROPAGATE LEGACY: If global .expected_focal_distance_mm exists, copy to all missing transducer fields
%   2. GEOMETRIC FALLBACK: For remaining unset transducers, compute Euclidean distance between .trans_pos 
%      and .focus_pos in T1 space, scaled by voxel dimensions from NIfTI header
%   3. ERROR CHECKS: Validates T1 existence and required grid positions; fails gracefully with diagnostics
%
% NOTES:
%   - Requires subject_id in scope for T1 loading (closure variable)
%   - T1 grid computation: ||focus_pos - trans_pos|| × t1_grid_step_mm
%   - Multi-transducer aware: handles heterogeneous focal distances across array
%   - No modification to existing valid values (preserves transducer-specific overrides)
%
% EXAMPLE:
%   parameters.transducer(1).expected_focal_distance_mm = 70;  % Set explicitly for primary
%   parameters.transducer(2).trans_pos = [100 120 80];     % Will auto-compute for secondary
%   parameters.transducer(2).focus_pos = [100 120 30];
%   parameters = focal_distance_calculation(parameters);
%

% 1) Implement provided expected focal distance
if isfield(parameters, 'expected_focal_distance_mm') && ~isempty(parameters.expected_focal_distance_mm)
    for ti = 1:numel(parameters.transducer)
        if ~isfield(parameters.transducer(ti), 'expected_focal_distance_mm') || ...
                isempty(parameters.transducer(ti).expected_focal_distance_mm)
            parameters.transducer(ti).expected_focal_distance_mm = parameters.expected_focal_distance_mm;
        end
    end
end

% 2) T1-based geometric derivation for unspecified distances
needs_t1 = false;
for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);
    if ~isfield(tr, 'expected_focal_distance_mm') || isempty(tr.expected_focal_distance_mm)
        needs_t1 = true;
        break
    end
end

if needs_t1
    disp('Expected focal distance not specified for all transducers, trying to get it from positions on T1 grid')

    % Load T1 header once
    filename_t1 = dir(fullfile(parameters.data_path, sprintf(parameters.t1_path_template, parameters.subject_id)));
    if isempty(filename_t1)
        error('File does not exist for T1 (t1_path_template): %s', ...
              fullfile(parameters.data_path, sprintf(parameters.t1_path_template, parameters.subject_id)));
    end
    filename_t1 = fullfile(filename_t1(1).folder, filename_t1(1).name); % first match
    t1_info = niftiinfo(filename_t1);
    t1_grid_step_mm = t1_info.PixelDimensions(1);

    % Fill missing expected_focal_distance_mm per transducer
    for ti = 1:numel(parameters.transducer)
        tr = parameters.transducer(ti);
        if ~isfield(tr, 'expected_focal_distance_mm') || isempty(tr.expected_focal_distance_mm)
            if ~isfield(tr, 'trans_pos')  || isempty(tr.trans_pos) || ...
               ~isfield(tr, 'focus_pos') || isempty(tr.focus_pos)
                error('Transducer %d: trans_pos or focus_pos missing; cannot compute expected focal distance.', ti);
            end

            focal_distance_t1 = norm(tr.focus_pos - tr.trans_pos);
            parameters.transducer(ti).expected_focal_distance_mm = focal_distance_t1 * t1_grid_step_mm;
        end
    end

    clear filename_t1 t1_info t1_grid_step_mm focal_distance_t1
end
end
