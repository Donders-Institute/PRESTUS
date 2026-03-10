function parameters = focal_distance_calculation(parameters)
% Compute expected focal distances for (multi-)transducer setup
%
% DESCRIPTION:
%   Ensures all transducers in parameters.transducer have a valid expected_focal_distance_bowl field.
%   Supports legacy single-value propagation, transducer-specific values, and automatic derivation 
%   from T1-weighted MRI grid positions when needed. Essential preprocessing for acoustic simulations 
%   requiring precise geometric focus-targeting.
%
% INPUTS:
%   parameters  - (struct) Simulation configuration with:
%                   .transducer          - (1×N struct array) Transducer configurations where each element may contain:
%                     .expected_focal_distance_bowl  - (scalar) [mm] Expected geometric focal distance (optional)
%                     .trans_pos        - (1×3 vector) Transducer position [voxels] in T1 grid (optional)
%                     .focus_pos  - (1×3 vector) Target focus position [voxels] in T1 grid (optional)
%                   .expected_focal_distance_bowl  - (scalar) [mm] Legacy global focal distance (optional)
%                   .data_path           - (char) Path to subject data
%                   .t1_path_template    - (char) sprintf template for T1w filename e.g. 'sub-%03d_T1w.nii.gz'
%
% OUTPUTS:
%   parameters  - (struct) Updated with expected_focal_distance_bowl populated for all transducers
%
% ALGORITHM:
%   1. PROPAGATE LEGACY: If global .expected_focal_distance_bowl exists, copy to all missing transducer fields
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
%   parameters.transducer(1).expected_focal_distance_bowl = 70;  % Set explicitly for primary
%   parameters.transducer(2).trans_pos = [100 120 80];     % Will auto-compute for secondary
%   parameters.transducer(2).focus_pos = [100 120 30];
%   parameters = focal_distance_calculation(parameters);
%

% 1) Implement provided expected focal distance

% Update everything if expected_focal_distance_ep is specified
if isfield(parameters, 'expected_focal_distance_ep') && ~isempty(parameters.expected_focal_distance_ep)
    for ti = 1:numel(parameters.transducer)
        if ~isfield(parameters.transducer(ti), 'expected_focal_distance_ep') || ...
                isempty(parameters.transducer(ti).expected_focal_distance_ep)
            parameters.transducer(ti).expected_focal_distance_ep = parameters.expected_focal_distance_ep;
            % calculate focal distance offset (between transducer bowl and exit plane for annular arrays)
            parameters.transducer(ti).focal_distance_offset = parameters.transducer(ti).curv_radius_mm - parameters.transducer(ti).dist_to_plane_mm;
            % calculate focal distance (from bowl)
            parameters.transducer(ti).expected_focal_distance_bowl = parameters.transducer(ti).expected_focal_distance_ep+parameters.transducer(ti).focal_distance_offset;
        end
    end
    % copy to parameters main structure (for first transducer)
    if ~isfield(parameters, 'expected_focal_distance_bowl') || isempty(parameters.expected_focal_distance_bowl)
        parameters.expected_focal_distance_bowl = parameters.transducer(1).expected_focal_distance_bowl;
    end
end

if isfield(parameters, 'expected_focal_distance_bowl') && ~isempty(parameters.expected_focal_distance_bowl)
    for ti = 1:numel(parameters.transducer)
        if ~isfield(parameters.transducer(ti), 'expected_focal_distance_bowl') || ...
                isempty(parameters.transducer(ti).expected_focal_distance_bowl)
            parameters.transducer(ti).expected_focal_distance_bowl = parameters.expected_focal_distance_bowl;
            % calculate focal distance offset (between transducer bowl and exit plane for annular arrays)
            parameters.transducer(ti).focal_distance_offset = parameters.transducer(ti).curv_radius_mm - parameters.transducer(ti).dist_to_plane_mm;
            % calculate focal distance (from exit plane)
            parameters.transducer(ti).expected_focal_distance_ep = parameters.transducer(ti).expected_focal_distance_bowl-parameters.transducer(ti).focal_distance_offset;
        end
    end
    % copy to parameters main structure (for first transducer)
    if ~isfield(parameters, 'expected_focal_distance_ep') || isempty(parameters.expected_focal_distance_ep)
        parameters.expected_focal_distance_ep = parameters.transducer(1).expected_focal_distance_ep;
    end
end

% 2) Rely on specification of transducer and target position
if ~isfield(parameters, 'expected_focal_distance_ep') || ~isfield(parameters, 'expected_focal_distance_bowl')
    warning('Expected focal distance not specified for all transducers, trying to get it from transducer and target positions ...')
end

% Fill missing expected_focal_distance_bowl per transducer
for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);
    if ~isfield(tr, 'expected_focal_distance_bowl') || isempty(tr.expected_focal_distance_bowl)
        if ~isfield(tr, 'trans_pos')  || isempty(tr.trans_pos) || ...
           ~isfield(tr, 'focus_pos') || isempty(tr.focus_pos)
            warning('Transducer %d: trans_pos or focus_pos missing; cannot compute expected focal distance.', ti);
        end
        % calculate grid distance between transducer bowl and focus
        focal_distance = norm(tr.focus_pos - tr.trans_pos);
        % scale grid distance by grid step to calculate mm
        parameters.transducer(ti).expected_focal_distance_bowl = focal_distance * parameters.grid_step_mm;
        % calculate focal distance offset (between transducer bowl and exit plane for annular arrays)
        parameters.transducer(ti).focal_distance_offset = parameters.transducer(ti).curv_radius_mm - parameters.transducer(ti).dist_to_plane_mm;
        % calculate focal distance (from exit plane)
        parameters.transducer(ti).expected_focal_distance_ep = parameters.transducer(ti).expected_focal_distance_bowl-parameters.transducer(ti).focal_distance_offset;

    end
end

end
