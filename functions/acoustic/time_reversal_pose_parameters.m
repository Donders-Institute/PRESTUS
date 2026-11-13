function [parameters, info] = time_reversal_pose_parameters(parameters, pose, options)
% TIME_REVERSAL_POSE_PARAMETERS  Run parameters that place the transducer at a pose of a pose search
%
% Turns one row of a pose table (<prefix>_timereversal_poses_T1.csv, written
% by TIME_REVERSAL_STAGE with time_reversal.placement_search) into the
% parameters of a forward run at that pose: a manual placement in RAS mm with
% the bowl apex at the pose's apex and the beam aimed at the target
% ('geometric') or at the aberration-corrected aim point ('corrected').
%
% The focal depth follows the pose (focal_distance_ep = apex-to-aim distance
% minus the bowl-to-exit-plane offset). The ring phases follow it too:
%   - drive resolved from the calibration library (transducer.serial): the
%     library phases (and amplitudes, if they came from the library) are
%     cleared, so the run resolves them again for the new focal depth;
%   - inline phases: replaced by geometric steering to the aim point in
%     water (ANNULAR_GEOMETRIC_PHASES), the drive the pose scores assume;
%   - single element: no phases to set.
% Time reversal and the pose workflow are switched off for the run.
%
% Use as:
%   [p, info] = time_reversal_pose_parameters(parameters, pose)
%   [p, info] = time_reversal_pose_parameters(parameters, pose, 'Aim', 'geometric', 'Affix', '_pose012_geo')
%
% Input:
%   parameters - PRESTUS config of the planned run (one curved annular transducer)
%   pose       - one row of a pose table (table or struct) with apex_ras_*,
%                aim_ras_*, target_ras_* [RAS mm]
%   options    - name-value pairs:
%     'Aim'   'corrected' (default) | 'geometric'
%     'Affix' appended to io.output_affix (default '')
%
% Output:
%   parameters - config of the run at that pose
%   info       - struct: apex_ras, focus_ras, target_ras, aim, focal_distance_ep,
%                phase_source ('library' | 'geometric' | 'none')
%
% See also: TIME_REVERSAL_STAGE, POSE_WORKFLOW_PIPELINE, ANNULAR_GEOMETRIC_PHASES

arguments
    parameters (1,1) struct
    pose
    options.Aim (1,:) char {mustBeMember(options.Aim, {'corrected', 'geometric'})} = 'corrected'
    options.Affix (1,:) char = ''
end

    if istable(pose)
        if height(pose) ~= 1
            error('PRESTUS:posePlacement:oneRow', 'Pass one row of the pose table.');
        end
        pose = table2struct(pose);
    end
    need = {'apex_ras_x', 'apex_ras_y', 'apex_ras_z', 'aim_ras_x', 'aim_ras_y', 'aim_ras_z', ...
            'target_ras_x', 'target_ras_y', 'target_ras_z'};
    missing = need(~isfield(pose, need));
    if ~isempty(missing)
        error('PRESTUS:posePlacement:noRas', ...
            ['The pose has no RAS positions (%s). They are written when the search run has the ' ...
             'T1 mapping (layered medium).'], strjoin(missing, ', '));
    end
    if numel(parameters.transducer) ~= 1
        error('PRESTUS:posePlacement:multiTransducer', ...
            'Pose placement supports configs with one transducer (found %d).', numel(parameters.transducer));
    end
    tr = parameters.transducer(1);
    if ~strcmp(char(tr.type), 'annular') || ~isfinite(tr.annular.curv_radius_mm)
        error('PRESTUS:posePlacement:transducerType', ...
            'Pose placement supports curved annular transducers (the pose search does too).');
    end

    apex   = [pose.apex_ras_x, pose.apex_ras_y, pose.apex_ras_z];
    target = [pose.target_ras_x, pose.target_ras_y, pose.target_ras_z];
    if strcmp(options.Aim, 'corrected')
        focus = [pose.aim_ras_x, pose.aim_ras_y, pose.aim_ras_z];
    else
        focus = target;
    end

    tr.trans_pos = apex;
    tr.focus_pos = focus;
    offset = 0;
    if isfield(tr.annular, 'dist_geom_ep_mm') && ~isempty(tr.annular.dist_geom_ep_mm)
        offset = annular_bowl_depth_mm(tr.annular);
    end
    tr.focal_distance_ep   = norm(focus - apex) - offset;
    tr.focal_distance_bowl = [];

    n_rings = numel(tr.annular.elem_od_mm);
    lib = struct('phases', false, 'amp', false);
    if isfield(parameters, 'state') && isfield(parameters.state, 'library_drive') && ...
            ~isempty(parameters.state.library_drive) && isstruct(parameters.state.library_drive{1})
        lib = parameters.state.library_drive{1};
    end
    if n_rings == 1
        phase_source = 'none';
    elseif lib.phases
        tr.annular.elem_phase_deg = [];
        tr.annular.elem_phase_rad = [];
        if lib.amp
            tr.annular.elem_p_amp = [];
        end
        phase_source = 'library';
    else
        deg = annular_geometric_phases(tr.annular, apex, focus, tr.freq_hz, ...
            parameters.medium_properties.water.sound_speed);
        tr.annular.elem_phase_deg = deg;
        tr.annular.elem_phase_rad = deg2rad(deg);
        phase_source = 'geometric';
    end
    parameters.transducer = tr;

    parameters.placement.mode = 'manual';
    parameters.placement.manual.coordinate_space = 'ras';
    parameters.modules.run_transducer_placement = 1;   % converts the RAS positions

    parameters.time_reversal.enabled = false;
    parameters.time_reversal.only    = false;
    parameters.time_reversal.placement_search.enabled = false;
    parameters.time_reversal.pose_workflow.enabled    = false;
    if isfield(parameters, 'state')
        for f = {'pose_workflow', 'time_reversal', 'library_drive'}
            if isfield(parameters.state, f{1}), parameters.state = rmfield(parameters.state, f{1}); end
        end
    end
    parameters.io.output_affix = [char(parameters.io.output_affix) options.Affix];

    info = struct('apex_ras', apex, 'focus_ras', focus, 'target_ras', target, 'aim', options.Aim, ...
        'focal_distance_ep', tr.focal_distance_ep, 'phase_source', phase_source);
end
