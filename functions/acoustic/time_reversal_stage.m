function [parameters, source, source_labels, tr_results] = time_reversal_stage(parameters, kgrid, kwave_medium, ...
        medium_masks, sensor, src_weights, trans_pos, focus_pos, source, source_labels, planimg, fov_state)
% TIME_REVERSAL_STAGE  Refocusing phases, aberration metrics and pose search from one reverse run per target
%
% Optional pipeline stage between the free-water baseline and the acoustic
% simulation (time_reversal.enabled). Per target, one k-Wave run with a
% point source at the target (REVERSE_POINT_SOURCE_FIELD) is recorded at the
% transducer surface and, with time_reversal.placement_search.enabled, over
% the region of the candidate poses. From it, for every transducer:
%
%   refocusing (TIME_REVERSAL_REFOCUS)
%     element transfers to the target, time-reversal phases, coherence of
%     the arriving contributions, the target-pressure gain of refocusing
%     and the predicted focal shift with the current and refocused drive.
%     Reported for every transducer; what can be applied depends on the
%     drive: matrix (per element), annular rings (per ring), single
%     element (nothing; the metrics still describe the aberration).
%     time_reversal.apply_phases: true replaces the element phases by the
%     time-reversal phases for the acoustic run (not for a single element).
%     The free-water baseline has already run with the geometric phases,
%     so ISPPA targeting keeps the element drive amplitudes of the
%     geometric (calibrated) drive.
%
%   pose search (RECIPROCITY_POSE_CANDIDATES, RECIPROCITY_SCORE_POSES)
%     curved annular transducers only, 3D: poses on an orbit around the
%     target, scored by predicted target pressure after aberration-corrected
%     aim. Suggestions only: the acoustic run keeps the planned pose; verify
%     the best poses with forward runs (e.g. partial FOV, then one full run;
%     POSE_WORKFLOW_PIPELINE does this). With an acoustic FOV crop active
%     (fov_state), the reverse run of a pose search uses the full grid, so
%     that poses outside the beam FOV can be scored; the refocusing metrics
%     come from the same run.
%
% Outputs (prefix sub-XXX_<medium><affix>):
%   <dir_tabular>/<prefix>_timereversal_summary.csv    one row per transducer
%   <dir_tabular>/<prefix>_timereversal_elements.csv   one row per element
%   <dir_tabular>/<prefix>_timereversal_poses_T<i>.csv pose search, per transducer
%   <dir_cache>/<prefix>_timereversal.mat              all results (+ reverse
%                                                     fields with save_field)
%
% Use as:
%   [parameters, source, source_labels, tr_results] = time_reversal_stage(parameters, kgrid, ...
%       kwave_medium, medium_masks, sensor, src_weights, trans_pos, focus_pos, source, source_labels, ...
%       planimg, fov_state)
%
% Input:
%   parameters    - PRESTUS config with a time_reversal section
%   kgrid, kwave_medium, sensor - forward k-Wave grid, medium and sensor
%   medium_masks  - grid label volume (0 = water)
%   src_weights   - element weights from SOURCE_SENSOR_SETUP / SOURCE_CREATE
%   trans_pos, focus_pos - transducer and focus positions [grid voxels]
%   source, source_labels - forward source (replaced when phases are applied)
%   planimg       - grid <-> T1 mapping (inv_transf, t1_header, fov_offset_ac);
%                   optional: without it, poses are given in grid voxels only
%   fov_state     - from ACOUSTIC_FOV_CROP (optional): full-grid medium and
%                   positions for the reverse run of a pose search
%
% Output:
%   parameters    - with time-reversal phases when applied; state.time_reversal
%   source, source_labels - rebuilt with the applied phases, else unchanged
%   tr_results    - struct array, one entry per transducer
%
% See also: REVERSE_POINT_SOURCE_FIELD, TIME_REVERSAL_REFOCUS,
%           RECIPROCITY_POSE_CANDIDATES, RECIPROCITY_SCORE_POSES

arguments
    parameters    (1,1) struct
    kgrid         (1,1)
    kwave_medium  (1,1) struct
    medium_masks  {mustBeNumeric}
    sensor        (1,1) struct
    src_weights
    trans_pos     {mustBeNumeric}
    focus_pos     {mustBeNumeric}
    source        (1,1) struct
    source_labels {mustBeNumericOrLogical}
    planimg       (1,1) struct = struct()
    fov_state     (1,1) struct = struct('active', false)
end

    cfg = parameters.time_reversal;
    nT  = numel(parameters.transducer);
    nd  = numel(parameters.grid.dims);
    dx  = parameters.grid.resolution_mm * 1e-3;
    if isempty(src_weights)
        error('PRESTUS:timeReversal:noWeights', ...
            'Time reversal needs the element weights of the source (not available for axisymmetric grids).');
    end
    search = isfield(cfg, 'placement_search') && cfg.placement_search.enabled;
    if search && nd ~= 3
        warning('PRESTUS:timeReversal:search2D', 'The pose search needs a 3D grid; skipped.');
        search = false;
    end

    % Forward-grid positions (for rebuilding the source) and the context of
    % the reverse run: the forward grid, or the full grid when a pose search
    % meets an acoustic FOV crop
    [trans_rows, focus_rows] = resolve_rows(parameters, trans_pos, focus_pos);
    ctx = struct('parameters', parameters, 'kgrid', kgrid, 'medium', kwave_medium, ...
        'masks', medium_masks, 'sensor', sensor, 'w', src_weights, ...
        'trans', trans_rows, 'focus', focus_rows, 'planimg', planimg);
    if search && isfield(fov_state, 'active') && fov_state.active
        ctx = full_grid_context(parameters, kgrid, fov_state, planimg);
    end
    dims    = ctx.parameters.grid.dims;
    w       = ctx.w;
    offsets = [0, cumsum(w.elem_counts(1:end-1))];

    % Source point positions [m]
    sub = cell(1, nd);
    [sub{:}] = ind2sub(dims, w.mask_ind);
    pts_m = (cell2mat(sub) - 1) * dx;

    % Pose candidates per transducer (before the reverse runs: they set the sensor region)
    cand = cell(1, nT);
    if search
        ps = cfg.placement_search;
        for it = 1:nT
            tr = parameters.transducer(it);
            if strcmp(tr.type, 'annular') && isfinite(tr.annular.curv_radius_mm)
                cand{it} = reciprocity_pose_candidates(ctx.parameters, ctx.masks, tr, ...
                    ctx.trans(it, :), ctx.focus(it, :), 'MaxAngle', ps.max_angle_deg, ...
                    'AngleStep', ps.angle_step_deg, 'StandoffMm', ps.standoff_mm, ...
                    'JitterMm', ps.jitter_mm);
            else
                fprintf('[time_reversal] pose search skipped for transducer %d (curved annular only)\n', it);
            end
        end
    end

    % One reverse run per distinct target
    [targets, ~, grp] = unique(round(ctx.focus), 'rows', 'stable');
    tr_results = struct('transducer', num2cell(1:nT));
    rev_fields = cell(1, size(targets, 1));
    for g = 1:size(targets, 1)
        members = find(grp == g)';
        mask = false(dims);
        for it = members
            cols = offsets(it) + (1:w.elem_counts(it));
            mask(w.mask_ind(any(w.W(:, cols) ~= 0, 2))) = true;
            if ~isempty(cand{it}), mask = mask | cand{it}.sensor_mask; end
        end
        freq = parameters.transducer(members(1)).freq_hz;
        P = reverse_point_source_field(ctx.parameters, ctx.kgrid, ctx.medium, ctx.sensor, targets(g, :), mask, freq);
        P_vol = complex(nan(dims, 'single'));
        P_vol(mask) = P;
        if cfg.save_field, rev_fields{g} = P_vol; end

        c_t = ctx.medium.sound_speed;
        if ~isscalar(c_t), c_t = c_t(sub2ind_nd(dims, targets(g, :))); end
        k_target = 2 * pi * freq / c_t;
        k_water  = 2 * pi * freq / parameters.medium_properties.water.sound_speed;
        target_m = (targets(g, :) - 1) * dx;

        for it = members
            tr   = parameters.transducer(it);
            cols = offsets(it) + (1:w.elem_counts(it));
            C_elem = w.drive(cols);                      % the drive the source was built with
            beam_u = target_m - (ctx.trans(it, :) - 1) * dx;
            beam_u = beam_u / norm(beam_u);

            res = time_reversal_refocus(double(P_vol(w.mask_ind)), w.W(:, cols), ...
                C_elem, pts_m, target_m, beam_u, k_target, tr.type);
            res.applied = false;
            if cfg.apply_phases && ~strcmp(res.capability, 'none')
                if strcmp(tr.type, 'matrix')
                    % CREATE_MATRIX_KARRAY recomputes geometric phases unless overridden
                    parameters.transducer(it).matrix.elem_phase_override_deg = rad2deg(res.phase_tr_rad(:)');
                else
                    parameters.transducer(it).(tr.type).elem_phase_rad = res.phase_tr_rad(:)';
                    parameters.transducer(it).(tr.type).elem_phase_deg = rad2deg(res.phase_tr_rad(:)');
                end
                res.applied = true;
            end
            fprintf(['[time_reversal] transducer %d (%s, %d elements): coherence %.2f (elements) / %.2f (aperture), ' ...
                'refocus gain x%.2f pressure, predicted focal shift %.1f mm (axial %+.1f, lateral %.1f)%s\n'], ...
                it, res.capability, numel(cols), res.coherence_elements, res.coherence_aperture, res.refocus_gain, ...
                res.pred_current.shift_mm, res.pred_current.shift_axial_mm, res.pred_current.shift_lateral_mm, ...
                ternary(res.applied, ' — phases applied', ''));
            tr_results(it).refocus = res;

            if ~isempty(cand{it})
                T = reciprocity_score_poses(cand{it}, P_vol, tr, k_water, k_target, ...
                    'AimIterations', cfg.placement_search.aim_iterations);
                T = add_t1_positions(T, ctx.planimg, targets(g, :));
                tr_results(it).poses = T;
                best = T(1, :);
                fprintf(['[time_reversal] transducer %d pose search: best pose %.1f deg from planned, ' ...
                    'x%.2f target pressure (corrected aim, %.1f mm off target)\n'], ...
                    it, best.polar_deg, best.p_corrected_rel, best.aim_offset_mm);
            end
        end
    end

    write_outputs(parameters, tr_results, rev_fields, targets);
    parameters.state.time_reversal = struct('applied', any(arrayfun(@(r) r.refocus.applied, tr_results)));

    if parameters.state.time_reversal.applied
        disp('[time_reversal] Rebuilding the source with the time-reversal phases...')
        [source, source_labels] = source_create(parameters, kgrid, trans_rows, focus_rows);
    end
end

function ctx = full_grid_context(parameters, kgrid, fov_state, planimg)
% Grid, medium, sensor and element weights of the full (uncropped) grid, with
% the forward run's time axis, for a reverse run that covers poses outside the
% acoustic FOV. The source is built under its own affix and never cached.
    disp('[time_reversal] acoustic FOV active: the pose search runs its reverse simulation on the full grid');
    p = parameters;
    p.grid.dims = fov_state.full_dims;
    p = fov_restore_transducers(p, fov_state);
    p.io.output_affix         = [char(parameters.io.output_affix) '_trfull'];
    p.io.overwrite_files      = 'always';
    p.io.save_source_matrices = 0;
    c = fov_state.kwave_medium_full.sound_speed;
    [kg, ~, sens, ~, w] = source_sensor_setup(p, max(c(:)), fov_state.trans_pos_full, ...
        fov_state.focus_pos_full, kgrid.dt, min(c(:)));
    [trans_rows, focus_rows] = resolve_rows(p, fov_state.trans_pos_full, fov_state.focus_pos_full);
    planimg.fov_offset_ac = [1 1 1];
    ctx = struct('parameters', p, 'kgrid', kg, 'medium', fov_state.kwave_medium_full, ...
        'masks', fov_state.medium_masks_full, 'sensor', sens, 'w', w, ...
        'trans', trans_rows, 'focus', focus_rows, 'planimg', planimg);
end

function [trans_rows, focus_rows] = resolve_rows(parameters, trans_pos, focus_pos)
% Per-transducer positions as SOURCE_SENSOR_SETUP resolves them
    nT = numel(parameters.transducer);
    nd = numel(parameters.grid.dims);
    trans_rows = repmat(reshape(trans_pos(1, :), 1, nd), nT, 1);
    focus_rows = repmat(reshape(focus_pos(1, :), 1, nd), nT, 1);
    if nT > 1
        for ti = 1:nT
            tr = parameters.transducer(ti);
            if isfield(tr, 'trans_pos') && ~isempty(tr.trans_pos), trans_rows(ti, :) = tr.trans_pos(:)'; end
            if isfield(tr, 'focus_pos') && ~isempty(tr.focus_pos), focus_rows(ti, :) = tr.focus_pos(:)'; end
        end
    end
end

function T = add_t1_positions(T, planimg, target_vox)
% Apex, corrected aim point and target in T1 voxels and RAS mm, so that a pose
% can be set up as a placement (trans_pos = apex, focus_pos = aim point)
    if ~isfield(planimg, 'inv_transf') || isempty(planimg.inv_transf)
        return
    end
    off = [1 1 1];
    if isfield(planimg, 'fov_offset_ac') && ~isempty(planimg.fov_offset_ac), off = planimg.fov_offset_ac; end
    to_t1 = @(v) [v + off - 1, ones(size(v, 1), 1)] * planimg.inv_transf;
    ax = 'xyz';
    for nm = {'apex', 'aim'}
        g = [T.([nm{1} '_vox_x']), T.([nm{1} '_vox_y']), T.([nm{1} '_vox_z'])];
        t1 = to_t1(g); t1 = t1(:, 1:3);
        ras = grid_to_ras(t1, planimg.t1_header);
        for d = 1:3
            T.(sprintf('%s_t1_%s', nm{1}, ax(d)))  = t1(:, d);
            T.(sprintf('%s_ras_%s', nm{1}, ax(d))) = ras(:, d);
        end
    end
    t1 = to_t1(target_vox); ras = grid_to_ras(t1(1:3), planimg.t1_header);
    for d = 1:3
        T.(sprintf('target_ras_%s', ax(d))) = repmat(ras(d), height(T), 1);
    end
end

function lin = sub2ind_nd(dims, v)
    v = num2cell(round(v));
    lin = sub2ind(dims, v{:});
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function write_outputs(parameters, tr_results, rev_fields, targets)
    prefix = sprintf('sub-%03d_%s%s', parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix);
    summ = {}; elems = {};
    for it = 1:numel(tr_results)
        r = tr_results(it).refocus;
        pc = r.pred_current; pr = r.pred_refocused;
        if isempty(pr), pr = struct('shift_mm', NaN, 'shift_axial_mm', NaN, 'shift_lateral_mm', NaN); end
        summ{end+1} = struct('transducer', it, 'capability', r.capability, ...
            'coherence_elements', r.coherence_elements, 'coherence_aperture', r.coherence_aperture, ...
            'refocus_gain_pressure', r.refocus_gain, 'refocus_gain_intensity', r.refocus_gain^2, ...
            'shift_mm', pc.shift_mm, 'shift_axial_mm', pc.shift_axial_mm, 'shift_lateral_mm', pc.shift_lateral_mm, ...
            'shift_refocused_mm', pr.shift_mm, 'shift_refocused_axial_mm', pr.shift_axial_mm, ...
            'shift_refocused_lateral_mm', pr.shift_lateral_mm, 'phases_applied', r.applied); %#ok<AGROW>
        n_el = numel(r.G);
        elems{end+1} = table(repmat(it, n_el, 1), (1:n_el)', abs(r.G) / max(abs(r.G)), ...
            rad2deg(r.phase_current_rad(:)), rad2deg(r.phase_tr_rad(:)), r.phase_error_deg(:), ...
            'VariableNames', {'transducer', 'element', 'transfer_rel', 'phase_current_deg', ...
                              'phase_timereversal_deg', 'arrival_phase_error_deg'}); %#ok<AGROW>
        if isfield(tr_results, 'poses') && ~isempty(tr_results(it).poses)
            writetable(tr_results(it).poses, fullfile(parameters.io.dir_tabular, ...
                sprintf('%s_timereversal_poses_T%d.csv', prefix, it)));
        end
    end
    writetable(struct2table([summ{:}]), fullfile(parameters.io.dir_tabular, [prefix '_timereversal_summary.csv']));
    writetable(vertcat(elems{:}), fullfile(parameters.io.dir_tabular, [prefix '_timereversal_elements.csv']));
    mat_file = fullfile(parameters.io.dir_cache, [prefix '_timereversal.mat']);
    if any(~cellfun(@isempty, rev_fields))
        save(mat_file, 'tr_results', 'rev_fields', 'targets', '-v7.3');
    else
        save(mat_file, 'tr_results', 'targets', '-v7.3');
    end
    fprintf('[time_reversal] results: %s\n', fullfile(parameters.io.dir_tabular, [prefix '_timereversal_summary.csv']));
end
