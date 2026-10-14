function jobs = pose_workflow_pipeline(parameters, options)
% POSE_WORKFLOW_PIPELINE  Transducer pose search, forward verification and final run in one call
%
% time_reversal.pose_workflow.enabled: true turns a run into a placement
% workflow for its (single, curved annular) transducer:
%
%   1. search   one run with time_reversal.placement_search and
%               time_reversal.only: a reverse simulation from the target on
%               the full grid scores poses around the planned one (no
%               forward simulation).
%               -> <prefix><affix>_posesearch_timereversal_poses_T1.csv
%   2. verify   forward runs of the planned pose (reference) and the n_top
%               best poses (plus pose_workflow.poses), each with geometric
%               and/or corrected aim (pose_workflow.aim), on a partial FOV
%               by default (pose_workflow.verify_fov); acoustics only.
%               -> <prefix><affix>_poseverify_selection.csv
%   3. final    target intensity and focus position of every verification
%               run (TIME_REVERSAL_POSE_METRICS), ranked (pose_workflow.rank_by);
%               then one full-grid run of the best pose with the config's
%               modules (thermal etc.) when pose_workflow.final_run.
%               -> <prefix><affix>_poseverify_results.csv, run <affix>_posebest
%
% On MATLAB the stages run one after the other. On SLURM/qsub each stage is
% submitted; a small CPU job after the search jobs submits the verification
% runs, and one after those collects the results and submits the final run
% (compute nodes must be allowed to submit jobs). A stage whose outputs exist
% is not run again unless io.overwrite_files is 'always', so a call after a
% partial failure continues where it stopped.
%
% Use as:
%   jobs = pose_workflow_pipeline(parameters)
%   (normally via PRESTUS_PIPELINE_START or PRESTUS_PIPELINE)
%
% Input:
%   parameters - PRESTUS config with time_reversal.pose_workflow.enabled;
%                parameters.state.pose_workflow (set by this function) selects
%                the stage of a submitted follow-up job
%   options    - pipeline options (unused; for the PRESTUS_PIPELINE interface)
%
% Output:
%   jobs - struct of the scheduler job IDs submitted by this call (HPC)
%
% See also: TIME_REVERSAL_STAGE, TIME_REVERSAL_POSE_PARAMETERS,
%           TIME_REVERSAL_POSE_METRICS, PRESTUS_PIPELINE_START

arguments
    parameters (1,1) struct
    options struct = struct() %#ok<INUSA>
end

    jobs = struct();
    if ~isfield(parameters, 'platform') || strcmp(parameters.platform, 'auto')
        parameters.platform = hpc_detect_system();
    end
    on_hpc = ~strcmp(parameters.platform, 'matlab');

    stage = 'search';
    if isfield(parameters, 'state') && isfield(parameters.state, 'pose_workflow')
        wf = parameters.state.pose_workflow;
        stage = wf.stage;
        parameters.simulation.code_type = wf.code_type;   % undo the dispatcher job's CPU settings
        parameters.hpc                  = wf.hpc;
        parameters.io.output_affix      = wf.base_affix;
        parameters.state = rmfield(parameters.state, 'pose_workflow');
    else
        validate_workflow(parameters);
    end
    cfg = parameters.time_reversal.pose_workflow;
    f   = workflow_files(parameters);

    if strcmp(stage, 'search')
        fprintf('[pose workflow] 1/3 pose search\n');
        if isfile(f.poses) && ~overwrite_all(parameters)
            fprintf('[pose workflow] pose table exists, search skipped: %s\n', f.poses);
            ids = [];
        else
            ids = run_or_submit(search_parameters(parameters, cfg), on_hpc, true);
        end
        jobs.search = ids;
        if on_hpc
            jobs.verify_dispatch = dispatch(parameters, 'verify', ids, 'afterok');
            return
        end
        stage = 'verify';
    end

    if strcmp(stage, 'verify')
        fprintf('[pose workflow] 2/3 forward verification\n');
        if ~isfile(f.poses)
            error('PRESTUS:poseWorkflow:noPoses', 'No pose table (did the search run fail?): %s', f.poses);
        end
        sel = select_poses(readtable(f.poses, 'TextType', 'string'), cfg, parameters.io.output_affix);
        writetable(sel, f.selection);
        fprintf('[pose workflow] %d verification runs: %s\n', height(sel), f.selection);
        ids = [];
        for r = 1:height(sel)
            if isfile(targeting_csv(parameters, sel.affix(r))) && ~overwrite_all(parameters)
                fprintf('[pose workflow] %s exists, skipped\n', sel.affix(r));
                continue
            end
            p = time_reversal_pose_parameters(parameters, pose_struct(sel(r, :)), 'Aim', 'corrected');
            p.io.output_affix = char(sel.affix(r));
            p = verify_settings(p, cfg);
            ids = [ids, run_or_submit(p, on_hpc, false)]; %#ok<AGROW>
        end
        jobs.verify = ids;
        if on_hpc
            jobs.final_dispatch = dispatch(parameters, 'final', ids, 'afterany');
            return
        end
        stage = 'final';
    end

    if strcmp(stage, 'final')
        fprintf('[pose workflow] 3/3 results and final run\n');
        res = collect_results(parameters, readtable(f.selection, 'TextType', 'string'), cfg);
        writetable(res, f.results);
        fprintf('[pose workflow] results: %s\n', f.results);
        best = res(res.rank == 1, :);
        if isempty(best)
            warning('PRESTUS:poseWorkflow:noResults', 'No verification run has results; no final run.');
            return
        end
        if cfg.final_run
            affix = [char(parameters.io.output_affix) '_posebest'];
            if isfile(targeting_csv(parameters, affix)) && ~overwrite_all(parameters)
                fprintf('[pose workflow] final run %s exists, skipped\n', affix);
            else
                p = time_reversal_pose_parameters(parameters, pose_struct(best), 'Aim', 'corrected');
                p.io.output_affix = affix;
                p = set_fov(p, 'full');
                fprintf('[pose workflow] final full-grid run of pose %d (%s aim): %s\n', ...
                    best.pose, best.aim, affix);
                jobs.final = run_or_submit(p, on_hpc, false);
            end
        end
    end
end

%% ── Stages ───────────────────────────────────────────────────────────────

function validate_workflow(parameters)
    if numel(parameters.transducer) ~= 1
        error('PRESTUS:poseWorkflow:multiTransducer', ...
            'The pose workflow supports one transducer (found %d).', numel(parameters.transducer));
    end
    tr = parameters.transducer(1);
    if ~strcmp(char(tr.type), 'annular') || ~isfinite(tr.annular.curv_radius_mm)
        error('PRESTUS:poseWorkflow:transducerType', 'The pose workflow needs a curved annular transducer.');
    end
    if ~contains(char(parameters.simulation.medium), 'layered')
        error('PRESTUS:poseWorkflow:medium', ...
            'The pose workflow needs a head model (simulation.medium: layered) for RAS positions.');
    end
    if is_multi_isppa_mode(parameters) || ...
            (isfield(parameters.simulation, 'uncertainty') && parameters.simulation.uncertainty)
        error('PRESTUS:poseWorkflow:mode', ...
            'The pose workflow cannot be combined with multi-ISPPA or uncertainty runs.');
    end
end

function p = search_parameters(parameters, cfg)
    p = parameters;
    p.io.output_affix = [char(parameters.io.output_affix) '_posesearch'];
    p.time_reversal.enabled      = true;
    p.time_reversal.only         = true;
    p.time_reversal.apply_phases = false;
    p.time_reversal.placement_search.enabled = true;
    p.time_reversal.pose_workflow.enabled    = false;
    p = set_fov(p, 'full');                 % the candidate poses lie outside a beam FOV
    p.modules.generate_report = 0;
    if ~isempty(cfg.search_memorylimit)
        p.hpc.memorylimit = cfg.search_memorylimit;
    end
    p.hpc.job_name = sprintf('PRESTUS-posesearch_sub-%03d', p.subject_id);
end

function p = verify_settings(p, cfg)
    p = set_fov(p, cfg.verify_fov);
    p.modules.run_acoustic_sims      = 1;
    p.modules.run_acoustic_analysis  = 1;
    p.modules.run_nifti_creation     = 1;
    p.modules.run_targeting_analysis = 1;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.hpc.job_name = sprintf('PRESTUS-poseverify_sub-%03d', p.subject_id);
end

function sel = select_poses(T, cfg, base_affix)
% Planned pose (reference, geometric aim) and the n_top best poses by predicted
% target pressure after aim correction, plus pose_workflow.poses
    if ismember('complete', T.Properties.VariableNames)
        T = T(T.complete == 1, :);
    end
    planned = strcmp(T.kind, 'orbit') & T.polar_deg == 0;
    if ~any(planned)
        error('PRESTUS:poseWorkflow:noPlanned', 'The pose table has no planned pose (orbit, polar 0).');
    end
    others = sortrows(T(~planned, :), 'p_corrected_rel', 'descend');
    pick   = others(1:min(cfg.n_top, height(others)), :);
    if ~isempty(cfg.poses)
        extra = T(ismember(T.pose, cfg.poses) & ~ismember(T.pose, pick.pose) & ~planned, :);
        pick  = [pick; extra];
        if any(~ismember(cfg.poses, T.pose))
            warning('PRESTUS:poseWorkflow:unknownPose', 'Poses not in the table (or incomplete): %s', ...
                num2str(setdiff(cfg.poses(:)', T.pose(:)')));
        end
    end
    switch cfg.aim
        case 'corrected', aims = {'corrected'};
        case 'geometric', aims = {'geometric'};
        case 'both',      aims = {'geometric', 'corrected'};
        otherwise
            error('PRESTUS:poseWorkflow:aim', 'pose_workflow.aim must be corrected, geometric or both.');
    end
    rows = {sel_row(T(find(planned, 1), :), 'geometric', true, base_affix)};
    if ismember('corrected', aims)
        rows{end+1} = sel_row(T(find(planned, 1), :), 'corrected', true, base_affix);
    end
    for i = 1:height(pick)
        for a = 1:numel(aims)
            rows{end+1} = sel_row(pick(i, :), aims{a}, false, base_affix); %#ok<AGROW>
        end
    end
    sel = vertcat(rows{:});
end

function r = sel_row(t, aim, is_ref, base_affix)
    tag = 'geo';
    if strcmp(aim, 'corrected'), tag = 'aim'; end
    if strcmp(aim, 'corrected')
        focus  = [t.aim_ras_x, t.aim_ras_y, t.aim_ras_z];
        p_rel  = t.p_corrected_rel;
        shift  = t.residual_shift_mm;
    else
        focus  = [t.target_ras_x, t.target_ras_y, t.target_ras_z];
        p_rel  = t.p_geometric_rel;
        shift  = hypot(t.shift_axial_mm, t.shift_lateral_mm);
    end
    r = table(t.pose, string(t.kind), t.polar_deg, t.azimuth_deg, string(aim), is_ref, ...
        string(sprintf('%s_pose%03d_%s', base_affix, t.pose, tag)), ...
        t.apex_ras_x, t.apex_ras_y, t.apex_ras_z, focus(1), focus(2), focus(3), ...
        t.target_ras_x, t.target_ras_y, t.target_ras_z, p_rel, p_rel^2, shift, ...
        'VariableNames', {'pose', 'kind', 'polar_deg', 'azimuth_deg', 'aim', 'reference', 'affix', ...
        'apex_ras_x', 'apex_ras_y', 'apex_ras_z', 'focus_ras_x', 'focus_ras_y', 'focus_ras_z', ...
        'target_ras_x', 'target_ras_y', 'target_ras_z', 'predicted_p_rel', ...
        'predicted_intensity_rel', 'predicted_shift_mm'});
end

function res = collect_results(parameters, sel, cfg)
    n = height(sel);
    M = cell(n, 1);
    for r = 1:n
        apex   = [sel.apex_ras_x(r), sel.apex_ras_y(r), sel.apex_ras_z(r)];
        target = [sel.target_ras_x(r), sel.target_ras_y(r), sel.target_ras_z(r)];
        M{r} = time_reversal_pose_metrics(parameters, char(sel.affix(r)), target, apex, cfg.target_radius_mm);
    end
    res = [sel, struct2table([M{:}])];
    if ~any(res.found)
        res.intensity_rel = nan(n, 1);
        res.rank = nan(n, 1);
        return
    end
    ref = find(res.reference & res.aim == "geometric" & res.found, 1);
    if isempty(ref)
        warning('PRESTUS:poseWorkflow:noReference', ...
            'The planned-pose run has no results; intensity_rel (relative to it) is NaN.');
        res.intensity_rel = nan(n, 1);
    else
        res.intensity_rel = res.intensity_target / res.intensity_target(ref);
    end

    ok = res.found & isfinite(res.intensity_target);
    if ~isempty(cfg.max_focus_error_mm)
        ok = ok & res.centroid_error_mm <= cfg.max_focus_error_mm;
    end
    switch cfg.rank_by
        case 'target_intensity', key = -res.intensity_target;
        case 'focus_error',      key = res.centroid_error_mm;
        otherwise
            error('PRESTUS:poseWorkflow:rankBy', 'pose_workflow.rank_by must be target_intensity or focus_error.');
    end
    key(~ok) = Inf;
    [~, order] = sort(key);
    res.rank = nan(n, 1);
    res.rank(order(1:nnz(ok))) = 1:nnz(ok);

    fprintf('[pose workflow] %d of %d runs with results; ranked by %s\n', nnz(res.found), n, cfg.rank_by);
    fprintf('  %4s %5s %-9s %7s %9s %9s %9s\n', 'rank', 'pose', 'aim', 'polar', 'I_pred', 'I_sim', 'err_mm');
    for i = order(1:min(nnz(ok), 8))'
        fprintf('  %4d %5d %-9s %7.1f %9.2f %9.2f %9.1f%s\n', res.rank(i), res.pose(i), res.aim(i), ...
            res.polar_deg(i), res.predicted_intensity_rel(i), res.intensity_rel(i), res.centroid_error_mm(i), ...
            ternary(res.reference(i), '  (planned)', ''));
    end
    found = res.found & isfinite(res.intensity_rel);
    if nnz(found) >= 3
        fprintf('[pose workflow] rank agreement predicted vs simulated intensity: Spearman %.2f (n = %d)\n', ...
            spearman(res.predicted_intensity_rel(found), res.intensity_rel(found)), nnz(found));
    end
end

%% ── Helpers ──────────────────────────────────────────────────────────────

function f = workflow_files(parameters)
    out = get_output_dir(parameters);
    prefix = sprintf('sub-%03d_%s%s', parameters.subject_id, char(parameters.simulation.medium), ...
        char(parameters.io.output_affix));
    f.poses     = fullfile(out, [prefix '_posesearch_timereversal_poses_T1.csv']);
    f.selection = fullfile(out, [prefix '_poseverify_selection.csv']);
    f.results   = fullfile(out, [prefix '_poseverify_results.csv']);
end

function file = targeting_csv(parameters, affix)
    file = fullfile(get_output_dir(parameters), sprintf('sub-%03d_%s%s_targeting.csv', ...
        parameters.subject_id, char(parameters.simulation.medium), char(affix)));
end

function tf = overwrite_all(parameters)
    tf = isfield(parameters.io, 'overwrite_files') && strcmp(char(parameters.io.overwrite_files), 'always');
end

function p = set_fov(p, mode)
    keys = {'acoustic_fov_length_mm', 'acoustic_fov_diameter_mm'};
    switch mode
        case 'partial'
            for k = keys
                v = [];
                if isfield(p.grid, k{1}), v = p.grid.(k{1}); end
                if isempty(v), p.grid.(k{1}) = 'auto'; end
            end
        case 'full'
            for k = keys, p.grid.(k{1}) = []; end
        case 'config'
        otherwise
            error('PRESTUS:poseWorkflow:fov', 'pose_workflow.verify_fov must be partial, full or config.');
    end
end

function s = pose_struct(row)
% Pose of a selection/results row as TIME_REVERSAL_POSE_PARAMETERS reads it,
% with its focus (target or corrected aim) as the aim point
    s = struct('apex_ras_x', row.apex_ras_x, 'apex_ras_y', row.apex_ras_y, 'apex_ras_z', row.apex_ras_z, ...
        'aim_ras_x', row.focus_ras_x, 'aim_ras_y', row.focus_ras_y, 'aim_ras_z', row.focus_ras_z, ...
        'target_ras_x', row.target_ras_x, 'target_ras_y', row.target_ras_y, 'target_ras_z', row.target_ras_z);
end

function ids = run_or_submit(p, on_hpc, required)
% Run inline (MATLAB) or submit; returns the numeric job ID ([] inline)
    ids = [];
    p.time_reversal.pose_workflow.enabled = false;
    if on_hpc
        p.hpc.wait_for_job = false;
        jid = prestus_pipeline_start(p);
        ids = str2double(regexp(char(string(jid)), '^\d+', 'match', 'once'));
        if ~isfinite(ids), ids = []; end
        return
    end
    try
        prestus_pipeline(p);
    catch ME
        if required, rethrow(ME); end
        warning('PRESTUS:poseWorkflow:runFailed', 'Run %s failed: %s', p.io.output_affix, ME.message);
    end
end

function job_id = dispatch(parameters, stage, dep_ids, dep_type)
% Small CPU job that runs the next workflow stage once dep_ids have ended
    p = parameters;
    p.state.pose_workflow = struct('stage', stage, 'base_affix', char(parameters.io.output_affix), ...
        'code_type', parameters.simulation.code_type, 'hpc', parameters.hpc);
    p.io.output_affix = [char(parameters.io.output_affix) '_poseworkflow'];
    p.simulation.code_type = 'matlab_cpu';
    for fn = {'gpu', 'partition', 'depend_job_id', 'depend_job_ids'}
        if isfield(p.hpc, fn{1}), p.hpc = rmfield(p.hpc, fn{1}); end
    end
    p.hpc.timelimit    = '01:00:00';
    p.hpc.memorylimit  = 16;
    p.hpc.wait_for_job = false;
    p.hpc.job_name     = sprintf('PRESTUS-pose%s_sub-%03d', stage, p.subject_id);
    if ~isempty(dep_ids)
        p.hpc.depend_job_ids = dep_ids;
        p.hpc.depend_type    = dep_type;
    end
    job_id = prestus_pipeline_start(p);
    fprintf('[pose workflow] %s stage: job %s after %d job(s)\n', stage, num2str(job_id), numel(dep_ids));
end

function rho = spearman(a, b)
    rho = corr_pearson(tied_rank(a), tied_rank(b));
end

function r = tied_rank(x)
    [xs, i] = sort(x(:));
    r = zeros(numel(x), 1);
    k = 1;
    while k <= numel(xs)
        j = k;
        while j < numel(xs) && xs(j + 1) == xs(k), j = j + 1; end
        r(i(k:j)) = (k + j) / 2;
        k = j + 1;
    end
end

function c = corr_pearson(a, b)
    a = a - mean(a); b = b - mean(b);
    c = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
