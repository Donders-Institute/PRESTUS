function parameters = preproc_transducer_placement(parameters)
% PREPROC_TRANSDUCER_PLACEMENT  Resolve transducer and focus positions
%
% Dispatches transducer placement to one of six modes controlled by
% parameters.placement.mode:
%
%   'manual'    Use trans_pos / focus_pos already in parameters (default).
%               placement.manual.coordinate_space sets their space:
%               'voxel' (default; 1-based T1 voxel indices), 'ras' (subject
%               RAS+ mm of the T1) or 'mni' (MNI mm; transducer via the 12dof
%               transform, focus via the nonlinear warp, no scalp snapping).
%               ras/mni are converted to T1 voxels here, per transducer.
%
%   'localite'  Derive trans_pos / focus_pos from a recorded Localite
%               instrument pose (NEURONAV_LOCALITE_POSE): either
%               placement.localite.file (single XML), or the ingested
%               markers under path.localite_post (ingested from
%               path.localite_raw on first use).
%
%   'heuristic' Run the sphere-expansion heuristic search to find the optimal
%               transducer position on the skull surface for a specified MNI
%               target. Requires placement.heuristic.mni_target_mm and
%               placement.heuristic.target_name to be set. On MATLAB platform
%               the positioning runs synchronously and the resolved positions
%               are written back into parameters. On HPC platforms
%               transducer_positioning_start is submitted as a job and the
%               pipeline exits so the user can re-run after the job finishes
%               (same pattern as segmentation_only).
%
%   'mni'       Place the transducer and focus directly from MNI coordinates.
%               Requires placement.mni.trans_pos_mm and placement.mni.focus_pos_mm
%               ([1x3] MNI mm). The focus is converted with the nonlinear (nonl)
%               SimNIBS transform (like the heuristic target); the transducer
%               point is converted with the linear (12dof) transform because it
%               lies outside the brain where the nonlinear warp is unreliable.
%               The transducer is then snapped to the nearest scalp voxel and
%               offset outward so its exit plane lies placement.skin_gap_mm off
%               the scalp, as in the heuristic (no search).
%
%   'plantus'   Call the external PlanTUS tool (Lueckel et al., Mainz) to
%               optimise transducer placement against five geometric objectives:
%               beam–target overlap, angle at skin, skin–skull angle, skull
%               thickness, and distance to target. Requires
%               placement.plantus.script_path, placement.plantus.env_path,
%               placement.plantus.focal_distance_list, placement.plantus.flhm_list,
%               and placement.plantus.target_name. The target can be specified
%               as placement.plantus.mni_target_mm (MNI space) or via a
%               pre-set transducer.focus_pos. PlanTUS must be installed
%               separately (see https://github.com/mlueckel/PlanTUS).
%
%   'tus_entry' Call the external TUS_entry R toolbox (Atkinson et al.), which
%               picks the scalp entry with the most even transducer-rim to
%               scalp distance (flattest contact) within a scalp-target depth
%               range. Requires placement.tus_entry.script_path and
%               placement.tus_entry.env_path (conda env with Rscript); the
%               target is placement.tus_entry.mni_target_mm, a target mask or
%               transducer.focus_pos. The head mask defaults to the SimNIBS
%               segmentation. See POSITION_TRANSDUCER_TUS_ENTRY.
%
% The resolved positions are written into parameters.transducer(ti).trans_pos
% and parameters.transducer(ti).focus_pos for every configured transducer.
%
% heuristic, plantus and tus_entry write everything they produce to their
% placement run folder, <placement dir>/<mode>/<target_name>/
% (PLACEMENT_RUN_DIR), and finish with a run record (placement_run.json
% sidecar: resolved positions + the settings they depend on). placement.reuse
% decides whether a recorded run is used instead of running the mode again,
% independent of io.overwrite_files (which governs the within-run caches):
%   auto   (default) reuse a record whose settings match the current ones;
%          otherwise run again and replace the record. Under
%          io.overwrite_files: never a mismatch is handled like 'always'.
%   always require a record; a settings mismatch is handled per
%          io.cache_signature_check (error | warn | off).
%   never  always run the mode.
% Records written as placement_run.mat by earlier versions are still read.
%
% Use as:
%   parameters = preproc_transducer_placement(parameters)
%
% Input:
%   parameters - (1,1) simulation configuration struct
%
% Output:
%   parameters - updated struct; on HPC heuristic mode the function does not
%                return (pipeline exits after job submission)
%
% See also: POSITION_TRANSDUCER_LOCALITE, POSITION_TRANSDUCER_PLANTUS,
%           POSITION_TRANSDUCER_TUS_ENTRY, TRANSDUCER_POSITIONING_START, NEURONAV_LOCALITE_POSE,
%           PLACEMENT_RUN_DIR, GET_PLACEMENT_DIR

arguments
    parameters (1,1) struct
end

    mode = 'manual';
    if isfield(parameters, 'placement') && isfield(parameters.placement, 'mode') ...
            && ~isempty(parameters.placement.mode)
        mode = lower(char(parameters.placement.mode));
    end

    fprintf('Transducer placement mode: %s\n', upper(mode));

    switch mode

        % ── MANUAL ────────────────────────────────────────────────────────
        case 'manual'
            parameters = run_manual_placement(parameters);

        % ── LOCALITE ──────────────────────────────────────────────────────
        case 'localite'
            parameters = run_localite_placement(parameters);

        % ── HEURISTIC ─────────────────────────────────────────────────────
        case 'heuristic'
            parameters = run_with_record(parameters, mode, @run_heuristic_placement);

        % ── MNI ───────────────────────────────────────────────────────────
        case 'mni'
            parameters = run_mni_placement(parameters);

        % ── PLANTUS ───────────────────────────────────────────────────────
        case 'plantus'
            parameters = run_with_record(parameters, mode, @run_plantus_placement);

        % ── TUS_ENTRY ─────────────────────────────────────────────────────
        case 'tus_entry'
            parameters = run_with_record(parameters, mode, @run_tus_entry_placement);

        otherwise
            error(['Unknown placement.mode ''%s''. ' ...
                   'Use ''manual'', ''localite'', ''heuristic'', ''mni'', ''plantus'', ' ...
                   'or ''tus_entry''.'], mode);
    end

    % Resolved positions for other tools, e.g. the Nexus viewer (never errors)
    write_placement_positions(parameters, mode);

    % Poses for other tools: Brainsight, 3D Slicer, k-Plan (placement.export_formats; never errors)
    export_placement(parameters);

    % Stand-alone placement report (never errors; see generate_placement_report)
    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'report') ...
            || isempty(parameters.placement.report) || logical(parameters.placement.report)
        generate_placement_report(parameters, mode);
    end

end

% ── Placement run record ─────────────────────────────────────────────────
function parameters = run_with_record(parameters, mode, run_fn)
% Run a placement mode that has a run folder and record its result there
% (placement_run.json). placement.reuse decides when a recorded run is used
% instead of running the mode again; see PLACEMENT_REUSE_MODE.

    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, mode)
        parameters = run_fn(parameters);    % reports the missing config itself
        return
    end
    run_dir   = placement_run_dir(parameters, mode);
    signature = placement_signature(parameters, mode);
    reuse     = placement_reuse_mode(parameters);
    remedy    = sprintf(['Set placement.reuse: never to run the %s placement again, or point ' ...
        'path.placement (or placement.%s.target_name) at another run folder'], mode, mode);

    rec = [];
    if ~strcmp(reuse, 'never')
        [rec, record_file] = read_placement_record(run_dir);
    end
    if isempty(rec) && strcmp(reuse, 'always')
        error('PRESTUS:placement:noRecord', ...
            ['placement.reuse: always, but there is no recorded %s placement in\n  %s\n' ...
             'Run the placement first (placement.reuse: auto or never).'], mode, run_dir);
    end

    if ~isempty(rec)
        strict = strcmp(reuse, 'always') || overwrite_never(parameters);
        if strict
            use_rec = check_sidecar_signature(record_file, rec.signature, signature, parameters, remedy) ...
                || strcmp(cache_check_mode(parameters), 'warn');
        else
            % auto: a record made with other settings is replaced, not an error
            quiet = parameters;
            quiet.io.cache_signature_check = 'warn';
            w = warning('off', 'prestus:cache:signatureMismatch');
            restore = onCleanup(@() warning(w));
            [use_rec, diffs] = check_sidecar_signature(record_file, rec.signature, signature, quiet, remedy);
            clear restore
            if ~use_rec
                fprintf(['The %s placement recorded in\n  %s\nwas made with different settings ' ...
                    '(%s); running it again (placement.reuse: auto).\n'], ...
                    mode, record_file, strjoin(diffs, ', '));
            end
        end
        if use_rec
            fprintf('Reusing the %s placement recorded in:\n  %s\n', mode, record_file);
            for ti = 1:numel(parameters.transducer)
                k = min(ti, size(rec.trans_pos, 1));
                parameters.transducer(ti).trans_pos = rec.trans_pos(k, :);
                parameters.transducer(ti).focus_pos = rec.focus_pos(k, :);
                if ~isempty(rec.trans_pos_ras)
                    parameters.transducer(ti).trans_pos_ras = rec.trans_pos_ras(k, :);
                    parameters.transducer(ti).focus_pos_ras = rec.focus_pos_ras(k, :);
                end
            end
            if ~isempty(rec.target_ras)
                parameters.placement.(mode).target_ras = rec.target_ras;
            end
            plot_placement_t1_overlay(parameters, ...
                parameters.transducer(1).trans_pos, parameters.transducer(1).focus_pos, mode);
            return
        end
    end

    parameters = run_fn(parameters);
    write_placement_record(parameters, mode, run_dir, signature);
end

function reuse = placement_reuse_mode(parameters)
% placement.reuse: 'auto' (default) reuses a recorded run whose settings
% match and runs the mode again otherwise; 'always' requires a recorded run
% (settings checked per io.cache_signature_check); 'never' always runs it.
    reuse = 'auto';
    if isfield(parameters.placement, 'reuse') && ~isempty(parameters.placement.reuse)
        reuse = parameters.placement.reuse;
        if islogical(reuse) || isnumeric(reuse)
            % unquoted YAML: true/1 -> always, false/0 -> never
            if logical(reuse); reuse = 'always'; else; reuse = 'never'; end
        end
        reuse = lower(char(reuse));
    end
    if ~any(strcmp(reuse, {'auto', 'always', 'never'}))
        error('PRESTUS:placement:badReuse', ...
            'placement.reuse must be ''auto'', ''always'' or ''never'' (got ''%s'').', reuse);
    end
end

function tf = overwrite_never(parameters)
    tf = isfield(parameters, 'io') && isfield(parameters.io, 'overwrite_files') && ...
        strcmp(char(parameters.io.overwrite_files), 'never');
end

function mode = cache_check_mode(parameters)
    mode = 'error';
    if isfield(parameters, 'io') && isfield(parameters.io, 'cache_signature_check') && ...
            (ischar(parameters.io.cache_signature_check) || isstring(parameters.io.cache_signature_check))
        mode = lower(char(parameters.io.cache_signature_check));
    end
end

% ── Manual placement dispatcher ──────────────────────────────────────────
function parameters = run_manual_placement(parameters)
% Convert manually specified positions to T1 voxel indices if they were
% given in RAS+ or MNI mm. Each transducer keeps its own positions (unlike
% the other modes, which write one resolved pair into all transducers).

    cs = manual_coordinate_space(parameters);
    if strcmp(cs, 'voxel')
        disp('Using manually specified transducer / focus positions (T1 voxels).');
        return;
    end

    % Voxel indices refer to the T1 grid, so convert via the T1 affine
    t1_glob = fullfile(parameters.path.anat, ...
        sprintf(parameters.path.t1_pattern, parameters.subject_id));
    t1_hits = dir(t1_glob);
    if isempty(t1_hits)
        error('PRESTUS:placement:t1NotFound', 'T1 file not found matching: %s', t1_glob);
    end
    t1_header = niftiinfo(fullfile(t1_hits(1).folder, t1_hits(1).name));

    if strcmp(cs, 'mni')
        m2m = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
    end

    n_tr = numel(parameters.transducer);
    input_trans = nan(n_tr, 3);
    input_focus = nan(n_tr, 3);
    for ti = 1:n_tr
        tr = parameters.transducer(ti);
        if ~isfield(tr, 'trans_pos') || numel(tr.trans_pos) ~= 3 || ...
                ~isfield(tr, 'focus_pos') || numel(tr.focus_pos) ~= 3
            error('PRESTUS:placement:missingPositions', ...
                ['Transducer %d: trans_pos and focus_pos must be [1x3] ' ...
                 'for placement.manual.coordinate_space = ''%s''.'], ti, cs);
        end
        input_trans(ti,:) = tr.trans_pos(:)';
        input_focus(ti,:) = tr.focus_pos(:)';

        if strcmp(cs, 'mni')
            % Transducer lies outside the brain: linear 12dof (nonl unreliable).
            % Focus lies inside the brain: nonlinear warp, as the heuristic.
            trans_ras = mni2subject_coords_LDfix(input_trans(ti,:), m2m, parameters, '12dof');
            focus_ras = mni2subject_coords_LDfix(input_focus(ti,:), m2m, parameters, 'nonl');
        else
            trans_ras = input_trans(ti,:);
            focus_ras = input_focus(ti,:);
        end

        parameters.transducer(ti).trans_pos = ras_to_grid(trans_ras, t1_header);
        parameters.transducer(ti).focus_pos = ras_to_grid(focus_ras, t1_header);

        fprintf('Transducer %d (%s mm -> T1 voxels):\n', ti, upper(cs));
        fprintf('  trans_pos = [%g %g %g] -> [%d %d %d]\n', input_trans(ti,:), parameters.transducer(ti).trans_pos);
        fprintf('  focus_pos = [%g %g %g] -> [%d %d %d]\n', input_focus(ti,:), parameters.transducer(ti).focus_pos);
    end

    % Positions are voxel indices from here on; keep the inputs for the record
    parameters.placement.manual.input_trans_pos = input_trans;
    parameters.placement.manual.input_focus_pos = input_focus;
    parameters.placement.manual.coordinate_space = 'voxel';

    % The distance fallback was skipped at load time (positions were mm);
    % fill any missing focal distances now that positions are voxels
    parameters = focal_distance_calculation(parameters);
end

% ── Localite placement dispatcher ────────────────────────────────────────
function parameters = run_localite_placement(parameters)
% Resolve the recorded Localite pose(s) and write positions back into parameters.
%
% placement.localite.target / .position may list one entry per transducer
% (e.g. target: [L_MD, R_MD] for a bilateral setup): transducer i then uses
% entry i. A single entry is used for every transducer.

    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'localite')
        error('PRESTUS:localite:missingConfig', ...
            'placement.localite sub-struct is missing. See config_default.yaml for required fields.');
    end

    n_tr = numel(parameters.transducer);
    lc = parameters.placement.localite;
    targets = per_transducer_list(lc, 'target', n_tr);
    positions = per_transducer_list(lc, 'position', n_tr);

    % Resolve every pose first (a missing recording fails before the T1 is read)
    poses = cell(1, n_tr);
    params_ti = cell(1, n_tr);
    for ti = 1:n_tr
        p = parameters;
        if ~isempty(targets);   p.placement.localite.target   = targets{ti};   end
        if ~isempty(positions); p.placement.localite.position = positions{ti}; end
        [poses{ti}, info] = neuronav_localite_pose(p);
        params_ti{ti} = p;
        fprintf('Localite pose for transducer %d (%s mode) from:\n  %s\n', ti, info.mode, info.source);
    end

    % Load T1 header for voxel ↔ RAS coordinate transform
    filename_t1 = fullfile(parameters.path.anat, ...
        sprintf(parameters.path.t1_pattern, parameters.subject_id));
    t1_header = niftiinfo(filename_t1);

    for ti = 1:n_tr
        [trans_pos, focus_pos] = localite_matrix_to_positions(poses{ti}, t1_header, params_ti{ti});
        parameters.transducer(ti).trans_pos = trans_pos;
        parameters.transducer(ti).focus_pos = focus_pos;
    end
end

function vals = per_transducer_list(lc, field, n_tr)
% {} when the field holds one value (use it as is for all transducers),
% else one cell per transducer.
    vals = {};
    if ~isfield(lc, field); return; end
    v = lc.(field);
    if iscell(v)
        items = v(:)';
    elseif isstring(v) && numel(v) > 1
        items = cellstr(v(:)');
    elseif isnumeric(v) && numel(v) > 1
        items = num2cell(v(:)');
    else
        return
    end
    if numel(items) ~= n_tr
        error('PRESTUS:localite:listLength', ...
            'placement.localite.%s lists %d entries but %d transducer(s) are configured.', ...
            field, numel(items), n_tr);
    end
    vals = cellfun(@(x) char_or_num(x), items, 'UniformOutput', false);
end

function x = char_or_num(x)
    if isstring(x); x = char(x); end
end

% ── Heuristic placement dispatcher ───────────────────────────────────────
function parameters = run_heuristic_placement(parameters)
% Run the sphere-expansion heuristic (transducer_positioning_start).
% On 'matlab' platform: runs synchronously and reads the output CSV to
% update parameters. On HPC: submits a job and exits the pipeline (the
% user must re-run after the job finishes with the resolved positions).

    % Validate required heuristic config fields
    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'heuristic')
        error('PRESTUS:heuristic:missingConfig', ...
            'placement.heuristic sub-struct is missing. See config_default.yaml for required fields.');
    end
    validate_required_fields(parameters.placement.heuristic, ...
        {'mni_target_mm', 'target_name'}, ...
        'PRESTUS:heuristic:missingConfig', 'placement.heuristic');

    mni_coords  = parameters.placement.heuristic.mni_target_mm(:)';
    target_name = char(parameters.placement.heuristic.target_name);

    % Ensure io.dir_output is resolved (derived from path.sim when params were
    % built before the pipeline ran, e.g. placement-only benchmark scripts).
    parameters = resolve_io_dirs(parameters);

    % Build mni_targets struct expected by transducer_positioning
    % (it reads mni_targets.(target_name))
    mni_targets.(target_name) = mni_coords;

    % Build pn path-names struct from parameters.path
    pn.seg_path  = parameters.path.seg;
    pn.data_path = parameters.path.anat;
    pn.sim_path  = parameters.path.sim;

    % Expected output file from tp_select_heuristic_position (rows2vars txt)
    tpos_file = fullfile(placement_run_dir(parameters, 'heuristic'), ...
        sprintf('sub-%03d_%s.txt', parameters.subject_id, target_name));

    % Detect platform
    if strcmp(parameters.platform, 'auto')
        platform = hpc_detect_system();
        parameters.platform = platform;
    else
        platform = parameters.platform;
    end

    if strcmp(platform, 'matlab')
        % Run synchronously
        transducer_positioning(parameters, pn, target_name, mni_targets);

        % Read result txt and write positions back into parameters
        if ~isfile(tpos_file)
            error(['Heuristic positioning completed but output file not found:\n  %s\n' ...
                   'Check transducer_positioning output.'], tpos_file);
        end
        parameters = read_tpos_file(tpos_file, parameters);
        plot_placement_t1_overlay(parameters, ...
            parameters.transducer(1).trans_pos, ...
            parameters.transducer(1).focus_pos, 'heuristic');

    else
        % Submit HPC job — positions not yet available; exit and let user re-run
        fprintf(['Heuristic positioning submitted as HPC job.\n' ...
                 'Re-run the pipeline with placement.mode = ''manual'' after the job\n' ...
                 'completes and update trans_pos / focus_pos from:\n  %s\n'], tpos_file);
        transducer_positioning_start(parameters, pn, target_name, mni_targets);
        error('PRESTUS:placement:hpcHeuristic', ...
              'Pipeline halted: re-run after heuristic positioning job completes.');
    end
end

% ── rows2vars txt reader ──────────────────────────────────────────────────
function parameters = read_tpos_file(tpos_file, parameters)
% Read trans_pos and focus_pos from the tab-separated file written by
% tp_select_heuristic_position via rows2vars + writetable(WriteVariableNames=false).
% The file has two columns (no header): variable name | value.
% Keys of interest: trans_x, trans_y, trans_z, targ_x, targ_y, targ_z.

    T = readtable(tpos_file, 'Delimiter', '\t', 'ReadVariableNames', false, ...
                  'TextType', 'char');
    % Column 1 = variable names, Column 2 = values
    keys   = T{:,1};
    values = T{:,2};
    lookup = containers.Map(keys, values);

    required = {'trans_x','trans_y','trans_z','targ_x','targ_y','targ_z'};
    missing  = required(~cellfun(@(k) isKey(lookup, k), required));
    if ~isempty(missing)
        error(['Missing keys in heuristic positioning file:\n  %s\n' ...
               'Missing: %s'], tpos_file, strjoin(missing, ', '));
    end

    trans_pos = [lookup('trans_x'), lookup('trans_y'), lookup('trans_z')];
    focus_pos = [lookup('targ_x'),  lookup('targ_y'),  lookup('targ_z')];

    fprintf('Heuristic placement resolved:\n');
    fprintf('  trans_pos = [%d %d %d]\n', trans_pos);
    fprintf('  focus_pos = [%d %d %d]\n', focus_pos);

    for ti = 1:numel(parameters.transducer)
        parameters.transducer(ti).trans_pos = trans_pos;
        parameters.transducer(ti).focus_pos = focus_pos;
    end
end

% ── MNI placement dispatcher ─────────────────────────────────────────────
function parameters = run_mni_placement(parameters)
% Place the transducer and focus directly from MNI coordinates.
%
% The focus (an in-brain target) is converted with the nonlinear (nonl)
% SimNIBS transform, identical to the heuristic. The transducer point lies
% outside the brain, so it is converted with the linear (12dof) transform
% (where the nonlinear warp is unreliable), then snapped to the nearest scalp
% voxel and offset outward by placement.skin_gap_mm using the heuristic's
% standoff geometry (transducer_scalp_geometry). No search is performed.

    % Validate required config
    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'mni')
        error('PRESTUS:placement:missingConfig', ...
            'placement.mni sub-struct is missing. See config_default.yaml for required fields.');
    end
    mni = parameters.placement.mni;
    is_mni_xyz = @(f) isfield(mni, f) && isnumeric(mni.(f)) && numel(mni.(f)) == 3 ...
        && all(isfinite(mni.(f)));
    if ~is_mni_xyz('trans_pos_mm')
        error(['placement.mni.trans_pos_mm must be a finite [1x3] MNI coordinate in mm ' ...
               'for mode=''mni''.']);
    end
    if ~is_mni_xyz('focus_pos_mm')
        error(['placement.mni.focus_pos_mm must be a finite [1x3] MNI coordinate in mm ' ...
               'for mode=''mni''.']);
    end
    skin_gap_mm = 5;
    if isfield(parameters.placement, 'skin_gap_mm') && ~isempty(parameters.placement.skin_gap_mm)
        skin_gap_mm = parameters.placement.skin_gap_mm;
    end

    % Load segmentation image + header (same source as the heuristic)
    m2m = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
    seg_file   = fullfile(m2m, 'final_tissues.nii.gz');
    img        = niftiread(seg_file);
    img_header = niftiinfo(seg_file);
    pixel_size = mean(img_header.PixelDimensions(1:3));

    % ── Coordinate conversion (per-point transform type) ──────────────────
    % Transducer scalp point -> LINEAR 12dof (outside brain; nonl unreliable).
    % transform_coordinates hardcodes nonl in its 'mni' branch, so do the
    % MNI->RAS+ step explicitly with '12dof', then RAS+->grid via the affine.
    trans_ras = mni2subject_coords_LDfix(mni.trans_pos_mm(:)', m2m, parameters, '12dof');
    trans_vox = transform_coordinates(parameters, trans_ras, 'ras_plus', 'grid', img_header);

    % Focus -> NONLINEAR nonl, identical to the heuristic target.
    focus_pos = transform_coordinates(parameters, mni.focus_pos_mm(:)', 'mni', 'grid', img_header);
    if any(~isfinite(trans_vox)) || any(~isfinite(focus_pos))
        error('PRESTUS:placement:mniUnmapped', ['MNI placement: SimNIBS could not map %s to this subject ' ...
            '(outside the warp field of %s). Check the MNI coordinate.'], ...
            strjoin([repmat({'placement.mni.trans_pos_mm'}, 1, any(~isfinite(trans_vox))), ...
                     repmat({'placement.mni.focus_pos_mm'}, 1, any(~isfinite(focus_pos)))], ' and '), m2m);
    end

    % ── Snap transducer to scalp, then apply standoff geometry ────────────
    scalp  = tp_scalp_boundary(img);                 % [N x 3] outer-boundary voxels
    [~, j] = min(vecnorm(scalp - trans_vox, 2, 2));  % nearest scalp voxel
    trans_pos = transducer_scalp_geometry(scalp(j,:), focus_pos, ...
                    parameters.transducer(1), pixel_size, skin_gap_mm);
    trans_pos = round(trans_pos);

    % Bowl-to-focus distance this placement asks for. A user-picked MNI
    % scalp/focus pair can land far from the transducer's configured focal
    % distance, which would otherwise go unnoticed.
    bowl_focus_mm = norm(trans_pos - focus_pos) * pixel_size;

    fprintf('MNI placement resolved (skin_gap_mm = %.1f):\n', skin_gap_mm);
    fprintf('  trans_pos = [%d %d %d]\n', trans_pos);
    fprintf('  focus_pos = [%d %d %d]\n', focus_pos);
    fprintf('  bowl-to-focus distance = %.1f mm\n', bowl_focus_mm);
    if isfield(parameters.transducer(1), 'focal_distance_bowl') ...
            && isnumeric(parameters.transducer(1).focal_distance_bowl) ...
            && isscalar(parameters.transducer(1).focal_distance_bowl)
        fd_bowl = parameters.transducer(1).focal_distance_bowl;
        if abs(bowl_focus_mm - fd_bowl) > 5
            warn(['MNI placement: bowl-to-focus distance %.1f mm differs from ' ...
                  'transducer.focal_distance_bowl %.1f mm by more than 5 mm.'], ...
                  bowl_focus_mm, fd_bowl);
        end
    end

    for ti = 1:numel(parameters.transducer)
        parameters.transducer(ti).trans_pos = trans_pos;
        parameters.transducer(ti).focus_pos = focus_pos;
    end
    plot_placement_t1_overlay(parameters, ...
        parameters.transducer(1).trans_pos, ...
        parameters.transducer(1).focus_pos, 'mni');
end

% ── PlanTUS placement dispatcher ─────────────────────────────────────────
function parameters = run_plantus_placement(parameters)
% Call position_transducer_plantus and write results back into parameters.

    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'plantus')
        error('PRESTUS:plantus:missingConfig', ...
            'placement.plantus sub-struct is missing. See config_default.yaml for required fields.');
    end
    % Validate before reading the T1, so config errors surface first
    % (position_transducer_plantus re-checks the same fields).
    validate_required_fields(parameters.placement.plantus, ...
        {'script_path', 'env_path', 'focal_distance_list', 'flhm_list', 'target_name'}, ...
        'PRESTUS:plantus:missingConfig', 'placement.plantus');

    filename_t1 = fullfile(parameters.path.anat, ...
        sprintf(parameters.path.t1_pattern, parameters.subject_id));
    t1_header = niftiinfo(filename_t1);

    [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras, target_ras] = ...
        position_transducer_plantus(parameters, t1_header);

    for ti = 1:numel(parameters.transducer)
        parameters.transducer(ti).trans_pos     = trans_pos;
        parameters.transducer(ti).focus_pos     = focus_pos;
        parameters.transducer(ti).trans_pos_ras = trans_pos_ras;
        parameters.transducer(ti).focus_pos_ras = focus_pos_ras;
    end
    parameters.placement.plantus.target_ras = target_ras;
    plot_placement_t1_overlay(parameters, ...
        parameters.transducer(1).trans_pos, ...
        parameters.transducer(1).focus_pos, 'plantus');

    pt_target_name = char(parameters.placement.plantus.target_name);
    traj_path = fullfile(placement_run_dir(parameters, 'plantus'), ...
        sprintf('sub-%03d_%s_trajectory.txt', parameters.subject_id, pt_target_name));
    export_brainsight_trajectory(trans_pos_ras, focus_pos_ras, pt_target_name, traj_path);
end

function parameters = run_tus_entry_placement(parameters)
% Call position_transducer_tus_entry and write results back into parameters.

    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'tus_entry')
        error('PRESTUS:tus_entry:missingConfig', ...
            'placement.tus_entry sub-struct is missing. See config_default.yaml for required fields.');
    end

    filename_t1 = fullfile(parameters.path.anat, ...
        sprintf(parameters.path.t1_pattern, parameters.subject_id));
    t1_header = niftiinfo(filename_t1);

    [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
        position_transducer_tus_entry(parameters, t1_header);

    for ti = 1:numel(parameters.transducer)
        parameters.transducer(ti).trans_pos     = trans_pos;
        parameters.transducer(ti).focus_pos     = focus_pos;
        parameters.transducer(ti).trans_pos_ras = trans_pos_ras;
        parameters.transducer(ti).focus_pos_ras = focus_pos_ras;
    end
    plot_placement_t1_overlay(parameters, ...
        parameters.transducer(1).trans_pos, ...
        parameters.transducer(1).focus_pos, 'tus_entry');

    te = parameters.placement.tus_entry;
    tus_target_name = 'target';
    if isfield(te, 'target_name') && ~isempty(te.target_name)
        tus_target_name = char(te.target_name);
    end
    traj_path = fullfile(placement_run_dir(parameters, 'tus_entry'), ...
        sprintf('sub-%03d_%s_trajectory.txt', parameters.subject_id, tus_target_name));
    export_brainsight_trajectory(trans_pos_ras, focus_pos_ras, tus_target_name, traj_path);
end
