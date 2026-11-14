function [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras, target_ras] = ...
    position_transducer_plantus(parameters, t1_header)
% POSITION_TRANSDUCER_PLANTUS  Optimise transducer placement using PlanTUS.
%
% Calls the external PlanTUS tool (Lueckel et al., Mainz) to score candidate
% transducer positions on the scalp against five geometric objectives:
%   - Maximise overlap between focal beam and target
%   - Minimise beam angle at skin surface
%   - Minimise skin–skull angle
%   - Minimise skull thickness at entry
%   - Minimise transducer-to-target distance
%
% The function:
%   1. Builds a PlanTUS transducer configuration YAML from PRESTUS parameters.
%   2. Creates a spherical target mask around the focal point.
%   3. Runs a Python subprocess — mode depends on placement.plantus.interactive
%      ([] = follow parameters.simulation.interactive):
%        interactive = 0 (default / HPC):
%          prestus_plantus_launcher.py — headless, automatically selects the
%          optimal skin vertex via weighted composite scoring.
%        interactive = 1 (desktop):
%          prestus_plantus_gui_launcher.py — computes the same scoring maps then
%          opens Connectome Workbench so the user can inspect them visually,
%          click a vertex, see the transducer there and confirm or choose
%          again in a dialog. Falls back to headless without a display.
%   4. Loads the output *Localite.mat (position_matrix) written by PlanTUS.
%   5. Converts the 4×4 RAS matrix to voxel positions via
%      localite_matrix_to_positions (shared with the 'localite' mode).
%
% The target position can be supplied either as:
%   (a) parameters.placement.plantus.mni_target_mm  — MNI coordinates [mm];
%       converted to subject space via SimNIBS mni2subject_coords.
%   (b) parameters.transducer(1).focus_pos in voxel space (manual pre-set) —
%       converted to RAS via the T1 header affine for mask creation.
%
% Use as:
%   [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
%       position_transducer_plantus(parameters, t1_header)
%
% Input:
%   parameters  - Simulation config struct with placement.plantus sub-struct
%                 (see config_default.yaml for all fields).
%   t1_header   - NIfTI header with Transform field (from niftiinfo), used for
%                 voxel ↔ RAS coordinate transforms.
%
% Output:
%   trans_pos      - [1×3] bowl rear-centre position in voxel coordinates (ijk)
%   focus_pos      - [1×3] acoustic focus position in voxel coordinates (ijk)
%   trans_pos_ras  - [1×3] bowl rear-centre in RAS space (mm)
%   focus_pos_ras  - [1×3] acoustic focus in RAS space (mm)
%   target_ras     - [1×3] intended target in RAS space (mm)
%
% Required configuration (parameters.placement.plantus):
%   script_path        - Path to PlanTUS installation dir (must contain code/PlanTUS.py)
%   env_path   - Path to SimNIBS installation (must contain simnibs_env/); empty =
%                the Python of the install in startup.simnibs_bin_path
%   focal_distance_list - [1×N] calibrated focal distances [mm]
%   flhm_list           - [1×N] full-width at half-max focal lengths [mm]
%   target_name         - Label string for output filenames
%
% Optional configuration (parameters.placement.plantus):
%   mni_target_mm           - [1×3] MNI target coordinate [mm]
%   mask_radius_mm          - Sphere radius for target mask [mm] (default 2)
%   max_angle_deg           - Max transducer tilt [degrees] (default 10)
%   min_distance_mm         - Skin–target distance window and optimum [mm]
%   max_distance_mm           used for eligibility and the distance criterion.
%   optimal_distance_mm       Defaults: min/max of focal_distance_list -/+
%                             max(flhm_list), and mean(focal_distance_list).
%                             Set them for a steerable array, where the
%                             default window grows with the steering range.
%   weight_skin_target_distances    - Objective weight (default 0.2)
%   weight_skin_target_angles       - Objective weight (default 0.2)
%   weight_skin_target_intersections- Objective weight (default 0.2)
%   weight_skin_skull_angles        - Objective weight (default 0.2)
%   weight_skull_thickness          - Objective weight (default 0.2); bone path
%                                     length along the skin normal, from
%                                     m2m_sub-XXX/final_tissues.nii.gz
%   interactive             - 1/true: choose the vertex in Connectome Workbench;
%                             0/false: headless; [] (default): follow
%                             simulation.interactive
%   connectome_wb_path      - Path to Connectome Workbench bin/ (default '')
%   save_metrics            - Write per-vertex criteria, composite score and
%                             selected vertex to <work_dir>/plantus_metrics.mat
%                             (default true; both modes)
%   evaluate                - Evaluate the placement from plantus_metrics.mat
%                             (weight sweep, diagnostics, figures; see
%                             plantus_evaluation.py) for the placement report
%                             (default true; needs save_metrics)
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, LOCALITE_MATRIX_TO_POSITIONS,
%           POSITION_TRANSDUCER_LOCALITE

    arguments
        parameters (1,1) struct
        t1_header  (1,1) struct
    end

    pt = parameters.placement.plantus;

    % ── Validate required fields ──────────────────────────────────────────

    validate_required_fields(pt, ...
        {'script_path', 'focal_distance_list', 'flhm_list', 'target_name'}, ...
        'PRESTUS:plantus:missingConfig', 'placement.plantus');

    script_path      = char(pt.script_path);
    % env_path empty: the Python of the SimNIBS install in startup.simnibs_bin_path
    env_path = '';
    if isfield(pt, 'env_path') && ~isempty(pt.env_path)
        env_path = char(pt.env_path);
    end
    target_name      = char(pt.target_name);
    focal_dist_list  = pt.focal_distance_list(:)';
    flhm_list        = pt.flhm_list(:)';

    if ~isfolder(script_path)
        error('PRESTUS:plantus:notFound', ...
            'placement.plantus.script_path does not exist:\n  %s', script_path);
    end
    if ~isfile(fullfile(script_path, 'code', 'PlanTUS.py'))
        error('PRESTUS:plantus:notFound', ...
            ['PlanTUS.py not found under placement.plantus.script_path/code/:\n  %s\n' ...
             'Ensure PlanTUS is installed at that path (git clone https://github.com/mlueckel/PlanTUS).'], ...
            script_path);
    end

    % ── Log PlanTUS version ───────────────────────────────────────────────
    % run_command merges stderr into stdout (like system()'s combined
    % output), so a failing git call (not a repo, no exact-match tag) must
    % be judged by its exit status, not just by cmdout being non-empty --
    % otherwise a git error message would be mistaken for a version string.
    [hash_status, plantus_hash] = run_command({'git', '-C', script_path, 'rev-parse', '--short', 'HEAD'});
    if hash_status ~= 0
        plantus_hash = '';
    end
    plantus_hash = strtrim(plantus_hash);
    [tag_status, plantus_tag] = run_command({'git', '-C', script_path, 'describe', '--tags', '--exact-match'});
    if tag_status ~= 0
        plantus_tag = '';
    end
    plantus_tag = strtrim(plantus_tag);
    if ~isempty(plantus_tag)
        plantus_ver = plantus_tag;
    elseif ~isempty(plantus_hash)
        plantus_ver = plantus_hash;
    else
        plantus_ver = 'unknown';
    end
    fprintf('PlanTUS version: %s\n', plantus_ver);

    % ── Resolve target position in RAS (mm) ──────────────────────────────

    target_ras = resolve_target_ras(parameters, t1_header);
    fprintf('PlanTUS: target position (RAS mm): [%.1f %.1f %.1f]\n', target_ras);
    fprintf(['PlanTUS citation: Lueckel et al. (2025) Brain Stimulation 18(5):1563-1565. ' ...
             'https://doi.org/10.1016/j.brs.2025.08.013\n']);

    % ── Paths ─────────────────────────────────────────────────────────────

    filename_t1 = fullfile(parameters.path.anat, ...
        sprintf(parameters.path.t1_pattern, parameters.subject_id));

    m2m_dir  = fullfile(parameters.path.seg, ...
        sprintf('m2m_sub-%03d', parameters.subject_id));
    mesh_path = fullfile(m2m_dir, sprintf('sub-%03d.msh', parameters.subject_id));

    if ~isfile(mesh_path)
        error('PRESTUS:plantus:notFound', ...
            ['SimNIBS mesh not found at expected path:\n  %s\n' ...
             'Check parameters.path.seg and subject_id.'], mesh_path);
    end

    work_dir  = placement_run_dir(parameters, 'plantus');
    if ~isfolder(work_dir), mkdir(work_dir); end

    mask_path     = fullfile(work_dir, sprintf('sub-%03d_%s_mask.nii.gz', ...
                        parameters.subject_id, target_name));
    tx_config_path = fullfile(work_dir, 'PlanTUSTxConfig.yaml');

    % ── Build target mask (sphere in T1 space) ────────────────────────────

    mask_radius_mm = 2;
    if isfield(pt, 'mask_radius_mm') && ~isempty(pt.mask_radius_mm)
        mask_radius_mm = pt.mask_radius_mm;
    end

    write_target_mask(filename_t1, target_ras, mask_path, mask_radius_mm);
    fprintf('PlanTUS: target mask written:\n  %s\n', mask_path);

    % ── Build PlanTUS transducer config YAML ─────────────────────────────

    tx_yaml = build_plantus_yaml(parameters, focal_dist_list, flhm_list, target_name, pt);
    write_yaml(tx_config_path, tx_yaml);
    fprintf('PlanTUS: transducer config written:\n  %s\n', tx_config_path);

    % ── Call PlanTUS ──────────────────────────────────────────────────────

    % Probe several layouts in order of specificity:
    %   1. Official installer:  <root>/simnibs_env/bin/python
    %   2. Windows installer:   <root>/simnibs_env/Scripts/python.exe
    %   3. Conda env (bin/ given directly): <path>/python
    %   4. Conda env (env root given):      <path>/bin/python
    %   (env_path empty: the interpreter charm of startup.simnibs_bin_path runs)
    if isempty(env_path)
        candidates = {''};
        if isfield(parameters, 'startup') && isfield(parameters.startup, 'simnibs_bin_path') ...
                && ~isempty(parameters.startup.simnibs_bin_path)
            info = simnibs_install_info(char(parameters.startup.simnibs_bin_path));
            candidates = {info.python};
        end
    else
        candidates = { ...
            fullfile(env_path, 'simnibs_env', 'bin',     'python'), ...
            fullfile(env_path, 'simnibs_env', 'Scripts', 'python.exe'), ...
            fullfile(env_path, 'python'), ...
            fullfile(env_path, 'python.exe'), ...
            fullfile(env_path, 'bin', 'python') ...
        };
    end
    python_bin = '';
    for i_cand = 1:numel(candidates)
        if isfile(candidates{i_cand})
            python_bin = candidates{i_cand};
            break;
        end
    end
    if isempty(python_bin)
        error('PRESTUS:plantus:notFound', ...
            ['Python binary not found. Tried:\n' ...
             '  %s\n' ...
             'Leave placement.plantus.env_path empty to use the SimNIBS install of\n' ...
             'startup.simnibs_bin_path, or set it to the SimNIBS installation\n' ...
             'root (e.g. /opt/SimNIBS-4.1) or to the conda env bin/ directory\n' ...
             '(e.g. ~/.conda/envs/simnibs_v4.6.0/bin).'], ...
            strjoin(candidates, '\n  '));
    end

    % Choose headless (automatic) or interactive (GUI) launcher.
    % Interactive mode opens Connectome Workbench so the user can inspect
    % the five scoring maps and manually confirm a transducer position.
    % Headless mode automatically selects the best vertex via composite score.
    %
    % placement.plantus.launcher_py_override (undocumented): when set, overrides the
    % launcher path directly — used by the test suite to inject mock scripts.
    if isfield(pt, 'interactive') && ~isempty(pt.interactive)
        use_gui = logical(pt.interactive);
    else
        use_gui = isfield(parameters, 'simulation') && ...
                  isfield(parameters.simulation, 'interactive') && ...
                  logical(parameters.simulation.interactive);
    end
    if use_gui && isunix && ~ismac && isempty(getenv('DISPLAY')) && isempty(getenv('WAYLAND_DISPLAY'))
        warning('PRESTUS:plantus:noDisplay', ...
            'No display for Connectome Workbench (e.g. a cluster job); PlanTUS runs headless.');
        use_gui = false;
    end

    if isfield(pt, 'launcher_py_override') && ~isempty(pt.launcher_py_override)
        launcher = char(pt.launcher_py_override);
        if ~isfile(launcher)
            error('PRESTUS:plantus:notFound', ...
                'Override launcher not found:\n  %s', launcher);
        end
    else
        fn_dir = fileparts(mfilename('fullpath'));

        if use_gui
            launcher = fullfile(fn_dir, 'prestus_plantus_gui_launcher.py');
            if ~isfile(launcher)
                error('PRESTUS:plantus:notFound', ...
                    ['GUI shim not found:\n  %s\n' ...
                     'Ensure prestus_plantus_gui_launcher.py exists alongside ' ...
                     'position_transducer_plantus.m.'], launcher);
            end
            fprintf(['PlanTUS: interactive mode — Connectome Workbench will open.\n' ...
                     '  Click a skin vertex (click again to change it).\n' ...
                     '  Close wb_view to place the last selected vertex, then use it\n' ...
                     '  or choose another in the review dialog.\n']);
        else
            launcher = fullfile(fn_dir, 'prestus_plantus_launcher.py');
            if ~isfile(launcher)
                error('PRESTUS:plantus:notFound', ...
                    ['Headless launcher not found:\n  %s\n' ...
                     'Ensure prestus_plantus_launcher.py exists alongside ' ...
                     'position_transducer_plantus.m.'], launcher);
            end
        end
    end

    plantus_argv = {python_bin, launcher, filename_t1, mesh_path, mask_path, ...
        tx_config_path, '--plantus_root', script_path, '--output_dir', work_dir};

    fprintf('PlanTUS: running placement optimisation...\n');
    % run_command invokes python directly via argv, with no shell at all --
    % every path here (config-derived, but still worth removing the shell
    % rather than trusting quoting) reaches the process as literal
    % argument values regardless of their contents.
    %
    % In interactive mode the launcher's output is echoed live, so the user
    % sees the wb_view instructions and each selected vertex while the
    % window is open. The child's stdin is never connected, so the launcher
    % is told to confirm in a review dialog (Qt, shipped with SimNIBS) rather
    % than by a console prompt; without Qt the last selected vertex is used.
    if use_gui
        plantus_argv = [plantus_argv, {'--selection', 'review'}];
    end
    % Snapshot outputs from earlier runs so only a file this run wrote is used.
    previous_outputs = dir(fullfile(work_dir, '**', '*Localite.mat'));
    [status, cmdout] = run_command(plantus_argv, 'Echo', use_gui);

    if status ~= 0
        error('PRESTUS:plantus:failed', ...
            'PlanTUS subprocess exited with status %d.\nOutput:\n%s', status, cmdout);
    end
    if ~use_gui
        fprintf('%s\n', cmdout);
    end

    % ── Load output position matrix ───────────────────────────────────────

    % Both launchers write into --output_dir (work_dir, one per target), with
    % PlanTUS adding a vtx<N>/ subfolder per placed vertex.
    localite_mat = find_plantus_output(work_dir, previous_outputs);
    fprintf('PlanTUS: loading result:\n  %s\n', localite_mat);

    result = load(localite_mat);
    if ~isfield(result, 'position_matrix')
        error('PRESTUS:plantus:badOutput', ...
            'Expected ''position_matrix'' variable in PlanTUS output:\n  %s', localite_mat);
    end

    coord_matrix = double(result.position_matrix);   % [4×4] RAS mm

    % ── Convert 4×4 RAS matrix → voxel positions (shared with localite) ──
    %
    % PlanTUS places column 4 (origin) at the bowl rear centre — it applies
    % plane_offset itself when building the matrix.  Tell localite_matrix_to_positions
    % that the origin is already at the bowl rear centre so it does not apply
    % the offset a second time (tracker_to_bowl_mm = 0).
    parameters_for_pos = parameters;
    if ~isfield(parameters_for_pos, 'placement')
        parameters_for_pos.placement = struct();
    end
    if ~isfield(parameters_for_pos.placement, 'localite')
        parameters_for_pos.placement.localite = struct();
    end
    parameters_for_pos.placement.localite.tracker_to_bowl_mm = 0;

    [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
        localite_matrix_to_positions(coord_matrix, t1_header, parameters_for_pos);

    fprintf('PlanTUS placement resolved:\n');
    fprintf('  trans_pos (vox): [%d %d %d]\n', trans_pos);
    fprintf('  focus_pos (vox): [%d %d %d]\n', focus_pos);
end

% ── Resolve target in RAS space ───────────────────────────────────────────
function target_ras = resolve_target_ras(parameters, t1_header)
% Return target (focus) position in RAS mm. Tries, in order:
%   1. placement.plantus.mni_target_mm  → subject RAS via SimNIBS
%   2. transducer(1).focus_pos (voxel)  → RAS via T1 header affine

    pt = parameters.placement.plantus;

    if isfield(pt, 'mni_target_mm') && ~isempty(pt.mni_target_mm)
        mni_mm = pt.mni_target_mm(:)';
        target_ras = transform_coordinates(parameters, mni_mm, 'mni', 'ras_plus', t1_header);
        return;
    end

    % Fall back to existing focus_pos (voxel → RAS)
    td = parameters.transducer(1);
    if isfield(td, 'focus_pos') && ~isempty(td.focus_pos)
        target_ras = transform_coordinates(parameters, td.focus_pos(:)', 'grid', 'ras_plus', t1_header);
        return;
    end

    error('PRESTUS:plantus:noTarget', ...
        ['Cannot determine target position for PlanTUS mask. ' ...
         'Set placement.plantus.mni_target_mm or provide transducer.focus_pos.']);
end

% ── Write spherical target mask as NIfTI ─────────────────────────────────
function write_target_mask(t1_path, target_ras, mask_path, radius_mm)
% Create a binary NIfTI sphere centred on target_ras with given radius.

    info = niftiinfo(t1_path);
    vol  = zeros(info.ImageSize, 'single');

    center_vox = ras_to_grid(target_ras(:)', info);
    ps = info.PixelDimensions(1:3);

    [X, Y, Z] = ndgrid(1:info.ImageSize(1), 1:info.ImageSize(2), 1:info.ImageSize(3));
    dist2 = ((X - center_vox(1)) .* ps(1)).^2 + ...
            ((Y - center_vox(2)) .* ps(2)).^2 + ...
            ((Z - center_vox(3)) .* ps(3)).^2;
    vol(dist2 <= radius_mm^2) = 1;

    % The header comes from the T1, whose datatype (often int16/uint8) and
    % intensity scaling must not carry over: niftiwrite requires the volume
    % class to match info.Datatype, and scl_slope/inter would rescale the mask.
    info.Datatype              = 'single';
    info.BitsPerPixel          = 32;
    info.MultiplicativeScaling = 1;
    info.AdditiveOffset        = 0;

    nii_path = regexprep(mask_path, '\.nii\.gz$', '.nii');
    niftiwrite(vol, nii_path, info);
    gzip(nii_path);
    delete(nii_path);
end

% ── Build PlanTUS YAML struct ─────────────────────────────────────────────
function tx_yaml = build_plantus_yaml(parameters, focal_dist_list, flhm_list, target_name, pt)
% Map PRESTUS transducer parameters to the flat PlanTUS YAML format.

    td = parameters.transducer(1);

    % Transducer diameter: outermost element outer diameter (annular) or
    % explicit parameters.transducer(1).diameter_mm if set.
    if isfield(td, 'diameter_mm') && ~isempty(td.diameter_mm)
        tx_diameter = td.diameter_mm;
    elseif isfield(td, 'annular') && isfield(td.annular, 'elem_od_mm')
        tx_diameter = max(td.annular.elem_od_mm);
    else
        error('PRESTUS:plantus:noDiameter', ...
            ['Cannot determine transducer diameter for PlanTUS config. ' ...
             'Set transducer.diameter_mm or ensure transducer.annular.elem_od_mm is loaded.']);
    end

    % plane_offset: gap between radiating surface and exit plane
    % = curv_radius_mm - dist_geom_ep_mm (same geometry used in localite mode)
    plane_offset = 0;
    if isfield(td, 'annular')
        plane_offset = annular_bowl_depth_mm(td.annular);
    end

    additional_offset = 5;      % gap between skin and exit plane (placement.skin_gap_mm)
    if isfield(parameters.placement, 'skin_gap_mm') && ~isempty(parameters.placement.skin_gap_mm)
        additional_offset = parameters.placement.skin_gap_mm;
    end

    max_angle = 10;
    if isfield(pt, 'max_angle_deg') && ~isempty(pt.max_angle_deg)
        max_angle = pt.max_angle_deg;
    end

    connectome_wb_path = '';
    if isfield(pt, 'connectome_wb_path') && ~isempty(pt.connectome_wb_path)
        connectome_wb_path = char(pt.connectome_wb_path);
    end

    save_metrics = true;
    if isfield(pt, 'save_metrics') && ~isempty(pt.save_metrics)
        save_metrics = logical(pt.save_metrics);
    end
    evaluate = true;
    if isfield(pt, 'evaluate') && ~isempty(pt.evaluate)
        evaluate = logical(pt.evaluate);
    end

    weights = struct( ...
        'skin_target_distances',     get_weight(pt, 'weight_skin_target_distances',     0.2), ...
        'skin_target_angles',        get_weight(pt, 'weight_skin_target_angles',        0.2), ...
        'skin_target_intersections', get_weight(pt, 'weight_skin_target_intersections', 0.2), ...
        'skin_skull_angles',         get_weight(pt, 'weight_skin_skull_angles',         0.2), ...
        'skull_thickness',           get_weight(pt, 'weight_skull_thickness',           0.2)  ...
    );

    % Older GUI versions saved the weights as placement.plantus.weights.<name>,
    % which was never read. Honour them, since they are what the user set.
    if isfield(pt, 'weights') && isstruct(pt.weights)
        legacy = intersect(fieldnames(pt.weights), fieldnames(weights));
        for k = 1:numel(legacy)
            v = pt.weights.(legacy{k});
            if ~isempty(v)
                weights.(legacy{k}) = v;
            end
        end
        if ~isempty(legacy)
            warning('position_transducer_plantus:legacyWeights', ...
                ['Using placement.plantus.weights.* (saved by an older GUI) over ', ...
                 'placement.plantus.weight_*. Rename them to weight_<name> in the config.']);
        end
    end

    % min_distance: allow up to max(flhm_list) below the closest focal depth so
    % that the min_distance filter does not exclude all vertices when focal depths
    % are tightly clustered (e.g. single-element transducer with one focal depth).
    flhm_margin = max(flhm_list);
    if isnan(flhm_margin) || flhm_margin <= 0
        flhm_margin = 10;   % 10 mm fallback when FLHM unavailable
    end
    min_distance     = get_weight(pt, 'min_distance_mm',     max(0, min(focal_dist_list) - flhm_margin));
    max_distance     = get_weight(pt, 'max_distance_mm',     max(focal_dist_list) + flhm_margin);
    optimal_distance = get_weight(pt, 'optimal_distance_mm', mean(focal_dist_list));
    tx_yaml = struct( ...
        'min_distance',                    min_distance, ...
        'max_distance',                    max_distance, ...
        'optimal_distance',                optimal_distance, ...
        'transducer_diameter',             tx_diameter, ...
        'max_angle',                       max_angle, ...
        'plane_offset',                    plane_offset, ...
        'additional_offset',               additional_offset, ...
        'focal_distance_list',             focal_dist_list, ...
        'flhm_list',                       flhm_list, ...
        'IDTarget',                        target_name, ...
        'connectome_wb_path',              connectome_wb_path, ...
        'weight_skin_target_distances',    weights.skin_target_distances, ...
        'weight_skin_target_angles',       weights.skin_target_angles, ...
        'weight_skin_target_intersections',weights.skin_target_intersections, ...
        'weight_skin_skull_angles',        weights.skin_skull_angles, ...
        'weight_skull_thickness',          weights.skull_thickness, ...
        'save_metrics',                    save_metrics, ...
        'evaluate',                        evaluate, ...
        'bUseGenericTransducerModel',      false ...
    );
end

function w = get_weight(pt, field, default_val)
    if isfield(pt, field) && ~isempty(pt.(field))
        w = pt.(field);
    else
        w = default_val;
    end
end

% ── Write minimal YAML ────────────────────────────────────────────────────
function write_yaml(fname, s)
% Write a scalar struct to a flat YAML file.
% Uses yaml.dump (MATLAB R2023a+) when available; falls back to manual write.

    if exist('yaml', 'class') == 8
        yaml.dumpFile(fname, s);
        return;
    end

    fid = fopen(fname, 'w');
    if fid == -1
        error('PRESTUS:plantus:ioError', 'Cannot write YAML to:\n  %s', fname);
    end
    cleanup = onCleanup(@() fclose(fid));

    fields = fieldnames(s);
    for i = 1:numel(fields)
        k = fields{i};
        v = s.(k);
        if islogical(v)
            fprintf(fid, '%s: %s\n', k, mat2str(v));
        elseif ischar(v) || isstring(v)
            fprintf(fid, '%s: "%s"\n', k, char(v));
        elseif isscalar(v) && ~endsWith(k, '_list')
            fprintf(fid, '%s: %g\n', k, v);
        else
            % Vector → YAML inline list
            fprintf(fid, '%s: [%s]\n', k, strjoin(arrayfun(@(x) num2str(x), v, 'UniformOutput', false), ', '));
        end
    end
end

% ── Locate PlanTUS output mat file ────────────────────────────────────────
function mat_path = find_plantus_output(work_dir, previous)
% Search work_dir recursively for the *Localite.mat written by this run.
% Files listed in previous (a dir() snapshot taken before the launch) are
% ignored unless they were modified since, so a result left over from an
% earlier run is never picked up in place of a missing one. Errors if no
% new file is found; takes the most recent if there are several.

    hits = dir(fullfile(work_dir, '**', '*Localite.mat'));

    if ~isempty(hits) && ~isempty(previous)
        prev_keys = strcat({previous.folder}, filesep, {previous.name}, '|', ...
                           arrayfun(@(h) num2str(h.datenum, 12), previous, 'UniformOutput', false)');
        hit_keys  = strcat({hits.folder}, filesep, {hits.name}, '|', ...
                           arrayfun(@(h) num2str(h.datenum, 12), hits, 'UniformOutput', false)');
        hits = hits(~ismember(hit_keys, prev_keys));
    end

    if isempty(hits)
        error('PRESTUS:plantus:noOutput', ...
            ['PlanTUS did not produce a new *Localite.mat file in:\n  %s\n' ...
             'Check the PlanTUS output above for errors.'], work_dir);
    end

    if numel(hits) > 1
        % Take the most recently modified one
        [~, idx] = max([hits.datenum]);
        hits = hits(idx);
        warning('PRESTUS:plantus:multipleOutputs', ...
            'Multiple new *Localite.mat files found; using the most recent:\n  %s', ...
            fullfile(hits.folder, hits.name));
    end

    mat_path = fullfile(hits.folder, hits.name);
end
