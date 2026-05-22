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
%   3. Runs a Python subprocess — mode depends on parameters.simulation.interactive:
%        interactive = 0 (default / HPC):
%          prestus_plantus_launcher.py — headless, automatically selects the
%          optimal skin vertex via weighted composite scoring.
%        interactive = 1 (desktop):
%          prestus_plantus_gui_launcher.py — computes the same scoring maps then
%          opens Connectome Workbench so the user can inspect them visually,
%          click a vertex and confirm placement interactively.
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
%   script_path        - Path to PlanTUS installation dir (must contain PlanTUS_wrapper.py)
%   env_path   - Path to SimNIBS installation (must contain simnibs_env/)
%   focal_distance_list - [1×N] calibrated focal distances [mm]
%   flhm_list           - [1×N] full-width at half-max focal lengths [mm]
%   target_name         - Label string for output filenames
%
% Optional configuration (parameters.placement.plantus):
%   mni_target_mm           - [1×3] MNI target coordinate [mm]
%   mask_radius_mm          - Sphere radius for target mask [mm] (default 2)
%   additional_offset_mm    - Standoff between skin and exit plane [mm] (default 0)
%   max_angle_deg           - Max transducer tilt [degrees] (default 10)
%   weight_skin_target_distances    - Objective weight (default 0.2)
%   weight_skin_target_angles       - Objective weight (default 0.2)
%   weight_skin_target_intersections- Objective weight (default 0.2)
%   weight_skin_skull_angles        - Objective weight (default 0.2)
%   weight_skull_thickness          - Objective weight (default 0.2)
%   connectome_wb_path      - Path to Connectome Workbench bin/ (default '')
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, LOCALITE_MATRIX_TO_POSITIONS,
%           POSITION_TRANSDUCER_LOCALITE

    arguments
        parameters (1,1) struct
        t1_header  (1,1) struct
    end

    pt = parameters.placement.plantus;

    % ── Validate required fields ──────────────────────────────────────────

    required = {'script_path', 'env_path', 'focal_distance_list', ...
                'flhm_list', 'target_name'};
    for i = 1:numel(required)
        f = required{i};
        if ~isfield(pt, f) || isempty(pt.(f))
            error('PRESTUS:plantus:missingConfig', ...
                'placement.plantus.%s must be set for plantus placement.', f);
        end
    end

    script_path      = char(pt.script_path);
    env_path = char(pt.env_path);
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
    [~, plantus_hash] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', script_path));
    plantus_hash = strtrim(plantus_hash);
    [~, plantus_tag]  = system(sprintf('git -C "%s" describe --tags --exact-match 2>/dev/null', script_path));
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

    work_dir  = fullfile(parameters.path.seg, 'PlanTUS', ...
        sprintf('sub-%03d', parameters.subject_id), target_name);
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
    candidates = { ...
        fullfile(env_path, 'simnibs_env', 'bin',     'python'), ...
        fullfile(env_path, 'simnibs_env', 'Scripts', 'python.exe'), ...
        fullfile(env_path, 'python'), ...
        fullfile(env_path, 'python.exe'), ...
        fullfile(env_path, 'bin', 'python') ...
    };
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
             'Set placement.plantus.env_path to the SimNIBS installation\n' ...
             'root (e.g. /opt/SimNIBS-4.1) or to the conda env bin/ directory\n' ...
             '(e.g. ~/.conda/envs/simnibs_v4.6.0/bin).'], ...
            strjoin(candidates, '\n  '));
    end

    % Choose headless (automatic) or interactive (GUI) launcher.
    % Interactive mode opens Connectome Workbench so the user can inspect
    % the five scoring maps and manually confirm a transducer position.
    % Headless mode automatically selects the best vertex via composite score.
    %
    % placement.plantus._launcher_py (undocumented): when set, overrides the
    % launcher path directly — used by the test suite to inject mock scripts.
    if isfield(pt, '_launcher_py') && ~isempty(pt.('_launcher_py'))
        launcher = char(pt.('_launcher_py'));
        if ~isfile(launcher)
            error('PRESTUS:plantus:notFound', ...
                'Override launcher not found:\n  %s', launcher);
        end
    else
        use_gui = isfield(parameters, 'simulation') && ...
                  isfield(parameters.simulation, 'interactive') && ...
                  parameters.simulation.interactive;

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
                     '  Click a skin vertex → type "yes" to confirm placement.\n' ...
                     '  Close wb_view when done.\n']);
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

    cmd = sprintf('"%s" "%s" "%s" "%s" "%s" "%s" --plantus_root "%s" --output_dir "%s"', ...
        python_bin, launcher, filename_t1, mesh_path, mask_path, ...
        tx_config_path, script_path, work_dir);

    fprintf('PlanTUS: running placement optimisation...\n  %s\n', cmd);
    [status, cmdout] = system(cmd);

    if status ~= 0
        error('PRESTUS:plantus:failed', ...
            'PlanTUS subprocess exited with status %d.\nOutput:\n%s', status, cmdout);
    end
    fprintf('%s\n', cmdout);

    % ── Load output position matrix ───────────────────────────────────────

    % PlanTUS names its output folder after the mask file stem and places it
    % directly under the subject PlanTUS/ dir, not inside work_dir.
    plantus_root_dir = fileparts(work_dir);
    localite_mat = find_plantus_output(plantus_root_dir);
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
    if isfield(td, 'annular') && isfield(td.annular, 'curv_radius_mm') && ...
            isfield(td.annular, 'dist_geom_ep_mm')
        plane_offset = td.annular.curv_radius_mm - td.annular.dist_geom_ep_mm;
    end

    additional_offset = 0;
    if isfield(pt, 'additional_offset_mm') && ~isempty(pt.additional_offset_mm)
        additional_offset = pt.additional_offset_mm;
    end

    max_angle = 10;
    if isfield(pt, 'max_angle_deg') && ~isempty(pt.max_angle_deg)
        max_angle = pt.max_angle_deg;
    end

    connectome_wb_path = '';
    if isfield(pt, 'connectome_wb_path') && ~isempty(pt.connectome_wb_path)
        connectome_wb_path = char(pt.connectome_wb_path);
    end

    weights = struct( ...
        'skin_target_distances',     get_weight(pt, 'weight_skin_target_distances',     0.2), ...
        'skin_target_angles',        get_weight(pt, 'weight_skin_target_angles',        0.2), ...
        'skin_target_intersections', get_weight(pt, 'weight_skin_target_intersections', 0.2), ...
        'skin_skull_angles',         get_weight(pt, 'weight_skin_skull_angles',         0.2), ...
        'skull_thickness',           get_weight(pt, 'weight_skull_thickness',           0.2)  ...
    );

    % min_distance: allow up to max(flhm_list) below the closest focal depth so
    % that the min_distance filter does not exclude all vertices when focal depths
    % are tightly clustered (e.g. single-element transducer with one focal depth).
    flhm_margin = max(flhm_list);
    if isnan(flhm_margin) || flhm_margin <= 0
        flhm_margin = 10;   % 10 mm fallback when FLHM unavailable
    end
    tx_yaml = struct( ...
        'min_distance',                    max(0, min(focal_dist_list) - flhm_margin), ...
        'max_distance',                    max(focal_dist_list) + flhm_margin, ...
        'optimal_distance',                mean(focal_dist_list), ...
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
function mat_path = find_plantus_output(work_dir)
% Search work_dir recursively for the *Localite.mat written by PlanTUS.
% Errors if zero or multiple files are found.

    hits = dir(fullfile(work_dir, '**', '*Localite.mat'));

    if isempty(hits)
        error('PRESTUS:plantus:noOutput', ...
            ['PlanTUS did not produce a *Localite.mat file in:\n  %s\n' ...
             'Check the PlanTUS output above for errors.'], work_dir);
    end

    if numel(hits) > 1
        % Take the most recently modified one
        [~, idx] = max([hits.datenum]);
        hits = hits(idx);
        warning('PRESTUS:plantus:multipleOutputs', ...
            'Multiple *Localite.mat files found; using the most recent:\n  %s', ...
            fullfile(hits.folder, hits.name));
    end

    mat_path = fullfile(hits.folder, hits.name);
end
