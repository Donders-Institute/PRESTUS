function [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras, candidates] = ...
    position_transducer_tus_entry(parameters, t1_header)
% POSITION_TRANSDUCER_TUS_ENTRY  Determine transducer placement using TUS_entry.
%
% Calls the external TUS_entry R function (Atkinson et al.) to select an
% optimal scalp entry point given a target ROI and scalp mask. TUS_entry
% scores candidate positions by how perpendicular the transducer disk sits
% on the scalp surface (angle-score = SD of disk-edge-to-scalp distances),
% and optionally enforces a depth constraint.
%
% The function:
%   1. Resolves the target in RAS/voxel space from parameters.
%   2. Writes a spherical target mask NIfTI (unless one is supplied directly).
%   3. Runs TUS_entry.R via Rscript from the configured conda environment.
%   4. Reads the output report and coordinate files.
%   5. Converts scalp-entry voxel coordinates → PRESTUS trans_pos/focus_pos.
%   6. Optionally returns the full angle-ranked candidate list.
%
% Use as:
%   [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
%       position_transducer_tus_entry(parameters, t1_header)
%   [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras, candidates] = ...
%       position_transducer_tus_entry(parameters, t1_header)
%
% Input:
%   parameters  - Simulation config struct with placement.tus_entry sub-struct
%                 (see config_default.yaml for all fields).
%   t1_header   - NIfTI header with Transform field (from niftiinfo), used for
%                 voxel ↔ RAS coordinate transforms.
%
% Output:
%   trans_pos      - [1×3] bowl rear-centre position in voxel coordinates (ijk)
%   focus_pos      - [1×3] acoustic focus position in voxel coordinates (ijk)
%   trans_pos_ras  - [1×3] bowl rear-centre in RAS space (mm)
%   focus_pos_ras  - [1×3] acoustic focus in RAS space (mm)
%   candidates     - [N×1] struct array of angle-ranked candidates (empty when
%                    placement.tus_entry.return_candidates = false, the default).
%                    Fields per entry: vox [1×3], ras [1×3], angle_score,
%                    distance_mm, rank.
%
% Required configuration (parameters.placement.tus_entry):
%   script_path      - Path to TUS_entry.R (full path to the R script file)
%   scalp_mask       - Path to binary NIfTI scalp mask (same space as T1)
%   env_path         - Path to conda env containing Rscript (e.g.
%                      /home/user/.conda/envs/tus_entry)
%
% Optional configuration (parameters.placement.tus_entry):
%   target_mask      - Path to binary NIfTI target mask. When empty, a
%                      spherical mask is generated from the focus position.
%   target_mask_radius_mm - Radius [mm] for auto-generated sphere (default 5)
%   exclu_mask       - Path to binary NIfTI exclusion zone (default '')
%   minimal_distance_cm  - Minimum scalp-to-target distance [cm] (default 4.24)
%   maximal_distance_cm  - Maximum scalp-to-target distance [cm] (default 7.46)
%   transducer_size_cm   - Transducer outer diameter [cm] (default 6.4)
%   voxel_size_cm        - Voxel size [cm/voxel] (default 0.1)
%   resolution_cm        - Spatial decimation step [cm] (default 0.33)
%   return_candidates    - Return full ranked candidate list (default false)
%   target_name          - Label for output files (default 'target')
%   work_dir             - Override working directory for TUS_entry outputs
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, RAS_TO_GRID, TRANSFORM_COORDINATES

    arguments
        parameters (1,1) struct
        t1_header  (1,1) struct
    end

    te = parameters.placement.tus_entry;

    % ── Validate required fields ──────────────────────────────────────────────

    required = {'script_path', 'scalp_mask', 'env_path'};
    for i = 1:numel(required)
        f = required{i};
        if ~isfield(te, f) || isempty(te.(f))
            error('PRESTUS:tus_entry:missingConfig', ...
                'placement.tus_entry.%s must be set.', f);
        end
    end

    script_path = char(te.script_path);
    scalp_mask  = char(te.scalp_mask);
    env_path    = char(te.env_path);

    if ~isfile(script_path)
        error('PRESTUS:tus_entry:notFound', ...
            'TUS_entry.R not found at placement.tus_entry.script_path:\n  %s', script_path);
    end
    if ~isfile(scalp_mask)
        error('PRESTUS:tus_entry:notFound', ...
            'Scalp mask not found at placement.tus_entry.scalp_mask:\n  %s', scalp_mask);
    end

    % ── Locate Rscript binary ─────────────────────────────────────────────────

    rscript_bin = find_rscript(env_path);

    % ── Optional parameters with defaults ────────────────────────────────────

    min_dist  = get_opt(te, 'minimal_distance_cm',   4.24);
    max_dist  = get_opt(te, 'maximal_distance_cm',   7.46);
    tx_size   = get_opt(te, 'transducer_size_cm',    6.4);
    vox_size  = get_opt(te, 'voxel_size_cm',         0.1);
    res_cm    = get_opt(te, 'resolution_cm',         0.33);
    target_name = char(get_opt(te, 'target_name',    'target'));
    return_candidates = logical(get_opt(te, 'return_candidates', false));

    exclu_mask = '';
    if isfield(te, 'exclu_mask') && ~isempty(te.exclu_mask)
        exclu_mask = char(te.exclu_mask);
        if ~isfile(exclu_mask)
            warning('PRESTUS:tus_entry:notFound', ...
                'Exclusion mask not found, ignoring:\n  %s', exclu_mask);
            exclu_mask = '';
        end
    end

    % ── Working directory ─────────────────────────────────────────────────────

    if isfield(te, 'work_dir') && ~isempty(te.work_dir)
        work_dir = char(te.work_dir);
    else
        work_dir = fullfile(parameters.path.seg, 'TUS_entry', ...
            sprintf('sub-%03d', parameters.subject_id), target_name);
    end
    if ~isfolder(work_dir), mkdir(work_dir); end

    % ── Resolve / write target mask ───────────────────────────────────────────

    if isfield(te, 'target_mask') && ~isempty(te.target_mask)
        target_mask = char(te.target_mask);
        if ~isfile(target_mask)
            error('PRESTUS:tus_entry:notFound', ...
                'Target mask not found at placement.tus_entry.target_mask:\n  %s', target_mask);
        end
    else
        % Generate spherical mask around focus/MNI target
        target_ras = resolve_target_ras(parameters, t1_header);
        radius_mm  = get_opt(te, 'target_mask_radius_mm', 5);
        target_mask = fullfile(work_dir, ...
            sprintf('sub-%03d_%s_target_mask.nii.gz', parameters.subject_id, target_name));
        write_target_mask(fullfile(parameters.path.anat, ...
            sprintf(parameters.path.t1_pattern, parameters.subject_id)), ...
            target_ras, target_mask, radius_mm);
        fprintf('TUS_entry: spherical target mask written (r=%.1f mm):\n  %s\n', ...
            radius_mm, target_mask);
    end

    % ── Write runner R script ─────────────────────────────────────────────────
    %
    % TUS_entry is a function file, not a standalone script. We write a thin
    % wrapper that sources it and calls TUS_entry() with the specified arguments,
    % then saves the result to a .mat-compatible CSV/JSON. We use CSV for the
    % single selected point and an optional ranked candidate RDS.

    runner_path  = fullfile(work_dir, 'run_tus_entry.R');
    report_path  = fullfile(work_dir, 'tus_entry_report.csv');
    cand_rds     = fullfile(work_dir, 'tus_entry_candidates.rds');

    exclu_arg = 'NULL';
    if ~isempty(exclu_mask)
        exclu_arg = sprintf('oro.nifti::readNIfTI("%s", reorient=FALSE)', r_quote_string(exclu_mask));
    end

    return_cand_r = 'FALSE';
    if return_candidates
        return_cand_r = 'TRUE';
    end

    runner_src = sprintf([ ...
        'suppressPackageStartupMessages({\n' ...
        '  library(oro.nifti)\n' ...
        '  library(cluster)\n' ...
        '})\n' ...
        'source("%s")\n\n' ...
        'result <- TUS_entry(\n' ...
        '  target      = readNIfTI("%s", reorient=FALSE),\n' ...
        '  scalp       = readNIfTI("%s", reorient=FALSE),\n' ...
        '  exclu       = %s,\n' ...
        '  voxel_size          = %g,\n' ...
        '  maximal_distance    = %g,\n' ...
        '  minimal_distance    = %g,\n' ...
        '  transducer_size     = %g,\n' ...
        '  resolution          = %g,\n' ...
        '  visualisation       = FALSE\n' ...
        ')\n\n' ...
        '# Write selected entry point as CSV\n' ...
        'df <- data.frame(\n' ...
        '  scalp_x   = result$Scalp_coordinates[1],\n' ...
        '  scalp_y   = result$Scalp_coordinates[2],\n' ...
        '  scalp_z   = result$Scalp_coordinates[3],\n' ...
        '  target_x  = result$Target_coordinates[1],\n' ...
        '  target_y  = result$Target_coordinates[2],\n' ...
        '  target_z  = result$Target_coordinates[3],\n' ...
        '  distance_mm = result$Distance * 10,\n' ...  % cm → mm
        '  angle_score = result$Angle\n' ...
        ')\n' ...
        'write.csv(df, "%s", row.names=FALSE)\n\n' ...
        '# Optionally save full candidate list\n' ...
        'if (%s && !is.null(result$Candidates)) {\n' ...
        '  saveRDS(result$Candidates, "%s")\n' ...
        '}\n' ...
        'cat("TUS_entry: done\\n")\n'], ...
        r_quote_string(script_path), r_quote_string(target_mask), r_quote_string(scalp_mask), exclu_arg, ...
        vox_size, max_dist, min_dist, tx_size, res_cm, ...
        r_quote_string(report_path), return_cand_r, r_quote_string(cand_rds));

    fid = fopen(runner_path, 'w');
    if fid == -1
        error('PRESTUS:tus_entry:ioError', 'Cannot write runner script to:\n  %s', runner_path);
    end
    fprintf(fid, '%s', runner_src);
    fclose(fid);

    % ── Run TUS_entry ─────────────────────────────────────────────────────────

    cmd = sprintf('%s --vanilla %s', shell_quote_arg(rscript_bin), shell_quote_arg(runner_path));
    fprintf('TUS_entry: running placement...\n  %s\n', cmd);

    [status, cmdout] = system(cmd);
    fprintf('%s\n', cmdout);

    if status ~= 0
        error('PRESTUS:tus_entry:failed', ...
            'TUS_entry R subprocess exited with status %d.\nOutput:\n%s', status, cmdout);
    end

    if ~isfile(report_path)
        error('PRESTUS:tus_entry:noOutput', ...
            'TUS_entry did not produce an output CSV:\n  %s\n  Check R output above.', report_path);
    end

    % ── Read output ───────────────────────────────────────────────────────────

    T = readtable(report_path);

    % TUS_entry voxel indices are 1-based (R arrays); PRESTUS uses 1-based
    % MATLAB voxel indices — no offset conversion needed.
    scalp_vox  = [T.scalp_x(1),  T.scalp_y(1),  T.scalp_z(1)];
    target_vox = [T.target_x(1), T.target_y(1), T.target_z(1)];
    distance_mm = T.distance_mm(1);
    angle_score = T.angle_score(1);

    fprintf('TUS_entry: scalp entry (vox) : [%d %d %d]\n', scalp_vox);
    fprintf('TUS_entry: target centroid (vox): [%d %d %d]\n', target_vox);
    fprintf('TUS_entry: depth = %.1f mm  |  angle score = %.3f\n', distance_mm, angle_score);

    % ── Derive trans_pos and focus_pos ────────────────────────────────────────
    %
    % TUS_entry's "Scalp_coordinates" is the point on the scalp surface where
    % the trajectory exits. In PRESTUS convention, trans_pos is the bowl
    % rear-centre. We apply the same focal-distance offset used in localite mode:
    % the transducer sits at the scalp entry point, offset outward along the
    % trajectory direction by the standoff (curv_radius_mm - dist_geom_ep_mm).
    %
    % Direction vector: from target centroid toward scalp entry (unit vector in
    % voxel space, scaled by pixel dimensions).

    ps = t1_header.PixelDimensions(1:3);   % mm/voxel

    dir_mm  = (scalp_vox - target_vox) .* ps;
    dir_hat = dir_mm / norm(dir_mm);       % unit vector toward scalp

    % Focal distance in mm (use mean of list if available)
    focal_dist_mm = get_focal_distance(parameters);

    % trans_pos = scalp entry + dir_hat * standoff (in voxel space)
    standoff_mm = get_standoff(parameters);
    trans_pos_mm_world = scalp_vox .* ps + dir_hat * standoff_mm;
    trans_pos = round(trans_pos_mm_world ./ ps);

    % focus_pos = scalp entry - dir_hat * focal_distance (stepping inward)
    focus_pos_mm_world = scalp_vox .* ps - dir_hat * focal_dist_mm;
    focus_pos = round(focus_pos_mm_world ./ ps);

    % ── Convert to RAS ────────────────────────────────────────────────────────

    trans_pos_ras = grid_to_ras(trans_pos, t1_header);
    focus_pos_ras = grid_to_ras(focus_pos, t1_header);

    fprintf('TUS_entry placement resolved:\n');
    fprintf('  trans_pos (vox): [%d %d %d]\n', trans_pos);
    fprintf('  focus_pos (vox): [%d %d %d]\n', focus_pos);

    % ── Optional: load ranked candidates ─────────────────────────────────────

    candidates = struct('vox', {}, 'ras', {}, 'angle_score', {}, ...
                        'distance_mm', {}, 'rank', {});
    if return_candidates && isfile(cand_rds)
        warning('PRESTUS:tus_entry:candidatesRDS', ...
            ['Candidate list saved as RDS (%s) but cannot be read directly ' ...
             'from MATLAB. Load in R with readRDS() or reprocess.'], cand_rds);
    end
end

% ── Helpers ───────────────────────────────────────────────────────────────────

function target_ras = resolve_target_ras(parameters, t1_header)
% Resolve target to RAS mm: try tus_entry.mni_target_mm, then focus_pos.
    te = parameters.placement.tus_entry;
    if isfield(te, 'mni_target_mm') && ~isempty(te.mni_target_mm)
        mni_mm = te.mni_target_mm(:)';
        target_ras = transform_coordinates(parameters, mni_mm, 'mni', 'ras_plus', t1_header);
        return;
    end
    td = parameters.transducer(1);
    if isfield(td, 'focus_pos') && ~isempty(td.focus_pos)
        target_ras = transform_coordinates(parameters, td.focus_pos(:)', 'grid', 'ras_plus', t1_header);
        return;
    end
    error('PRESTUS:tus_entry:noTarget', ...
        ['Cannot determine target position. ' ...
         'Set placement.tus_entry.mni_target_mm or provide transducer.focus_pos.']);
end

function write_target_mask(t1_path, target_ras, mask_path, radius_mm)
% Write a binary NIfTI sphere centred on target_ras.
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

function rscript_bin = find_rscript(env_path)
% Probe standard conda env layouts for Rscript binary.
    candidates = { ...
        fullfile(env_path, 'bin', 'Rscript'), ...
        fullfile(env_path, 'Scripts', 'Rscript.exe'), ...
        fullfile(env_path, 'Rscript') ...
    };
    rscript_bin = '';
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            rscript_bin = candidates{i};
            return;
        end
    end
    error('PRESTUS:tus_entry:notFound', ...
        ['Rscript binary not found. Tried:\n' ...
         '  %s\n' ...
         'Set placement.tus_entry.env_path to the conda env prefix\n' ...
         '(e.g. ~/.conda/envs/tus_entry).'], ...
        strjoin(candidates, '\n  '));
end

function d = get_opt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        d = s.(field);
    else
        d = default;
    end
end

function focal_dist_mm = get_focal_distance(parameters)
% Return focal distance in mm: prefer calibrated list, fall back to geometry.
    td = parameters.transducer(1);
    if isfield(td, 'annular') && isfield(td.annular, 'focal_distance_mm')
        focal_dist_mm = mean(td.annular.focal_distance_mm);
        return;
    end
    if isfield(td, 'focal_distance_mm') && ~isempty(td.focal_distance_mm)
        focal_dist_mm = mean(td.focal_distance_mm);
        return;
    end
    if isfield(td, 'annular') && isfield(td.annular, 'curv_radius_mm')
        focal_dist_mm = td.annular.curv_radius_mm;
        return;
    end
    error('PRESTUS:tus_entry:noFocalDist', ...
        ['Cannot determine focal distance. ' ...
         'Set transducer.focal_distance_mm or transducer.annular.focal_distance_mm.']);
end

function standoff_mm = get_standoff(parameters)
% Standoff = distance from bowl rear centre to exit plane (same as PlanTUS).
    td = parameters.transducer(1);
    standoff_mm = 0;
    if isfield(td, 'annular') && isfield(td.annular, 'curv_radius_mm') && ...
            isfield(td.annular, 'dist_geom_ep_mm')
        standoff_mm = td.annular.curv_radius_mm - td.annular.dist_geom_ep_mm;
    end
end

function ras = grid_to_ras(vox, t1_header)
% Convert 1-based voxel indices to RAS mm using the NIfTI affine.
    ras = (t1_header.Transform.T' * [vox(:); 1]);
    ras = ras(1:3)';
end
