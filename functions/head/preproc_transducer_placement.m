function parameters = preproc_transducer_placement(parameters)
% PREPROC_TRANSDUCER_PLACEMENT  Resolve transducer and focus positions
%
% Dispatches transducer placement to one of three modes controlled by
% parameters.placement.mode:
%
%   'manual'    Use trans_pos / focus_pos already in parameters (default).
%               Coordinates must be in voxel space relative to the planning image.
%
%   'localite'  Parse a Localite XML file (TriggerMarkers or GUMMarkers)
%               and derive trans_pos / focus_pos from the recorded instrument
%               matrix. Requires placement.localite.file (direct path) or a
%               Localite session under path.localite (automatic selection via
%               neuronav_select_localite).
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
%               offset outward by placement.mni.skin_gap_mm (default 5 mm) using
%               the same standoff geometry as the heuristic (no search).
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
% The resolved positions are written into parameters.transducer(ti).trans_pos
% and parameters.transducer(ti).focus_pos for every configured transducer.
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
%           TRANSDUCER_POSITIONING_START, NEURONAV_SELECT_LOCALITE

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
            disp('Using manually specified transducer / focus positions.');

        % ── LOCALITE ──────────────────────────────────────────────────────
        case 'localite'
            localite_file = resolve_localite_file(parameters);
            fprintf('Reading Localite file:\n  %s\n', localite_file);

            % Load T1 header for voxel ↔ RAS coordinate transform
            filename_t1 = fullfile(parameters.path.anat, ...
                sprintf(parameters.path.t1_pattern, parameters.subject_id));
            t1_header = niftiinfo(filename_t1);

            [trans_pos, focus_pos] = ...
                position_transducer_localite(localite_file, t1_header, parameters);

            % Write resolved positions into every configured transducer
            for ti = 1:numel(parameters.transducer)
                parameters.transducer(ti).trans_pos = trans_pos;
                parameters.transducer(ti).focus_pos = focus_pos;
            end

        % ── HEURISTIC ─────────────────────────────────────────────────────
        case 'heuristic'
            parameters = run_heuristic_placement(parameters);

        % ── MNI ───────────────────────────────────────────────────────────
        case 'mni'
            parameters = run_mni_placement(parameters);

        % ── PLANTUS ───────────────────────────────────────────────────────
        case 'plantus'
            parameters = run_plantus_placement(parameters);

        otherwise
            error(['Unknown placement.mode ''%s''. ' ...
                   'Use ''manual'', ''localite'', ''heuristic'', ''mni'', or ''plantus''.'], mode);
    end

end

% ── Localite file resolver ────────────────────────────────────────────────
function localite_file = resolve_localite_file(parameters)
% Return path to Localite XML, either from placement.localite.file
% (direct path) or via neuronav_select_localite (session-based search).

    lc = parameters.placement.localite;

    if isfield(lc, 'file') && ~isempty(lc.file) && isfile(lc.file)
        localite_file = char(lc.file);
        return;
    end

    % Fall back to automatic session selection
    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'localite') || isempty(parameters.path.localite)
        error(['placement.localite.file is empty and path.localite is not set. ' ...
               'Provide one of these to use Localite placement.']);
    end

    session  = 1;
    if isfield(lc, 'session') && ~isempty(lc.session), session = lc.session; end
    markertype = 'TriggerMarkers';
    if isfield(lc, 'markertype') && ~isempty(lc.markertype), markertype = char(lc.markertype); end
    position = 1;
    if isfield(lc, 'position') && ~isempty(lc.position), position = lc.position; end

    localite_file = neuronav_select_localite( ...
        parameters.path.localite, parameters.subject_id, ...
        session, markertype, position);
end

% ── Heuristic placement dispatcher ───────────────────────────────────────
function parameters = run_heuristic_placement(parameters)
% Run the sphere-expansion heuristic (transducer_positioning_start).
% On 'matlab' platform: runs synchronously and reads the output CSV to
% update parameters. On HPC: submits a job and exits the pipeline (the
% user must re-run after the job finishes with the resolved positions).

    % Validate required heuristic config fields
    if ~isfield(parameters.placement.heuristic, 'mni_target_mm') || ...
            isempty(parameters.placement.heuristic.mni_target_mm)
        error(['placement.heuristic.mni_target_mm must be set for heuristic placement. ' ...
               'Provide a [1×3] MNI coordinate in mm.']);
    end
    if ~isfield(parameters.placement.heuristic, 'target_name') || ...
            isempty(parameters.placement.heuristic.target_name)
        error(['placement.heuristic.target_name must be set for heuristic placement. ' ...
               'Provide a label string (e.g. ''DLPFC'').']);
    end

    mni_coords  = parameters.placement.heuristic.mni_target_mm(:)';
    target_name = char(parameters.placement.heuristic.target_name);

    % Build mni_targets struct expected by transducer_positioning
    % transducer_positioning.m:69 does: target_mni = mni_targets.(target_name)
    mni_targets.(target_name) = mni_coords;

    % Build pn path-names struct from parameters.path
    pn.seg_path  = parameters.path.seg;
    pn.data_path = parameters.path.anat;
    pn.sim_path  = parameters.path.sim;

    % Expected output file from tp_select_heuristic_position (rows2vars txt)
    tpos_file = fullfile(parameters.io.dir_output, ...
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
        error(['PRESTUS:placement:hpcHeuristic', ...
               'Pipeline halted: re-run after heuristic positioning job completes.']);
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
% voxel and offset outward by skin_gap_mm using the heuristic's standoff
% geometry (transducer_scalp_geometry). No search is performed.

    % Validate required config
    if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'mni')
        error('PRESTUS:placement:missingConfig', ...
            'placement.mni sub-struct is missing. See config_default.yaml for required fields.');
    end
    mni = parameters.placement.mni;
    if ~isfield(mni, 'trans_pos_mm') || numel(mni.trans_pos_mm) ~= 3 || isempty(mni.trans_pos_mm)
        error(['placement.mni.trans_pos_mm must be a [1x3] MNI coordinate in mm ' ...
               'for mode=''mni''.']);
    end
    if ~isfield(mni, 'focus_pos_mm') || numel(mni.focus_pos_mm) ~= 3 || isempty(mni.focus_pos_mm)
        error(['placement.mni.focus_pos_mm must be a [1x3] MNI coordinate in mm ' ...
               'for mode=''mni''.']);
    end
    skin_gap_mm = 5;
    if isfield(mni, 'skin_gap_mm') && ~isempty(mni.skin_gap_mm)
        skin_gap_mm = mni.skin_gap_mm;
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

    % ── Snap transducer to scalp, then apply standoff geometry ────────────
    scalp  = tp_scalp_boundary(img);                 % [N x 3] outer-boundary voxels
    [~, j] = min(pdist2(scalp, trans_vox));          % nearest scalp voxel
    trans_pos = transducer_scalp_geometry(scalp(j,:), focus_pos, ...
                    parameters.transducer(1), pixel_size, skin_gap_mm);
    trans_pos = round(trans_pos);

    fprintf('MNI placement resolved (skin_gap_mm = %.1f):\n', skin_gap_mm);
    fprintf('  trans_pos = [%d %d %d]\n', trans_pos);
    fprintf('  focus_pos = [%d %d %d]\n', focus_pos);

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
        error('PRESTUS:placement:missingConfig', ...
            'placement.plantus sub-struct is missing. See config_default.yaml for required fields.');
    end

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
end
