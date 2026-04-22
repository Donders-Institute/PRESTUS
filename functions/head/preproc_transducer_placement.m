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
% See also: POSITION_TRANSDUCER_LOCALITE, TRANSDUCER_POSITIONING_START,
%           NEURONAV_SELECT_LOCALITE

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

        otherwise
            error(['Unknown placement.mode ''%s''. ' ...
                   'Use ''manual'', ''localite'', or ''heuristic''.'], mode);
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
    if ~isfield(parameters.path, 'localite') || isempty(parameters.path.localite)
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
    tpos_file = fullfile(parameters.io.output_dir, ...
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
