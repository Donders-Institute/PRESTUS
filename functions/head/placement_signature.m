function sig = placement_signature(parameters, mode)
% PLACEMENT_SIGNATURE  Settings a placement run's result depends on
%
% Stored in the run record (placement_run.mat) and compared on reuse with
% CHECK_CACHE_SIGNATURE, like CACHE_SIGNATURE for the grid/medium caches:
% the mode's own settings (placement.<mode>, without the work folder, tool
% install paths and derived outputs) and the transducers without their
% positions (the result) and drive level; for modes that read the SimNIBS
% segmentation, its content (SEGMENTATION_FINGERPRINT).
%
% Use as:
%   sig = placement_signature(parameters, mode)
%
% Input:
%   parameters - PRESTUS parameter struct (before placement runs)
%   mode       - placement mode with a run folder (see PLACEMENT_RUN_DIR)
%
% Output:
%   sig - struct of the relevant settings
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, CACHE_SIGNATURE, CHECK_CACHE_SIGNATURE

arguments
    parameters (1,1) struct
    mode       {mustBeTextScalar}
end

    mode = char(mode);
    sig = struct();
    sig.signature_version = 2;      % 2: exit plane from dist_geom_ep_mm, placement.skin_gap_mm
    sig.mode = mode;
    sig.subject_id = parameters.subject_id;
    s = parameters.placement.(mode);
    for f = {'work_dir', 'target_ras', 'script_path', 'env_path'}
        if isfield(s, f{1}); s = rmfield(s, f{1}); end
    end
    % PlanTUS: interactive [] follows simulation.interactive, as before the
    % key existed, so records written then still match
    if isfield(s, 'interactive') && isempty(s.interactive)
        s = rmfield(s, 'interactive');
    end
    if isfield(s, 'vertex') && isempty(s.vertex)    % likewise: no vertex chosen
        s = rmfield(s, 'vertex');
    end
    if isfield(s, 'candidate') && isempty(s.candidate)  % heuristic: no candidate chosen
        s = rmfield(s, 'candidate');
    end
    sig.settings = s;
    if isfield(parameters.placement, 'skin_gap_mm')
        sig.skin_gap_mm = parameters.placement.skin_gap_mm;
    end
    t = parameters.transducer;
    for f = {'trans_pos', 'focus_pos', 'trans_pos_ras', 'focus_pos_ras', ...
             'target_isppa_wcm2', 'tx_calibration'}
        if isfield(t, f{1}); t = rmfield(t, f{1}); end
    end
    sig.transducer = t;
    % Modes that read the SimNIBS segmentation (scalp, mesh, MNI warps): a re-run
    % of charm changes the placement (field only when there is a segmentation)
    if any(strcmp(mode, {'heuristic', 'mni', 'plantus', 'tus_entry'}))
        fp = segmentation_fingerprint(parameters);
        if ~isempty(fp)
            sig.segmentation = fp;
        end
    end
end
