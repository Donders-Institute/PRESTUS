function [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, ...
          trans_pos, focus_pos] = load_grid_cache_file(filename, parameters)
% LOAD_GRID_CACHE_FILE  Load a PRESTUS grid cache and restore pipeline state.
%
% Encapsulates the shared cache-load sequence used in two places in
% prestus_pipeline:
%   1. run_grid_setup == 1, cache-hit branch (skip head preprocessing).
%   2. run_grid_setup == 0, optional cache-load for downstream stages.
%
% Both branches need identical logic:
%   - load the .mat file into a local scope (not the caller's workspace)
%   - restore the grid keys the grid stage computes (dims, default_dims,
%     resolution_mm, display_orientation) from the cache; other grid keys
%     keep the caller's settings
%   - conditionally restore parameters.transducer: keep the caller's array
%     when the cache was built with MORE transducers (async shared preproc),
%     otherwise use the cached array.
%   - recompute derived transducer fields via load_transducer_parameters
%   - refresh trans_pos / focus_pos: for async per-transducer jobs
%     (io.async_transducer_idx present) the grid-space positions are copied
%     from the correct cached transducer entry (set by grid_transducer_location
%     during Stage 0) so that downstream code (acoustic_grid_fov, source_create,
%     acoustic_analysis) uses PRESTUS-grid voxel coordinates, not the input-space
%     coordinates (e.g. T1 voxels) that the posthoc script injected.
%
% Use as:
%   [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, ...
%    trans_pos, focus_pos] = load_grid_cache_file(filename, parameters)
%
% Input:
%   filename   - path to the grid cache .mat file (must exist)
%   parameters - current pipeline parameters struct (parameters.transducer
%                reflects the caller's per-transducer selection)
%
% Output:
%   parameters   - updated: grid-stage keys of parameters.grid from cache;
%                  parameters.transducer preserved or restored (see above)
%   planimg      - planning image loaded from cache
%   medium_masks - tissue mask struct loaded from cache
%   segmentation - segmentation volume loaded from cache
%   bone_mask    - bone mask loaded from cache
%   pseudoCT     - pseudo-CT volume loaded from cache
%   trans_pos    - parameters.transducer(1).trans_pos (caller's transducer)
%   focus_pos    - parameters.transducer(1).focus_pos (caller's transducer)

transducer_before_load = parameters.transducer;

cached = load(filename);

% The cache supplies what the grid stage computed; every other grid key keeps
% the caller's value. A run that reuses another run's grid (io.preproc_affix)
% may differ in settings the grid stage does not use, e.g. a partial acoustic
% FOV (grid.acoustic_fov_*) or none for a job reading a back-projected field.
grid_settings = parameters.grid;
parameters.grid = cached.grid_cache_info.parameters.grid;
keep = setdiff(fieldnames(grid_settings), grid_stage_outputs());
for k = 1:numel(keep)
    parameters.grid.(keep{k}) = grid_settings.(keep{k});
end

% Preserve the caller's transducer array when the cache was built with more
% transducers (async shared preproc). Restoring the full array would
% overwrite the per-transducer selection and cause every per-transducer job
% to simulate transducer(1).
%
% However, the caller's transducer may carry positions in the original input
% space (e.g. T1 voxels from a posthoc Localite script) rather than the
% PRESTUS grid space set by grid_transducer_location during the preproc run.
% When parameters.io.async_transducer_idx is set, copy the grid-space
% trans_pos / focus_pos from the correct cached transducer entry so that
% all downstream geometry (acoustic_grid_fov, source_create, …) uses the
% right coordinate system.
if numel(cached.grid_cache_info.parameters.transducer) > numel(transducer_before_load)
    parameters.transducer = transducer_before_load;
    if isfield(parameters, 'io') && isfield(parameters.io, 'async_transducer_idx')
        ti = parameters.io.async_transducer_idx;
        cached_transducers = cached.grid_cache_info.parameters.transducer;
        if ti >= 1 && ti <= numel(cached_transducers)
            cached_tr = cached_transducers(ti);
            if isfield(cached_tr, 'trans_pos') && ~isempty(cached_tr.trans_pos)
                parameters.transducer(1).trans_pos = cached_tr.trans_pos;
            end
            if isfield(cached_tr, 'focus_pos') && ~isempty(cached_tr.focus_pos)
                parameters.transducer(1).focus_pos = cached_tr.focus_pos;
            end
        else
            warning('load_grid_cache_file:badTransducerIdx', ...
                ['async_transducer_idx (%d) is out of range for cached transducer ' ...
                 'array (length %d). Grid positions may be wrong.'], ...
                ti, numel(cached_transducers));
        end
    else
        warning('load_grid_cache_file:missingTransducerIdx', ...
            ['Cache has more transducers than caller but parameters.io.async_transducer_idx ' ...
             'is not set. Grid positions may be in the wrong coordinate space. ' ...
             'Ensure make_acoustic_params sets io.async_transducer_idx.']);
    end
elseif numel(cached.grid_cache_info.parameters.transducer) == numel(transducer_before_load)
    % Same transducers: keep the caller's settings (drive, phases, ISPPA
    % targets, ...) and take only what the grid stage computed, the grid-space
    % positions (GRID_TRANSDUCER_LOCATION). Restoring the whole cached array
    % would replace e.g. the targets of a run that reuses another run's grid.
    cached_transducers = cached.grid_cache_info.parameters.transducer;
    parameters.transducer = transducer_before_load;
    for ti = 1:numel(cached_transducers)
        parameters.transducer(ti).trans_pos = cached_transducers(ti).trans_pos;
        parameters.transducer(ti).focus_pos = cached_transducers(ti).focus_pos;
    end
else
    parameters.transducer = cached.grid_cache_info.parameters.transducer;
end

% Recompute derived transducer fields (handles caches from older PRESTUS
% versions that did not store all derived fields).
parameters = load_transducer_parameters(parameters);

% Refresh locals from parameters.transducer(1) so they reflect the
% caller's per-transducer selection, not whatever the cache stored.
trans_pos = parameters.transducer(1).trans_pos;
focus_pos = parameters.transducer(1).focus_pos;

% Pass through workspace variables saved by grid_tissue_setup / grid_transducer_location.
planimg      = cached.planimg;
medium_masks = cached.medium_masks;
segmentation = cached.segmentation;
bone_mask    = cached.bone_mask;
pseudoCT     = cached.pseudoCT;

% Backfill display_orientation for caches built before nifti_slice_orientation
% was wired into grid_tissue_setup.  Without it every plot_overlay call falls
% back to no-orient defaults and x/z-slice plots appear wrongly rotated.
if isfield(parameters, 'grid') && ~isfield(parameters.grid, 'display_orientation') && ...
        isfield(parameters.grid, 'mode') && strcmp(parameters.grid.mode, 'ras_plus') && ...
        isfield(planimg, 't1_header') && ~isempty(planimg.t1_header)
    ori = struct();
    for ax_lgcf = {'x', 'y', 'z'}
        [dt, fr, fc] = nifti_slice_orientation(planimg.t1_header, ax_lgcf{1});
        ori.(ax_lgcf{1}) = struct('do_transpose', dt, 'flip_rows', fr, 'flip_cols', fc);
    end
    parameters.grid.display_orientation = ori;
end

end

function names = grid_stage_outputs()
% Grid keys the grid stage sets (GRID_TISSUE_SETUP, PREPROC_HEAD,
% AUTOPAD_GRID_FOR_TRANSDUCERS): restored from the cache, not from the caller.
names = {'dims', 'default_dims', 'resolution_mm', 'display_orientation'};
end
