function [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, trans_pos, focus_pos] = ...
    grid_stage(parameters)
% GRID_STAGE  Build or load the head model grid cache.
%
% Determines whether to run full grid setup or load a cached result, then
% writes a planimg sidecar for downstream stages.  When grid setup is
% disabled via modules.run_grid_setup == 0, attempts to load any existing
% cache; falls back to safe empty values if none is present.
%
% Use as:
%   [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, ...
%    trans_pos, focus_pos] = grid_stage(parameters)

arguments
    parameters (1,1) struct
end

if isfield(parameters.io, 'preproc_affix')
    preproc_affix = parameters.io.preproc_affix;
else
    preproc_affix = parameters.io.output_affix;
end

filename_grid_cache = fullfile(parameters.io.dir_cache, ...
    sprintf('sub-%03d_%s_grid_cache%s.mat', ...
    parameters.subject_id, parameters.simulation.medium, preproc_affix));

filename_planimg_sidecar = fullfile(parameters.io.dir_cache, ...
    sprintf('sub-%03d_%s_planimg%s.mat', ...
    parameters.subject_id, parameters.simulation.medium, preproc_affix));

% Settings the grid cache depends on, taken before this stage changes them.
grid_signature = cache_signature(parameters, 'grid');

if ~isfield(parameters.modules, 'run_grid_setup') || parameters.modules.run_grid_setup == 1
    internal_reuse = reuse_upstream_cache(filename_grid_cache, parameters, 'grid');
    if internal_reuse || ~confirm_overwriting(filename_grid_cache, parameters)
        if ~internal_reuse
            % User-driven reuse (io.overwrite_files: never): the cache must match
            % this run's settings. Internal sub-jobs reuse their own stage-0 cache.
            check_cache_signature(filename_grid_cache, stored_signature(filename_grid_cache), ...
                grid_signature, parameters);
        end
        disp('Loading grid cache - skipping head preprocessing.')
        [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, ...
         trans_pos, focus_pos] = load_grid_cache_file(filename_grid_cache, parameters);
    else
        parameters = load_transducer_parameters(parameters);

        [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = ...
            grid_tissue_setup(parameters);

        [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = ...
            autopad_grid_for_transducers(parameters, medium_masks, segmentation, ...
                                         bone_mask, pseudoCT, planimg);

        [parameters] = grid_transducer_location(parameters, planimg);

        [parameters, segmentation, bone_mask, pseudoCT, medium_masks] = ...
            grid_axisymmetry(parameters, segmentation, bone_mask, pseudoCT, medium_masks);

        trans_pos = parameters.transducer(1).trans_pos;
        focus_pos = parameters.transducer(1).focus_pos;

        if should_save_output(parameters.io, 'save_grid_cache')
            grid_cache_info.parameters = parameters;
            grid_cache_info.signature  = grid_signature;
            save(filename_grid_cache, ...
                'planimg', 'medium_masks', 'segmentation', 'bone_mask', 'pseudoCT', ...
                'trans_pos', 'focus_pos', 'grid_cache_info', '-v7.3');
        end
        % Sidecar always written — lets downstream stages recover planimg
        % even if the full cache is later deleted.
        save(filename_planimg_sidecar, 'planimg', '-v7.3');
        fprintf('[pipeline] planimg sidecar written: %s\n', filename_planimg_sidecar);
    end
else
    disp('No grid setup requested...no simulations will be performed.')
    if isfile(filename_grid_cache)
        fprintf('[pipeline] Grid skipped — loading grid cache: %s\n', filename_grid_cache);
        [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, ...
         trans_pos, focus_pos] = load_grid_cache_file(filename_grid_cache, parameters);
    else
        planimg      = struct();
        medium_masks = struct();
        segmentation = [];
        bone_mask    = [];
        pseudoCT     = [];
        trans_pos    = [];
        focus_pos    = [];
    end
end

% Fallback: if trans_pos was never set (e.g. no cache and grid skipped),
% pull it from the transducer struct when available.
if isempty(trans_pos)
    if isfield(parameters, 'transducer') && numel(parameters.transducer) >= 1 && ...
            isfield(parameters.transducer(1), 'trans_pos') && ~isempty(parameters.transducer(1).trans_pos)
        trans_pos = parameters.transducer(1).trans_pos;
        focus_pos = parameters.transducer(1).focus_pos;
    end
end

end

function sig = stored_signature(filename)
% Signature saved with a grid cache; [] for caches written before signatures.
    sig = [];
    info = whos('-file', filename, 'grid_cache_info');
    if isempty(info)
        return
    end
    loaded = load(filename, 'grid_cache_info');
    if isfield(loaded.grid_cache_info, 'signature')
        sig = loaded.grid_cache_info.signature;
    end
end
