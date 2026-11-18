function [kwave_medium, medium_plus, planimg, medium_masks] = ...
    medium_stage(parameters, medium_masks, planimg, pseudoCT)
% MEDIUM_STAGE  Build or load the acoustic/thermal medium property cache.
%
% Maps tissue-class labels onto acoustic and thermal material properties.
% When medium setup is disabled, loads any available cache.  If the grid
% cache has been deleted but the medium cache exists, recovers planimg from
% the compact planimg sidecar so NIfTI export and thermal stages work.
% medium_masks may be updated from the cache when the grid cache is absent.
%
% Use as:
%   [kwave_medium, medium_plus, planimg, medium_masks] = ...
%       medium_stage(parameters, medium_masks, planimg, pseudoCT)

arguments
    parameters   (1,1) struct
    medium_masks
    planimg      (1,1) struct
    pseudoCT
end

if isfield(parameters.io, 'preproc_affix')
    preproc_affix = parameters.io.preproc_affix;
else
    preproc_affix = parameters.io.output_affix;
end

% Uncertainty variants share the stage-1 head/grid (io.preproc_affix) but each
% maps its own medium, so io.medium_affix keeps their medium caches apart.
medium_affix = preproc_affix;
if isfield(parameters.io, 'medium_affix')
    medium_affix = parameters.io.medium_affix;
end

filename_medium_cache = fullfile(parameters.io.dir_cache, ...
    sprintf('sub-%03d_%s_medium_cache%s.mat', ...
    parameters.subject_id, parameters.simulation.medium, medium_affix));

filename_planimg_sidecar = fullfile(parameters.io.dir_cache, ...
    sprintf('sub-%03d_%s_planimg%s.mat', ...
    parameters.subject_id, parameters.simulation.medium, preproc_affix));

medium_signature = cache_signature(parameters, 'medium');

if ~isfield(parameters.modules, 'run_medium_setup') || parameters.modules.run_medium_setup == 1
    internal_reuse = reuse_upstream_cache(filename_medium_cache, parameters, 'medium');
    if internal_reuse || ~confirm_overwriting(filename_medium_cache, parameters)
        loaded = load(filename_medium_cache);
        if ~internal_reuse
            % User-driven reuse (io.overwrite_files: never): see GRID_STAGE.
            stored = [];
            if isfield(loaded, 'medium_cache_signature'); stored = loaded.medium_cache_signature; end
            check_cache_signature(filename_medium_cache, stored, medium_signature, parameters);
        end
        disp('Loading medium cache - skipping medium property mapping.')
        kwave_medium = loaded.kwave_medium;
        medium_plus  = loaded.medium_plus;
        if isfield(loaded, 'medium_masks')
            medium_masks = loaded.medium_masks;
        end
        % Recover planimg from sidecar when grid cache was deleted.
        if ~isfield(planimg, 'transf') || isempty(planimg)
            if isfile(filename_planimg_sidecar)
                fprintf('[pipeline] Grid cache absent — loading planimg from sidecar: %s\n', ...
                    filename_planimg_sidecar);
                s = load(filename_planimg_sidecar, 'planimg');
                planimg = s.planimg;
            else
                warning('prestus_pipeline:missingPlanimg', ...
                    ['planimg not found in workspace or sidecar (%s). ' ...
                     'NIfTI exports and thermal transforms may fail. ' ...
                     'Re-run with run_grid_setup=1 to regenerate.'], ...
                    filename_planimg_sidecar);
                planimg = struct();
            end
        end
    else
        [kwave_medium, medium_plus] = medium_setup(parameters, medium_masks, planimg, pseudoCT);

        if should_save_output(parameters.io, 'save_medium_cache')
            % medium_masks saved so this cache is self-contained for thermal/NIfTI stages.
            medium_cache_signature = medium_signature; %#ok<NASGU> saved with the cache
            save(filename_medium_cache, 'kwave_medium', 'medium_plus', 'medium_masks', ...
                'medium_cache_signature', '-v7.3');
        end
    end
else
    disp('No medium mapping requested...no simulations will be performed.')
    if isfile(filename_medium_cache)
        fprintf('[pipeline] Medium skipped — loading medium cache: %s\n', filename_medium_cache);
        loaded       = load(filename_medium_cache);
        kwave_medium = loaded.kwave_medium;
        medium_plus  = loaded.medium_plus;
        if isfield(loaded, 'medium_masks')
            medium_masks = loaded.medium_masks;
        end
    else
        kwave_medium = struct();
        medium_plus  = struct();
    end
    if ~isfield(planimg, 'transf') || isempty(planimg)
        if isfile(filename_planimg_sidecar)
            fprintf('[pipeline] Medium skipped — loading planimg from sidecar: %s\n', ...
                filename_planimg_sidecar);
            s = load(filename_planimg_sidecar, 'planimg');
            planimg = s.planimg;
        else
            planimg = struct();
        end
    end
end

end
