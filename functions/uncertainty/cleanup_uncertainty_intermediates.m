function cleanup_uncertainty_intermediates(parameters)
% CLEANUP_UNCERTAINTY_INTERMEDIATES  Delete preprocessing cache files after uncertainty analysis
%
%   Called by uncertainty_pipeline (MATLAB platform) and by prestus_pipeline
%   (uncertainty report stage) when io.save_matrices is 0.
%
%   Files removed:
%     cache/sub-<id>_<medium>_rotated_scaled.mat
%     cache/sub-<id>_<medium>_cropped_smoothed.mat
%
%   Heating results (cache/*_heating_res*.mat) are kept: later analyses
%   read them.
%
% Use as:
%   cleanup_uncertainty_intermediates(parameters)
%
% Input:
%   parameters - PRESTUS config with subject_id and simulation.medium;
%                used by GET_OUTPUT_DIR to resolve the cache folder
%
% See also: GET_OUTPUT_DIR, UNCERTAINTY_PIPELINE

arguments
    parameters (1,1) struct
end

    medium    = parameters.simulation.medium;
    subj      = sprintf('sub-%03d', parameters.subject_id);
    cache_dir = fullfile(get_output_dir(parameters), 'cache');

    % Preprocessing cache intermediates
    delete_if_exists(fullfile(cache_dir, sprintf('%s_%s_rotated_scaled.mat', subj, medium)));
    delete_if_exists(fullfile(cache_dir, sprintf('%s_%s_cropped_smoothed.mat', subj, medium)));
end

function delete_if_exists(f)
    if isfile(f)
        delete(f);
        fprintf('[Cleanup] Removed: %s\n', f);
    end
end
