function cleanup_uncertainty_intermediates(parameters, affixes)
% CLEANUP_UNCERTAINTY_INTERMEDIATES  Delete large intermediate cache files after uncertainty analysis
%
%   Called automatically by uncertainty_pipeline (MATLAB platform) and by
%   prestus_pipeline (HPC platform, after the report job runs).
%
% Use as:
%   cleanup_uncertainty_intermediates(parameters, affixes)
%
% Input:
%   parameters - PRESTUS config with subject_id and simulation.medium;
%                used by GET_OUTPUT_DIR to resolve the cache folder
%   affixes    - struct with fields .default, .liberal, .conservative
%                (suffix strings for each uncertainty variant)
%
% See also: GET_OUTPUT_DIR, UNCERTAINTY_PIPELINE

arguments
    parameters (1,1) struct
    affixes    (1,1) struct
end
%
% Files removed (subject to io flags):
%   cache/<subj>_<medium>_heating_res<affix>.mat  — removed unless save_thermal_matrices=1
%   cache/<subj>_<medium>_after_rotating_and_scaling.mat
%   cache/<subj>_<medium>_after_cropping_and_smoothing.mat
%
% Files always kept:
%   output_table<affix>.csv, *_report<affix>.html, *.png,
%   parameters_<affix>.mat, uncertainty_report.html

    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = get_output_dir(parameters);
    subj       = sprintf('sub-%03d', subject_id);

    cache_dir = fullfile(output_dir, 'cache');

    % Heating result matrices for all three variants
    if ~should_save_output(parameters.io, 'save_thermal_matrices')
        variant_affixes = {affixes.default, affixes.liberal, affixes.conservative};
        for i = 1:numel(variant_affixes)
            af = variant_affixes{i};
            f  = fullfile(cache_dir, sprintf('%s_%s%s_heating_res.mat', subj, medium, af));
            delete_if_exists(f);
        end
    end

    % Preprocessing cache intermediates
    delete_if_exists(fullfile(cache_dir, sprintf('%s_%s_after_rotating_and_scaling.mat', subj, medium)));
    delete_if_exists(fullfile(cache_dir, sprintf('%s_%s_after_cropping_and_smoothing.mat', subj, medium)));
end

function delete_if_exists(f)
    if isfile(f)
        delete(f);
        fprintf('[Cleanup] Removed: %s\n', f);
    end
end
