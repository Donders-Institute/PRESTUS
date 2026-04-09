function cleanup_uncertainty_intermediates(parameters, affixes)
% CLEANUP_UNCERTAINTY_INTERMEDIATES  Delete large intermediate files after
% the uncertainty report (stage 5) when save_matrices = 0.
%
% Called automatically by uncertainty_pipeline (MATLAB platform) and by
% prestus_pipeline (HPC platform, after the report job runs).
%
% Files removed:
%   heating_res<affix>.mat  — thermal matrices for each variant (large)
%   debug/<subj>_after_rotating_and_scaling.mat
%   debug/<subj>_<medium>_after_cropping_and_smoothing.mat
%
% Files always kept:
%   output_table<affix>.csv, *_report<affix>.html, *.png,
%   parameters_<affix>.mat, uncertainty_report.html

    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = get_output_dir(parameters);
    subj       = sprintf('sub-%03d', subject_id);

    % Heating result matrices for all three variants
    variant_affixes = {affixes.default, affixes.liberal, affixes.conservative};
    for i = 1:numel(variant_affixes)
        af = variant_affixes{i};
        f  = fullfile(output_dir, sprintf('%s_%s_heating_res%s.mat', subj, medium, af));
        delete_if_exists(f);
    end

    % Debug intermediates
    debug_dir = fullfile(output_dir, 'debug');
    delete_if_exists(fullfile(debug_dir, sprintf('%s_after_rotating_and_scaling.mat', subj)));
    delete_if_exists(fullfile(debug_dir, sprintf('%s_%s_after_cropping_and_smoothing.mat', subj, medium)));
end

function delete_if_exists(f)
    if isfile(f)
        delete(f);
        fprintf('[Cleanup] Removed: %s\n', f);
    end
end
