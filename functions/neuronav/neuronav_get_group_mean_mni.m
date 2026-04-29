function [mean_mni_trans_mm, mean_mni_targ_mm, mean_mni_trans_vox, mean_mni_targ_vox] = ...
    neuronav_get_group_mean_mni(session_target, parameters)
% NEURONAV_GET_GROUP_MEAN_MNI  Compute group-level average MNI coordinates
%
% Reads per-subject session CSVs from parameters.path.sim, filters to
% non-interpolated rows, and returns mean transducer and target positions
% for the specified target label in both RAS mm and voxel MNI coordinates.
%
% Use as:
%   [mean_mni_trans_mm, mean_mni_targ_mm, mean_mni_trans_vox, mean_mni_targ_vox] = ...
%       neuronav_get_group_mean_mni(session_target, parameters)
%
% Input:
%   session_target - target label string to average over
%   parameters     - (1,1) simulation parameters struct with path.sim field
%
% Output:
%   mean_mni_trans_mm  - [1x3] mean transducer position in MNI RAS space [mm]
%   mean_mni_targ_mm   - [1x3] mean target position in MNI RAS space [mm]
%   mean_mni_trans_vox - [1x3] mean transducer position in MNI voxel space
%   mean_mni_targ_vox  - [1x3] mean target position in MNI voxel space
%
% See also: NEURONAV_EXPORT_SESSION_CSV, NEURONAV_CONVERT_NATIVE_TO_MNI

% Load all subject session CSVs
csv_files = dir(fullfile(parameters.path.sim, 'sub-*', '*_info.csv'));

expected_vars = { ...
    'trans_MNI_mm_x', 'trans_MNI_mm_y', 'trans_MNI_mm_z', ...
    'target_MNI_mm_x', 'target_MNI_mm_y', 'target_MNI_mm_z', ...
    'trans_MNI_vox_x', 'trans_MNI_vox_y', 'trans_MNI_vox_z', ...
    'target_MNI_vox_x', 'target_MNI_vox_y', 'target_MNI_vox_z', ...
    'target_label', 'interpolated'};

all_data = table();

for k = 1:numel(csv_files)
    try
        T = readtable(fullfile(csv_files(k).folder, csv_files(k).name));
        if all(ismember(expected_vars, T.Properties.VariableNames))
            all_data = [all_data; T(:, expected_vars)];
        else
            warn('Skipping file with missing fields: %s', csv_files(k).name);
        end
    catch
        warn('Could not read: %s', csv_files(k).name);
    end
end

% Initialize outputs
n_targets = numel(session_target);
mean_mni_trans_mm  = nan(n_targets, 3);
mean_mni_targ_mm   = nan(n_targets, 3);
mean_mni_trans_vox = nan(n_targets, 3);
mean_mni_targ_vox  = nan(n_targets, 3);

% Loop over each target and average mm + voxel values
for i = 1:n_targets
    name = strtrim(lower(session_target{i}));
    idx = strcmpi(strtrim(all_data.target_label), name) & all_data.interpolated == 0;

    if any(idx)
        % RAS mm
        mean_mni_trans_mm(i,:) = mean([ ...
            all_data.trans_MNI_mm_x(idx), ...
            all_data.trans_MNI_mm_y(idx), ...
            all_data.trans_MNI_mm_z(idx) ], 1, 'omitnan');

        mean_mni_targ_mm(i,:) = mean([ ...
            all_data.target_MNI_mm_x(idx), ...
            all_data.target_MNI_mm_y(idx), ...
            all_data.target_MNI_mm_z(idx) ], 1, 'omitnan');

        % Voxel indices
        mean_mni_trans_vox(i,:) = mean([ ...
            all_data.trans_MNI_vox_x(idx), ...
            all_data.trans_MNI_vox_y(idx), ...
            all_data.trans_MNI_vox_z(idx) ], 1, 'omitnan');

        mean_mni_targ_vox(i,:) = mean([ ...
            all_data.target_MNI_vox_x(idx), ...
            all_data.target_MNI_vox_y(idx), ...
            all_data.target_MNI_vox_z(idx) ], 1, 'omitnan');
    else
        warn('No data found for target %s with interpolated == 0', name);
    end
end

% Optional rounding for mm coordinates
mean_mni_trans_mm  = round(mean_mni_trans_mm);
mean_mni_targ_mm   = round(mean_mni_targ_mm);
mean_mni_trans_vox = round(mean_mni_trans_vox);
mean_mni_targ_vox  = round(mean_mni_targ_vox);

end
