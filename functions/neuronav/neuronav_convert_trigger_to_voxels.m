function [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = ...
    neuronav_convert_trigger_to_voxels(sub_id, positions, outputStruct, parameters, pn)
% NEURONAV_CONVERT_TRIGGER_TO_VOXELS - Convert Localite trigger positions to voxel coordinates.
%
% USAGE:
%   [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = ...
%       neuronav_convert_trigger_to_voxels(sub_id, positions, outputStruct, parameters, pn)
%
% INPUT:
%   sub_id       - Subject identifier string (e.g., 'sub-001')
%   positions    - Array of stimulation position indices (e.g., 1:2 for two sides)
%   outputStruct - Cell array of series structs from neuronav_compute_series_statistics;
%                  each cell contains a matrix4d_mean field ([1×4×4])
%   parameters   - Simulation config with transducer geometry
%   pn           - Path struct; pn.data_prelocalite must point to the folder
%                  containing the subject T1 NIfTI
%
% OUTPUT:
%   transducer_ras - [N×3] transducer positions in RAS space (mm), rounded
%   transducer_pos - [N×3] transducer positions in voxel space (ijk)
%   target_ras     - [N×3] acoustic focus positions in RAS space (mm), rounded
%   target_pos     - [N×3] acoustic focus positions in voxel space (ijk)
%   t1_image       - Loaded T1 MRI volume (3D matrix)
%
% COORDINATE SYSTEM NOTE:
%   All outputs are in native scanner/MRI space (RAS mm or voxel ijk). They are
%   NOT in MNI space. For standard-space coordinates apply a normalisation step
%   post hoc (e.g., via neuronav_convert_native_to_MNI).

    % --- Load T1 MRI ---
    t1_file   = dir(fullfile(pn.data_prelocalite, sprintf('%s_*T1*.nii*', sub_id)));
    t1_header = niftiinfo(fullfile(t1_file(1).folder, t1_file(1).name));
    t1_image  = niftiread(fullfile(t1_file(1).folder, t1_file(1).name));

    % --- Pre-allocate outputs ---
    transducer_ras = zeros(max(positions), 3);
    transducer_pos = zeros(max(positions), 3);
    target_ras     = zeros(max(positions), 3);
    target_pos     = zeros(max(positions), 3);

    % --- Convert each position using shared math ---
    for i = positions
        coord_matrix = reshape(squeeze(outputStruct{i}.matrix4d_mean), [4, 4])';

        [tp, fp, tr, fr] = localite_matrix_to_positions(coord_matrix, t1_header, parameters);

        transducer_pos(i, :) = tp;
        target_pos(i, :)     = fp;
        transducer_ras(i, :) = tr;
        target_ras(i, :)     = fr;
    end

    transducer_ras = round(transducer_ras);
    target_ras     = round(target_ras);
end
