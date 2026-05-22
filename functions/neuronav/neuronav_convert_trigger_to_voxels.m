function [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = ...
    neuronav_convert_trigger_to_voxels(sub_id, positions, outputStruct, parameters, pn)
% NEURONAV_CONVERT_TRIGGER_TO_VOXELS  Convert Localite trigger positions to voxel coordinates
%
% Extracts the mean 4×4 matrix from each requested series, derives
% transducer and focus positions via LOCALITE_MATRIX_TO_POSITIONS, and
% returns native-space RAS mm and voxel ijk coordinates. Outputs are in
% native scanner space; apply NEURONAV_CONVERT_NATIVE_TO_MNI for MNI.
%
% Use as:
%   [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = ...
%       neuronav_convert_trigger_to_voxels(sub_id, positions, outputStruct, parameters, pn)
%
% Input:
%   sub_id       - subject identifier string (e.g. 'sub-001')
%   positions    - array of series indices to process
%   outputStruct - cell array of series structs from NEURONAV_COMPUTE_SERIES_STATISTICS
%   parameters   - (1,1) simulation config struct with transducer geometry
%   pn           - (1,1) path names struct with data_prelocalite field
%
% Output:
%   transducer_ras - [Nx3] transducer positions in native RAS space [mm]
%   transducer_pos - [Nx3] transducer positions in voxel space (ijk)
%   target_ras     - [Nx3] acoustic focus positions in native RAS space [mm]
%   target_pos     - [Nx3] acoustic focus positions in voxel space (ijk)
%   t1_image       - [Nx x Ny x Nz] T1 MRI volume
%
% See also: NEURONAV_COMPUTE_SERIES_STATISTICS, LOCALITE_MATRIX_TO_POSITIONS

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
