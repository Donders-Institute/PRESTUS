function [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = ...
    neuronav_convert_trigger_to_voxels(sub_id, positions, outputStruct, parameters)
% NEURONAV_CONVERT_TRIGGER_TO_VOXELS  Convert Localite trigger positions to voxel coordinates
%
% Extracts the mean 4×4 matrix from each requested series, derives
% transducer and focus positions via LOCALITE_MATRIX_TO_POSITIONS, and
% returns native-space RAS mm and voxel ijk coordinates. Outputs are in
% native scanner space; apply NEURONAV_CONVERT_NATIVE_TO_MNI for MNI.
%
% Use as:
%   [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = ...
%       neuronav_convert_trigger_to_voxels(sub_id, positions, outputStruct, parameters)
%
% Input:
%   sub_id       - subject identifier string (e.g. 'sub-001')
%   positions    - array of series indices to process
%   outputStruct - cell array of series structs from NEURONAV_COMPUTE_SERIES_STATISTICS
%   parameters   - PRESTUS config struct with transducer geometry and
%                  path.localite_t1 (or localite_raw) for the planning T1
%                  (see NEURONAV_FIND_PLANNING_T1)
%
% Output:
%   transducer_ras - [Nx3] transducer positions in native RAS space [mm]
%   transducer_pos - [Nx3] transducer positions in voxel space (ijk)
%   target_ras     - [Nx3] acoustic focus positions in native RAS space [mm]
%   target_pos     - [Nx3] acoustic focus positions in voxel space (ijk)
%   t1_image       - [Nx x Ny x Nz] T1 MRI volume
%
% See also: NEURONAV_COMPUTE_SERIES_STATISTICS, LOCALITE_MATRIX_TO_POSITIONS

arguments
    sub_id       {mustBeTextScalar}
    positions    {mustBeNumeric}
    outputStruct cell
    parameters   (1,1) struct
end

    % --- Load the neuronavigation planning T1 ---
    t1_path   = neuronav_find_planning_t1(sub_id, parameters);
    t1_header = niftiinfo(t1_path);
    t1_image  = niftiread(t1_path);

    % --- Pre-allocate outputs ---
    transducer_ras = zeros(max(positions), 3);
    transducer_pos = zeros(max(positions), 3);
    target_ras     = zeros(max(positions), 3);
    target_pos     = zeros(max(positions), 3);

    % --- Convert each position using shared math ---
    for i = positions
        coord_matrix = reshape(squeeze(outputStruct{i}.matrix4d_mean), [4, 4]);

        [tp, fp, tr, fr] = localite_matrix_to_positions(coord_matrix, t1_header, parameters);

        transducer_pos(i, :) = tp;
        target_pos(i, :)     = fp;
        transducer_ras(i, :) = tr;
        target_ras(i, :)     = fr;
    end

    transducer_ras = round(transducer_ras);
    target_ras     = round(target_ras);
end
