function [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
    position_transducer_localite(localite_file, t1_header, parameters)
% POSITION_TRANSDUCER_LOCALITE  Determine transducer/focus positions from a Localite XML file.
%
% Parses a TriggerMarkers, GUMMarkers or InstrumentMarker XML file (type
% detected from its content), picks one series and converts its mean 4×4
% instrument matrix to voxel and RAS positions via LOCALITE_MATRIX_TO_POSITIONS.
% Series selection follows NEURONAV_LOCALITE_POSE: placement.localite.target
% (by description) if set, else placement.localite.position (default 1).
%
% Use as:
%   [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
%       position_transducer_localite(localite_file, t1_header, parameters)
%
% Input:
%   localite_file - Path to a Localite XML file.
%   t1_header     - NIfTI header with Transform field (from niftiinfo) of the
%                   planning T1 the Localite coordinates refer to
%   parameters    - Simulation config struct (transducer geometry, optional
%                   placement.localite.position / target / reference_distance_mm)
%
% Output:
%   trans_pos      - [1×3] transducer position in voxel coordinates (ijk)
%   focus_pos      - [1×3] focus position in voxel coordinates (ijk)
%   trans_pos_ras  - [1×3] transducer position in RAS space (mm)
%   focus_pos_ras  - [1×3] focus position in RAS space (mm)
%
% See also: NEURONAV_LOCALITE_POSE, LOCALITE_MATRIX_TO_POSITIONS

    arguments
        localite_file (1,:) char
        t1_header     (1,1) struct
        parameters    (1,1) struct
    end

    parameters.placement.localite.file = localite_file;
    coord_matrix = neuronav_localite_pose(parameters);

    [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
        localite_matrix_to_positions(coord_matrix, t1_header, parameters);

    fprintf('Localite: transducer position (vox): [%d %d %d]\n', trans_pos);
    fprintf('Localite: focus position     (vox): [%d %d %d]\n', focus_pos);
end
