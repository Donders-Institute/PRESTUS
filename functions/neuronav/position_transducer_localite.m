function [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
    position_transducer_localite(localite_file, t1_header, parameters)
% POSITION_TRANSDUCER_LOCALITE  Determine transducer/focus positions from a Localite XML file.
%
%   Parses a TriggerMarkers or GUMMarkers XML file using the neuronav pipeline
%   (neuronav_compute_series_statistics) and converts the resulting mean 4×4
%   transformation matrix to voxel and RAS positions via localite_matrix_to_positions.
%
%   When a TriggerMarkers file contains multiple series (e.g., repeated trigger
%   trains per position), only the first series is used. For multi-position or
%   multi-series workflows use neuronav_compute_series_statistics +
%   neuronav_convert_trigger_to_voxels directly.
%
% Input:
%   localite_file - Path to a Localite XML file (TriggerMarkers or GUMMarkers).
%                   Marker type is inferred from the filename.
%   t1_header     - NIfTI header with Transform field (from niftiinfo)
%   parameters    - Simulation config struct (transducer geometry + optional
%                   placement.localite.reference_distance_mm)
%
% Output:
%   trans_pos      - [1×3] transducer position in voxel coordinates (ijk)
%   focus_pos      - [1×3] focus position in voxel coordinates (ijk)
%   trans_pos_ras  - [1×3] transducer position in RAS space (mm)
%   focus_pos_ras  - [1×3] focus position in RAS space (mm)

    % --- Detect marker type from filename ---
    [~, fname] = fileparts(localite_file);
    if contains(fname, 'GUMMarkers', 'IgnoreCase', true)
        markertype = 'GUMMarkers';
    else
        markertype = 'TriggerMarkers';
    end

    % --- Parse using the robust neuronav XML parser ---
    localite = readstruct(localite_file);
    stats = neuronav_compute_series_statistics(localite, 1, [], markertype);

    if isempty(stats)
        error('PRESTUS:localite:noStats', ...
            'No valid marker series found in: %s', localite_file);
    end

    if length(stats) > 1
        warning('PRESTUS:localite:multipleSeries', ...
            ['%d series found in %s; using the first. ' ...
             'For multi-series workflows use neuronav_compute_series_statistics ' ...
             'and neuronav_convert_trigger_to_voxels directly.'], ...
            length(stats), fname);
    end

    % --- Extract mean 4×4 matrix from first series ---
    coord_matrix = reshape(squeeze(stats{1}.matrix4d_mean), [4, 4])';

    % --- Shared position computation ---
    [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
        localite_matrix_to_positions(coord_matrix, t1_header, parameters);

    fprintf('Localite: transducer position (vox): [%d %d %d]\n', trans_pos);
    fprintf('Localite: focus position     (vox): [%d %d %d]\n', focus_pos);
end
