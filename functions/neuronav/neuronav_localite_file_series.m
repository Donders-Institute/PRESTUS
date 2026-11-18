function [series, markertype] = neuronav_localite_file_series(file, seg_len)
% NEURONAV_LOCALITE_FILE_SERIES  The marker series of one Localite XML export.
%
% Reads a TriggerMarkers, GUMMarkers or InstrumentMarker file (detected from
% its content) and returns one entry per detected series, with its mean 4×4
% instrument matrix. Repeated triggers of a series are averaged by
% NEURONAV_COMPUTE_SERIES_STATISTICS. Placement from a single file
% (placement.localite.file, NEURONAV_LOCALITE_POSE) and the Nexus import
% (NEXUS_NEURONAV 'file_series') both use this list.
%
% Use as:
%   [series, markertype] = neuronav_localite_file_series(file, seg_len)
%
% Input:
%   file    - path of the Localite XML export
%   seg_len - expected number of pulses per series (TriggerMarkers only;
%             [] = no filter), placement.localite.expected_segment_length
%
% Output:
%   series     - struct array, one per series, with fields
%                .index        series number (placement.localite.position)
%                .description  marker description ('' if none;
%                              placement.localite.target)
%                .n_markers    number of markers averaged
%                .pos_dev_mm   mean distance of the markers from their mean [mm]
%                .is_pose      false for a planned EntryTargetPair (a point)
%                .tracked      false when the mean position is at the origin
%                .matrix       [4×4] mean instrument matrix (column 4 =
%                              position, column 1 = beam direction)
%   markertype - 'TriggerMarkers', 'GUMMarkers' or 'InstrumentMarker'
%
% See also: NEURONAV_LOCALITE_POSE, NEURONAV_COMPUTE_SERIES_STATISTICS,
%           LOCALITE_MATRIX_TO_POSITIONS

    localite = readstruct(file);
    if isfield(localite, 'TriggerMarker')
        markertype = 'TriggerMarkers';
    elseif isfield(localite, 'Element')
        markertype = 'GUMMarkers';
    elseif isfield(localite, 'InstrumentMarker')
        markertype = 'InstrumentMarker';
    else
        error('PRESTUS:localite:unknownFormat', ...
            'Not a TriggerMarkers, GUMMarkers or InstrumentMarker file: %s', file);
    end
    if ~strcmp(markertype, 'TriggerMarkers'), seg_len = []; end

    stats = neuronav_compute_series_statistics(localite, 1, seg_len, markertype);
    if isempty(stats)
        error('PRESTUS:localite:noStats', 'No valid marker series found in: %s', file);
    end

    series = repmat(struct('index', 0, 'description', '', 'n_markers', 0, ...
        'pos_dev_mm', NaN, 'is_pose', true, 'tracked', true, 'matrix', zeros(4)), ...
        1, numel(stats));
    for k = 1:numel(stats)
        r  = stats{k};
        mk = r.markers(1);
        M  = reshape(squeeze(r.matrix4d_mean), [4, 4]);
        series(k).index       = k;
        series(k).description = series_description(mk);
        series(k).n_markers   = r.N_pulses;
        series(k).pos_dev_mm  = r.position_dev_mm.mean;
        series(k).is_pose     = ~(isfield(mk, 'element_type') && strcmp(mk.element_type, 'EntryTargetPair'));
        series(k).tracked     = ~all(abs(M(1:3, 4)) < 1e-9);
        series(k).matrix      = M;
    end
end


function d = series_description(mk)
    d = '';
    if isfield(mk, 'description') && ~isempty(mk.description) && ~ismissing(string(mk.description))
        d = strtrim(char(string(mk.description)));
    end
end
