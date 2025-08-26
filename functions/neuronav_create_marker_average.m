function outputStruct = neuronav_create_marker_average(localite, results)
% CREATE_MARKER_AVERAGE - Build an averaged Localite-trigger structure for export.
%
% USAGE:
%   outputStruct = create_marker_average_struct(localite, results)
%
% INPUT:
%   localite - Original Localite struct with raw trigger markers
%   results  - Struct containing average coordinate fields for two series
%
% OUTPUT:
%   outputStruct - New Localite-compatible structure with 2 averaged trigger markers
%
% DESCRIPTION:
%   Using the position and matrix averages from `compute_series_statistics`,
%   it creates a simplified trigger structure with two static (average) matrices.

    fields = fieldnames(results.series1);
    outputStruct = localite;

    try
        outputStruct.TriggerMarker = repmat(outputStruct.TriggerMarker(1), 2, 1);
    catch
        warning('Failed to replicate trigger marker structure.');
        outputStruct = [];
        return
    end

    for f = 1:length(fields)
        field = fields{f};
        outputStruct.TriggerMarker(1).Matrix4D.(field) = results.series1.(field);
        outputStruct.TriggerMarker(2).Matrix4D.(field) = results.series2.(field);
    end
    outputStruct.TriggerMarker(1).recordingTime = '0';
    outputStruct.TriggerMarker(2).recordingTime = '1';
end
