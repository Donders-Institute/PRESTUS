function [results, n_series] = neuronav_compute_series_statistics(localite, parameters, voxel_size, expected_segment_length)
% NEURONAV_COMPUTE_SERIES_STATISTICS - Compute mean position and variability over stimulus trains
%
% This function processes Localite trigger marker data to analyze stimulus trains
% by segmenting the triggers into series based on time gaps and computing summary 
% statistics for each series. It calculates mean spatial coordinates per series,
% distances of each stimulus from the mean location in millimeters, and variability 
% metrics to quantify stimulus localization consistency.
%
% INPUT:
%   localite               - Struct with Localite trigger marker info, including timing and spatial coordinates
%   parameters             - Struct with transducer parameters (e.g., curvature radius)
%   voxel_size             - Voxel size of the planning image (mm)
%   expected_segment_length - (optional) Expected number of stimulations per train to filter segments
%
% OUTPUT:
%   results    - Struct containing series-wise average coordinates and distance statistics
%   n_series   - Cell array where each element is a struct array of stim data with computed distance metrics per pulse
%
% Processing steps:
%   1. Sort triggers by recording time
%   2. Detect large temporal gaps to segment triggers into stimulus trains
%   3. Filter segments by expected length if provided
%   4. Compute average coordinate per segment and Euclidean distances from the mean
%   5. Calculate mean and standard deviation of these distances per segment
%
% This allows quantifying the spatial stability and precision of stimulation localization over stimulus trains.
%
% USAGE:
%   [results, n_series] = neuronav_compute_series_statistics(localite, parameters, expected_segment_length)

    % ----------- Setup and Preprocessing ------------
    if nargin < 3
        expected_segment_length = [];  % optional filter for real use
    end

    ref_dist = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

    triggerMarkers = [localite.TriggerMarker];
    recordingTimes = [triggerMarkers.recordingTimeAttribute];

    % Sort just in case
    [recordingTimes, sortIdx] = sort(recordingTimes);
    triggerMarkers = triggerMarkers(sortIdx);

    % ----------- Detect Time Gaps to Segment Series ----------
    timeDiffs = diff(recordingTimes);
    gap_threshold = median(timeDiffs) * 3;  % typical heuristic: large gap = break
    jumpIdx = find(timeDiffs > gap_threshold);

    segmentStarts = [1, jumpIdx + 1];
    segmentEnds = [jumpIdx, numel(recordingTimes)];

    % Create segments from trigger markers
    segments = {};
    for i = 1:length(segmentStarts)
        range = segmentStarts(i):segmentEnds(i);
        if isempty(expected_segment_length) || abs(length(range) - expected_segment_length) <= 20
            segments{end+1} = triggerMarkers(range); %#ok<AGROW>
        end
    end

    if isempty(segments)
        warning("⚠ No stim series matched the expected_segment_length of %d", expected_segment_length);
    end

    % ----------- Compute Statistics for Each Series ------------
    field_names = fieldnames(triggerMarkers(1).Matrix4D);
    results = struct();
    n_series = {};

    for s = 1:length(segments)
        stim_data = [segments{s}.Matrix4D];
        series_name = sprintf('series%d', s);
        results.(series_name) = struct();

        % Mean coordinate per Matrix4D field (e.g., data03Attribute)
        for f = 1:length(field_names)
            field = field_names{f};
            results.(series_name).(field) = mean([stim_data.(field)]);
        end

        % Compute Euclidean distance to mean coordinate
        ref = results.(series_name);
        ref_voxel = [ref.data03Attribute, ref.data13Attribute, ref.data23Attribute];

        for i = 1:length(stim_data)
            pt_voxel = [stim_data(i).data03Attribute, stim_data(i).data13Attribute, stim_data(i).data23Attribute];
            dist = sqrt(sum((pt_voxel - ref_voxel).^2));
            dist_mm = dist / voxel_size;
            stim_data(i).euclideanDist = dist_mm;
        end

        % Add stat summaries to each pulse
        vals = [stim_data.euclideanDist];
        [stim_data.euclideanDist_M] = deal(mean(vals));
        [stim_data.euclideanDist_SD] = deal(std(vals));

        % Save output
        n_series{end+1} = stim_data;
    end
end
