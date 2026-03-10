function stats = neuronav_compute_series_statistics(localite, voxel_size, expected_segment_length, markertype)

    if nargin < 4 || isempty(markertype)
        markertype = 'TriggerMarkers';
    end
    if nargin < 3, expected_segment_length = []; end

    % Field lists
    coord_fields = {'data03Attribute', 'data13Attribute', 'data23Attribute'};
    rot_fields  = {'data00Attribute', 'data01Attribute', 'data02Attribute', ...
                   'data10Attribute', 'data11Attribute', 'data12Attribute', ...
                   'data20Attribute', 'data21Attribute', 'data22Attribute'};

    % ----------- Setup: Robust 4D Matrix Extraction ------------
    if strcmp(markertype, 'TriggerMarkers')
        if isfield(localite, 'TriggerMarker')
            triggerMarkers = [localite.TriggerMarker];
            recordingTimes = [triggerMarkers.recordingTimeAttribute];
        else
            error('TriggerMarkers expected but localite.TriggerMarker not found');
        end
    else  % GUMMarkers
        if ~isfield(localite, 'Element')
            error('GUMMarkers expected but localite.Element not found');
        end
        elements = localite.Element;
        if ~iscell(elements), elements = {elements}; end
        n_elements = numel(elements{:});
        triggerMarkers(n_elements) = struct();
        
        for i = 1:n_elements
            instr_marker = elements{:}(i).InstrumentMarker;
            matrix4d = instr_marker.Matrix4D;
            triggerMarkers(i).Matrix4D = matrix4d;
            
            % Build clean 4x4 matrix
            matrix_data = zeros(4,4);
            all_fields = fieldnames(matrix4d);
            for row = 0:3
                for col = 0:3
                    fieldname = sprintf('data%d%dAttribute', row, col);
                    if ismember(fieldname, all_fields)
                        val = matrix4d.(fieldname);
                        if isnumeric(val)
                            matrix_data(row+1, col+1) = val(1);
                        end
                    end
                end
            end
            triggerMarkers(i).Matrix4D_full = matrix_data;
            
            % Extract components
            for f = 1:3
                triggerMarkers(i).(coord_fields{f}) = matrix_data(f,4);
            end
            rot_flat = reshape(matrix_data(1:3,1:3), 1, 9);
            for f = 1:9
                triggerMarkers(i).(rot_fields{f}) = rot_flat(f);
            end
        end
        recordingTimes = (1:n_elements)';
    end

    % Process TriggerMarkers Matrix4D
    for i = 1:length(triggerMarkers)
        if isfield(triggerMarkers(i), 'Matrix4D') && ~isfield(triggerMarkers(i), 'Matrix4D_full')
            matrix4d = triggerMarkers(i).Matrix4D;
            matrix_data = zeros(4,4);
            all_fields = fieldnames(matrix4d);
            for row = 0:3
                for col = 0:3
                    fieldname = sprintf('data%d%dAttribute', row, col);
                    if ismember(fieldname, all_fields)
                        val = matrix4d.(fieldname);
                        if isnumeric(val)
                            matrix_data(row+1, col+1) = val(1);
                        end
                    end
                end
            end
            triggerMarkers(i).Matrix4D_full = matrix_data;
            for f = 1:3
                triggerMarkers(i).(coord_fields{f}) = matrix_data(f,4);
            end
            rot_flat = reshape(matrix_data(1:3,1:3), 1, 9);
            for f = 1:9
                triggerMarkers(i).(rot_fields{f}) = rot_flat(f);
            end
        end
    end

    % Sort by recording time
    [recordingTimes, sortIdx] = sort(recordingTimes);
    triggerMarkers = triggerMarkers(sortIdx);

    % Detect segments (TriggerMarkers only)
    if strcmp(markertype, 'TriggerMarkers') && length(recordingTimes) > 1
        timeDiffs = diff(recordingTimes);
        gap_threshold = median(timeDiffs) * 3;
        jumpIdx = find(timeDiffs > gap_threshold);
        segmentStarts = [1, jumpIdx + 1];
        segmentEnds = [jumpIdx, numel(recordingTimes)];
    else
        segmentStarts = 1;
        segmentEnds = length(triggerMarkers);
    end

    % Create segments
    segments = {};
    for i = 1:length(segmentStarts)
        range = segmentStarts(i):segmentEnds(i);
        if isempty(expected_segment_length) || abs(length(range) - expected_segment_length) <= 20
            segments{end+1} = triggerMarkers(range);
        end
    end

    if isempty(segments)
        warning("⚠ No valid %s segments found", markertype);
        stats = {};
        return;
    end

    % ----------- ONE CELL PER SERIES WITH Matrix4D_full AT CELL LEVEL -----------
    stats = {};  % 1×N cell array
    
    for s = 1:length(segments{:})
        stim_data = segments{1}(s);
        N = length(stim_data);
        
        % Extract all Matrix4D_full matrices for this series
        all_matrix4d = zeros(N, 4, 4);
        coords = zeros(N, 3);
        for i = 1:N
            all_matrix4d(i,:,:) = stim_data(i).Matrix4D_full;
            coords(i,:) = [stim_data(i).(coord_fields{1}), stim_data(i).(coord_fields{2}), stim_data(i).(coord_fields{3})];
        end
        
        mean_coords = mean(coords, 1);
        mean_matrix4d = mean(all_matrix4d, 1);  % [1×4×4] mean transform
        
        % Rotation submatrix for deviation calculation
        rotations = zeros(N, 9);
        for f = 1:9
            rotations(:, f) = [stim_data.(rot_fields{f})];
        end
        rotations = reshape(rotations, N, 3, 3);
        mean_rotation = mean(rotations, 1);
        
        euclideanDists = sqrt(sum((coords - mean_coords).^2, 2)) / voxel_size;
        rot_deviations = zeros(N, 1);
        for i = 1:N
            R_dev = squeeze(rotations(i,:,:)) - mean_rotation;
            rot_deviations(i) = sqrt(sum(R_dev(:).^2));
        end
        
        % CREATE COMPLETE SERIES CELL
        stats{s} = struct(...
            'series_id', s, ...
            'N_pulses', N, ...
            'voxel_size', voxel_size, ...
            'position_mean', mean_coords, ...           % [1×3]
            'matrix4d_mean', mean_matrix4d, ...         % [1×4×4] MEAN 4D MATRIX
            'matrix4d_all', all_matrix4d, ...           % [N×4×4] ALL 4D MATRICES
            'position_dev_mm', struct(...
                'mean', mean(euclideanDists), ...
                'std', std(euclideanDists), ...
                'all', euclideanDists'), ...
            'rotation_dev_rad', struct(...
                'mean', mean(rot_deviations), ...
                'std', std(rot_deviations), ...
                'all', rot_deviations'), ...
            'markers', stim_data);                      % Original enriched markers
    end
    
    fprintf('Created %d series cells (%s)\n', length(stats), markertype);
end
