function stats = neuronav_compute_series_statistics(localite, voxel_size, expected_segment_length, markertype)
% NEURONAV_COMPUTE_SERIES_STATISTICS  Compute per-series position and orientation statistics
%
% Segments a parsed Localite XML struct into series by timing gaps and
% computes mean position, mean 4×4 matrix, and position/rotation
% deviations for each series.
%
% Use as:
%   stats = neuronav_compute_series_statistics(localite, voxel_size)
%   stats = neuronav_compute_series_statistics(localite, voxel_size, expected_segment_length, markertype)
%
% Input:
%   localite                - struct from neuronav_select_localite or readstruct
%   voxel_size              - voxel size for converting distances to voxels [mm]
%   expected_segment_length - expected pulses per series ([] = accept all)
%   markertype              - 'TriggerMarkers' (default), 'GUMMarkers', or 'InstrumentMarker'
%
% Output:
%   stats - 1×N cell array; each cell is a struct with fields:
%           series_id, N_pulses, voxel_size, position_mean [1×3],
%           matrix4d_mean [1×4×4], matrix4d_all [N×4×4],
%           position_dev_mm, rotation_dev_rad, markers
%
% See also: POSITION_TRANSDUCER_LOCALITE, NEURONAV_CONVERT_TRIGGER_TO_VOXELS

    if nargin < 4 || isempty(markertype)
        markertype = 'TriggerMarkers';
    end
    if nargin < 3, expected_segment_length = []; end

    % Field lists (used for TriggerMarkers / InstrumentMarker Matrix4D paths)
    coord_fields = {'data03Attribute', 'data13Attribute', 'data23Attribute'};
    rot_fields  = {'data00Attribute', 'data01Attribute', 'data02Attribute', ...
                   'data10Attribute', 'data11Attribute', 'data12Attribute', ...
                   'data20Attribute', 'data21Attribute', 'data22Attribute'};

    % ----------- Setup: extract markers array and recording times -----------

    if strcmp(markertype, 'TriggerMarkers')
        % ---- TriggerMarkers: one marker per pulse, full 4×4 matrix --------
        if ~isfield(localite, 'TriggerMarker')
            error('TriggerMarkers expected but localite.TriggerMarker not found');
        end
        triggerMarkers = [localite.TriggerMarker];
        recordingTimes = [triggerMarkers.recordingTimeAttribute];

    elseif strcmp(markertype, 'GUMMarkers')
        % ---- GUMMarkers: Entry/Target pairs with ColVec3D (3-vector) ------
        % XML structure: Element → EntryTargetPair → Target → Marker → ColVec3D
        %                                           → Rotation → RotationReference → ColVec3D
        % No full 4×4 matrix is stored; we synthesise one from the target
        % position and the rotation reference axis.
        if ~isfield(localite, 'Element')
            error('GUMMarkers expected but localite.Element not found');
        end
        elements = localite.Element;
        if ~iscell(elements), elements = num2cell(elements); end
        elem_n = numel(elements);
        triggerMarkers(elem_n) = struct();

        for i = 1:elem_n
            etp = elements{i}.EntryTargetPair;

            % Target position (ColVec3D data0/1/2 = X/Y/Z in RAS mm)
            tgt_vec = etp.Target.Marker.ColVec3D;
            pos = [tgt_vec.data0Attribute, tgt_vec.data1Attribute, tgt_vec.data2Attribute];

            % Rotation reference axis (unit vector pointing away from head)
            rot_ref = etp.Rotation.RotationReference.ColVec3D;
            z_axis = [rot_ref.data0Attribute, rot_ref.data1Attribute, rot_ref.data2Attribute];
            z_axis = z_axis / max(norm(z_axis), 1e-9);

            % Build an arbitrary but consistent 4×4 matrix:
            % column 4 = position, column 3 = z_axis (approach direction)
            % columns 1-2 = orthonormal complement (arbitrary in-plane frame)
            x_axis = null(z_axis)';   % 2×3; take first row as x
            x_axis = x_axis(1,:) / norm(x_axis(1,:));
            y_axis = cross(z_axis, x_axis);
            matrix_data = eye(4);
            matrix_data(1:3, 1) = x_axis';
            matrix_data(1:3, 2) = y_axis';
            matrix_data(1:3, 3) = z_axis';
            matrix_data(1:3, 4) = pos';

            triggerMarkers(i).Matrix4D_full = matrix_data;
            for f = 1:3
                triggerMarkers(i).(coord_fields{f}) = matrix_data(f, 4);
            end
            rot_flat = reshape(matrix_data(1:3, 1:3), 1, 9);
            for f = 1:9
                triggerMarkers(i).(rot_fields{f}) = rot_flat(f);
            end

            % Store description if present
            if isfield(etp, 'descriptionAttribute')
                triggerMarkers(i).description = etp.descriptionAttribute;
            end
        end
        recordingTimes = (1:elem_n)';

    else  % InstrumentMarker
        % ---- InstrumentMarker: static positions with full 4×4 matrix ------
        % XML structure: InstrumentMarker → Marker → Matrix4D
        if ~isfield(localite, 'InstrumentMarker')
            error('InstrumentMarker expected but localite.InstrumentMarker not found');
        end
        instr_markers = localite.InstrumentMarker;
        if ~iscell(instr_markers), instr_markers = num2cell(instr_markers); end
        elem_n = numel(instr_markers);
        triggerMarkers(elem_n) = struct();

        for i = 1:elem_n
            matrix4d = instr_markers{i}.Marker.Matrix4D;
            triggerMarkers(i).Matrix4D = matrix4d;

            % Store description if present
            if isfield(instr_markers{i}.Marker, 'descriptionAttribute')
                triggerMarkers(i).description = instr_markers{i}.Marker.descriptionAttribute;
            end
        end
        recordingTimes = (1:elem_n)';
    end

    % Expand raw Matrix4D structs to Matrix4D_full + coord/rot fields
    % (needed for TriggerMarkers and InstrumentMarker; GUMMarkers are already expanded above)
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

    % Detect segments (TriggerMarkers only — GUMMarkers and InstrumentMarker have no timing)
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
