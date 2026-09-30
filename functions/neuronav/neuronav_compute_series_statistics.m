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
%   voxel_size              - voxel size [mm]; used for position_dev_vox
%   expected_segment_length - expected pulses per series ([] = accept all)
%   markertype              - 'TriggerMarkers' (default), 'GUMMarkers', or 'InstrumentMarker'
%
% Output:
%   stats - 1×N cell array; each cell is a struct with fields:
%           series_id, N_pulses, voxel_size, position_mean [1×3],
%           matrix4d_mean [1×4×4], matrix4d_all [N×4×4],
%           position_dev_mm   - distance of each pulse from the series mean
%                               position [mm] (.mean, .std, .all)
%           position_dev_vox  - the same distances divided by voxel_size
%           rotation_dev_rad  - rotation angle of each pulse relative to the
%                               mean orientation [rad] (.mean, .std, .all)
%           markers
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
        % Some Localite versions mix in Element → InstrumentMarker → Matrix4D
        % entries (a recorded instrument pose); those carry a full 4×4 matrix.
        if ~isfield(localite, 'Element')
            error('GUMMarkers expected but localite.Element not found');
        end
        elements = localite.Element;
        if ~iscell(elements), elements = num2cell(elements); end
        elem_n = numel(elements);
        triggerMarkers(elem_n) = struct();

        for i = 1:elem_n
            % readstruct gives every Element both fields in a mixed file and
            % fills the absent one with <missing>, so test for a struct.
            if has_struct(elements{i}, 'InstrumentMarker')
                im = elements{i}.InstrumentMarker;
                matrix_data = matrix4d_to_full(im.Matrix4D);
                triggerMarkers(i).Matrix4D_full = matrix_data;
                for f = 1:3
                    triggerMarkers(i).(coord_fields{f}) = matrix_data(f, 4);
                end
                rot_flat = reshape(matrix_data(1:3, 1:3), 1, 9);
                for f = 1:9
                    triggerMarkers(i).(rot_fields{f}) = rot_flat(f);
                end
                if isfield(im, 'descriptionAttribute')
                    triggerMarkers(i).description = im.descriptionAttribute;
                end
                triggerMarkers(i).element_type = 'InstrumentMarker';
                continue;
            end
            if ~has_struct(elements{i}, 'EntryTargetPair')
                error('GUMMarkers Element %d is neither an EntryTargetPair nor an InstrumentMarker', i);
            end
            etp = elements{i}.EntryTargetPair;
            triggerMarkers(i).element_type = 'EntryTargetPair';

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
    % Note: do NOT guard on ~isfield(...,'Matrix4D_full') — setting that field on one
    % element of a MATLAB struct array propagates it (as []) to all other elements,
    % so isfield returns true for elements that were never actually populated.
    for i = 1:length(triggerMarkers)
        if isfield(triggerMarkers(i), 'Matrix4D')
            matrix_data = matrix4d_to_full(triggerMarkers(i).Matrix4D);
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

    % Build segments: one cell per series.
    % TriggerMarkers: detect series by timing gaps; each cell = cluster of triggers to average.
    % GUMMarkers / InstrumentMarker: each element is its own series (no averaging across targets).
    if strcmp(markertype, 'TriggerMarkers') && length(recordingTimes) > 1
        timeDiffs = diff(recordingTimes);
        gap_threshold = median(timeDiffs) * 3;
        jumpIdx = find(timeDiffs > gap_threshold);
        segmentStarts = [1, jumpIdx + 1];
        segmentEnds = [jumpIdx, numel(recordingTimes)];
    else
        segmentStarts = 1:length(triggerMarkers);
        segmentEnds   = 1:length(triggerMarkers);
    end

    % Create segments — one cell per series
    segments = {};
    for i = 1:length(segmentStarts)
        range = segmentStarts(i):segmentEnds(i);
        if strcmp(markertype, 'TriggerMarkers') && ...
           ~isempty(expected_segment_length) && abs(length(range) - expected_segment_length) > 20
            continue;
        end
        segments{end+1} = triggerMarkers(range); %#ok<AGROW>
    end

    if isempty(segments)
        warning('neuronav_compute_series_statistics: no valid %s segments found', markertype);
        stats = {};
        return;
    end

    % ----------- ONE CELL PER SERIES WITH Matrix4D_full AT CELL LEVEL -----------
    stats = {};  % 1×N cell array

    for s = 1:numel(segments)
        stim_data = segments{s};
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
        R_mean = mean_orientation(squeeze(mean(rotations, 1)));

        euclideanDists = sqrt(sum((coords - mean_coords).^2, 2));   % mm
        rot_deviations = zeros(N, 1);
        for i = 1:N
            % Recorded matrices are not exactly orthonormal: project each too.
            rot_deviations(i) = rotation_angle(R_mean, ...
                mean_orientation(reshape(rotations(i,:,:), 3, 3)));
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
            'position_dev_vox', struct(...
                'mean', mean(euclideanDists) / voxel_size, ...
                'std', std(euclideanDists) / voxel_size, ...
                'all', euclideanDists' / voxel_size), ...
            'rotation_dev_rad', struct(...
                'mean', mean(rot_deviations), ...
                'std', std(rot_deviations), ...
                'all', rot_deviations'), ...
            'markers', stim_data);                      % Original enriched markers
    end
    
    fprintf('Created %d series cells (%s)\n', length(stats), markertype);
end


function matrix_data = matrix4d_to_full(matrix4d)
% Parsed <Matrix4D data00=… data33=…> attributes → 4×4 double.
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
end


function R = mean_orientation(M)
% Nearest proper rotation to the element-wise mean of rotation matrices.
% Degenerate input (e.g. an all-zero rotation block) is returned as-is.
    [U, ~, V] = svd(M);
    D = eye(3);
    D(3,3) = sign(det(U * V'));
    if D(3,3) == 0, R = M; return; end
    R = U * D * V';
end


function theta = rotation_angle(R_ref, R)
% Angle [rad] of the relative rotation R_ref' * R (0 for identical poses).
% Chord form ||R_ref - R||_F = 2*sqrt(2)*sin(theta/2); unlike acos of the
% trace it stays accurate for small angles.
    theta = 2 * asin(min(norm(R_ref - R, 'fro') / (2 * sqrt(2)), 1));
end


function tf = has_struct(s, field)
    tf = isfield(s, field) && isstruct(s.(field));
end
