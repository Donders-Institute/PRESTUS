function [transducer_pos, target_pos] = get_trans_pos_from_trigger_markers(trigger_markers_file, trigger_index, reference_to_transducer_distance, reference_to_target_distance)

% GET_TRANS_POS_FROM_TRIGGER_MARKERS  Compute transducer/target positions from a Localite XML file
%
% Calculates transducer and target positions from a specified trigger marker
% in a Localite XML file. The transducer position is derived by projecting
% a reference point along the approach direction by the given distance.
%
% Use as:
%   [transducer_pos, target_pos] = get_trans_pos_from_trigger_markers( ...
%       trigger_markers_file, trigger_index, reference_to_transducer_distance)
%   [transducer_pos, target_pos] = get_trans_pos_from_trigger_markers( ...
%       trigger_markers_file, trigger_index, reference_to_transducer_distance, ...
%       reference_to_target_distance)
%
% Input:
%   trigger_markers_file             - path to Localite XML file with trigger markers
%   trigger_index                    - index of the trigger marker to use
%   reference_to_transducer_distance - distance from reference point to transducer [mm]
%   reference_to_target_distance     - distance from reference point to target [mm] (default: 0)
%
% Output:
%   transducer_pos - [1x3] transducer position in Cartesian coordinates
%   target_pos     - [1x3] target position in Cartesian coordinates
%
% See also: POSITION_TRANSDUCER_LOCALITE, LOCALITE_MATRIX_TO_POSITIONS

    arguments
        trigger_markers_file             (1,:) char
        trigger_index                    (1,1) double
        reference_to_transducer_distance (1,1) double
        reference_to_target_distance     (1,1) double = 0
    end

    % Parse the XML file using `xml2struct`
    xml = xml2struct(trigger_markers_file);

    % Locate the specified trigger marker by its index
    trigger_counter = 0;
    for i = 1:length(xml.Children)
        cur_child = xml.Children(i);
        if strcmp(cur_child.Name, 'TriggerMarker')
            trigger_counter = trigger_counter + 1; % Increment counter for each TriggerMarker found
        end
        if trigger_counter == trigger_index
            break; % Stop once we reach the desired marker index
        end
    end

    % Extract transformation matrix from the marker's attributes
    coord_matrix = str2double({cur_child.Children(4).Attributes.Value});
    coord_matrix = reshape(coord_matrix', [4, 4])'; % Reshape into a 4x4 matrix

    % Extract reference position and direction vector from transformation matrix
    reference_pos = coord_matrix(:, 4); % Position of the reference point (translation vector)
    reference_center_to_head = coord_matrix(:, 1); % Direction vector from center of coil towards head

    % Compute transducer position based on distance from reference point
    transducer_pos = reference_pos + reference_to_transducer_distance * reference_center_to_head;
    transducer_pos = transducer_pos(1:3); % Extract x, y, z coordinates

    % Compute target position based on optional distance parameter
    target_pos = reference_pos + reference_to_target_distance * reference_center_to_head;
    target_pos = target_pos(1:3); % Extract x, y, z coordinates
end
