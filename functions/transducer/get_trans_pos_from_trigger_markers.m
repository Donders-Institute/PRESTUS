function [transducer_pos, target_pos] = get_trans_pos_from_trigger_markers(trigger_markers_file, trigger_index, reference_to_transducer_distance, reference_to_target_distance)

% GET_TRANS_POS_FROM_TRIGGER_MARKERS Computes transducer and target positions from Localite trigger marker files.
%
% This function calculates the positions of the transducer and target based on 
% a specified trigger marker in a Localite XML file. The transducer position is 
% computed relative to a reference point, and the target position is optionally 
% calculated based on an additional distance parameter.
%
% Input:
%   trigger_markers_file          - Path to the Localite XML file containing trigger markers.
%   trigger_index                 - Index of the trigger marker to use for position calculations.
%   reference_to_transducer_distance - Distance from the reference point to the transducer (in mm).
%   reference_to_target_distance  - Distance from the reference point to the target (optional, default: 0 mm).
%
% Output:
%   transducer_pos                - [1x3] array specifying the transducer position in Cartesian coordinates.
%   target_pos                    - [1x3] array specifying the target position in Cartesian coordinates.
%
% Notes:
%   - Assumes that the center of the radius of the transducer is the target.
%   - Requires `xml2struct` (https://github.com/joe-of-all-trades/xml2struct) for parsing XML files.

    arguments
        trigger_markers_file
        trigger_index
        reference_to_transducer_distance
        reference_to_target_distance = 0; % Default distance to target is 0 mm
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
