function [transducer_pos target_pos] = get_trans_pos_from_trigger_markers(trigger_markers_file, trigger_index, reference_to_transducer_distance, reference_to_target_distance )
    
    % Finds the position of the transducer and target based on the provided
    % localite trigger marker files.
    % Assumes that the center of the radius of the transducer is the target
    %   Meaning that a mask should be loaded in separately to perform
    %   calculations on simulation results in your actual target

    % requires https://github.com/joe-of-all-trades/xml2struct

    arguments
    trigger_markers_file
    trigger_index
    reference_to_transducer_distance
    reference_to_target_distance = 0;
    end
    xml = xml2struct(trigger_markers_file);

    % find i-th marker
    trigger_counter = 0;
    for i = 1:length(xml.Children)
        cur_child = xml.Children(i);
        if strcmp(cur_child.Name, 'TriggerMarker')
            trigger_counter=trigger_counter+1;
        end
        if trigger_counter==trigger_index
            break
        end
    end

    coord_matrix = str2double({cur_child.Children(4).Attributes.Value});
    coord_matrix = reshape(coord_matrix',[4,4])';

    reference_pos = coord_matrix(:,4); % Position of the reference
    reference_center_to_head = coord_matrix(:,1); % From the center of the coil towards head

    transducer_pos = reference_pos + reference_to_transducer_distance*reference_center_to_head ;
    transducer_pos = transducer_pos(1:3);
    
    target_pos = reference_pos + reference_to_target_distance*reference_center_to_head ;
    target_pos = target_pos(1:3);
end