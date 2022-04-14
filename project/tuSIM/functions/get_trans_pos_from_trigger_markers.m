function transducer_pos = get_trans_pos_from_trigger_markers(trigger_markers_file, trigger_index, reference_to_transducer_distance )

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
end