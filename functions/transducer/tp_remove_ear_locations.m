function [locs] = tp_remove_ear_locations(parameters, locs)
%% Exclude ear entries from heuristic transducer positions
%  Removes transducer positions within tp_ear_radius of ear centers
%  Only runs if all three ear parameters are specified in parameters struct

% Only run if all three ear parameters are specified
if isfield(parameters, {'tp_ear_radius', 'tp_left_ear_center', 'tp_right_ear_center'})
    % Validate complete specification
    if ~isempty(parameters.tp_ear_radius) && ~isempty(parameters.tp_left_ear_center) && ~isempty(parameters.tp_right_ear_center)
        
        % Define nogo zone spheres
        ear_radius = parameters.tp_ear_radius;
        left_ear_center = parameters.tp_left_ear_center;
        right_ear_center = parameters.tp_right_ear_center;

        % Calculate Euclidean distances from each transducer position to both ears
        locs.dist_to_left_ear = sqrt((locs.trans_x - left_ear_center(1)).^2 + ...
                                    (locs.trans_y - left_ear_center(2)).^2 + ...
                                    (locs.trans_z - left_ear_center(3)).^2);

        locs.dist_to_right_ear = sqrt((locs.trans_x - right_ear_center(1)).^2 + ...
                                     (locs.trans_y - right_ear_center(2)).^2 + ...
                                     (locs.trans_z - right_ear_center(3)).^2);

        % Keep only positions outside both ear exclusion zones
        locs = locs(locs.dist_to_left_ear > ear_radius & locs.dist_to_right_ear > ear_radius, :);
        
        fprintf('[TP_REMOVE] Removed %d/%d positions in ear exclusion zones\n', ...
                height(locs_original) - height(locs), height(locs_original));
    else
        fprintf('[TP_REMOVE] Ear exclusion zones specified but incomplete - skipping\n');
    end
else
    fprintf('[TP_REMOVE] Ear exclusion parameters missing - all positions kept\n');
end

end
