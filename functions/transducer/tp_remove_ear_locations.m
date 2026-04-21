function [locs] = tp_remove_ear_locations(parameters, locs)
% TP_REMOVE_EAR_LOCATIONS  Exclude candidate positions near ear centres
%
% Removes transducer candidates that fall within a specified radius of
% the left or right ear centre. Only runs when all three ear parameters
% (ear_radius, left_ear_center, right_ear_center) are set.
%
% Use as:
%   locs = tp_remove_ear_locations(parameters, locs)
%
% Input:
%   parameters - (1,1) simulation parameters struct with placement.heuristic fields
%   locs       - table of candidate positions (from TP_EVALUATE_CANDIDATE_POSITIONS)
%
% Output:
%   locs - table with ear-region candidates removed
%
% See also: TP_EVALUATE_CANDIDATE_POSITIONS, TP_SELECT_HEURISTIC_POSITION

% Only run if all three ear parameters are specified
if isfield(parameters, 'placement') && isfield(parameters.placement, 'heuristic') && ...
        isfield(parameters.placement.heuristic, 'ear_radius') && ...
        isfield(parameters.placement.heuristic, 'left_ear_center') && ...
        isfield(parameters.placement.heuristic, 'right_ear_center')
    % Validate complete specification
    if ~isempty(parameters.placement.heuristic.ear_radius) && ~isempty(parameters.placement.heuristic.left_ear_center) && ~isempty(parameters.placement.heuristic.right_ear_center)

        % keep copy of original positions
        locs_original = locs;

        % Define nogo zone spheres
        ear_radius = parameters.placement.heuristic.ear_radius;
        left_ear_center = parameters.placement.heuristic.left_ear_center;
        right_ear_center = parameters.placement.heuristic.right_ear_center;

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
