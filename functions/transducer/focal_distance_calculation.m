function parameters = focal_distance_calculation(parameters)
% FOCAL_DISTANCE_CALCULATION  Compute expected focal distances for all transducers
%
% Ensures every entry in parameters.transducer has valid focal_distance_ep
% and focal_distance_bowl fields. Derives missing values from the available
% counterpart using the transducer geometry offset (ep ↔ bowl), or falls
% back to the geometric distance between trans_pos and focus_pos.
%
% Use as:
%   parameters = focal_distance_calculation(parameters)
%
% Input:
%   parameters - (1,1) simulation configuration struct with .transducer array
%
% Output:
%   parameters - updated struct with focal_distance_ep and focal_distance_bowl
%                set on every transducer
%
% See also: TRANSDUCER_SETUP, TRANSDUCER_POSITIONING

arguments
    parameters (1,1) struct
end

for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);

    has_ep   = isfield(tr, 'focal_distance_ep')   && ~isempty(tr.focal_distance_ep);
    has_bowl = isfield(tr, 'focal_distance_bowl') && ~isempty(tr.focal_distance_bowl);

    % Compute the geometry offset (bowl → exit plane) from type sub-struct
    type_tr = tr.(tr.type);
    if isfield(type_tr, 'curv_radius_mm') && isfield(type_tr, 'dist_geom_ep_mm') && ...
            ~isempty(type_tr.curv_radius_mm) && ~isempty(type_tr.dist_geom_ep_mm) && ...
            isfinite(type_tr.curv_radius_mm)
        focal_distance_offset = type_tr.curv_radius_mm - type_tr.dist_geom_ep_mm;
    else
        focal_distance_offset = 0;
    end
    tr.focal_distance_offset = focal_distance_offset;

    if has_ep && ~has_bowl
        tr.focal_distance_bowl = tr.focal_distance_ep + focal_distance_offset;

    elseif has_bowl && ~has_ep
        tr.focal_distance_ep = tr.focal_distance_bowl - focal_distance_offset;

    elseif ~has_ep && ~has_bowl
        % Fall back: compute geometric distance from trans_pos / focus_pos
        disp('Expected focal distance (EP or bowl) not specified, trying to get it from transducer and target positions ...')

        has_pos = isfield(tr, 'trans_pos') && ~isempty(tr.trans_pos) && ...
                  isfield(tr, 'focus_pos') && ~isempty(tr.focus_pos);

        if ~has_pos
            warning('Transducer %d: trans_pos or focus_pos missing; cannot compute expected focal distance.', ti);
        else
            focal_distance = norm(tr.focus_pos - tr.trans_pos);

            try
                t1_file = dir(fullfile(parameters.data_path, sprintf(parameters.t1_path_template, parameters.subject_id)));
                t1_header = niftiinfo(fullfile(t1_file.folder, t1_file.name));
                t1_resolution_mm = round(mean(t1_header.PixelDimensions(1:3)));
                tr.focal_distance_bowl = focal_distance * t1_resolution_mm;
            catch
                warning('Could not load T1 header to determine resolution; assuming grid resolution.');
                tr.focal_distance_bowl = focal_distance * parameters.grid.resolution_mm;
            end

            tr.focal_distance_ep = tr.focal_distance_bowl - focal_distance_offset;
        end
    end

    parameters.transducer(ti) = tr;
end

end
