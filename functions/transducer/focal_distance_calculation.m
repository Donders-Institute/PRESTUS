function parameters = focal_distance_calculation(parameters)
% Compute expected focal distances for (multi-)transducer setup
%
% DESCRIPTION:
%   Ensures all transducers in parameters.transducer have valid focal_distance_ep
%   and focal_distance_bowl fields at the transducer top level.
%   Derives missing values from the counterpart (ep↔bowl) using the transducer
%   geometry offset, or falls back to geometric distance from trans_pos/focus_pos.
%
% INPUTS:
%   parameters  - (struct) Simulation configuration with:
%                   .transducer - (1×N struct array) each element may contain:
%                     .focal_distance_ep   - [mm] focal distance from exit plane (optional)
%                     .focal_distance_bowl - [mm] focal distance from bowl (optional)
%                     .trans_pos        - [voxels] transducer position in T1 grid (optional)
%                     .focus_pos        - [voxels] target position in T1 grid (optional)
%
% OUTPUTS:
%   parameters  - (struct) Updated with focal_distance_ep/bowl on all transducers.

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
