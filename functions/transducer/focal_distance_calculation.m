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

% Create the output fields on the whole array first: assigning an element
% that gained a new field back into a struct array would otherwise fail
% ("Subscripted assignment between dissimilar structures").
for f = {'focal_distance_offset', 'focal_distance_ep', 'focal_distance_bowl'}
    if ~isfield(parameters.transducer, f{1})
        [parameters.transducer.(f{1})] = deal([]);
    end
end

for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);

    has_ep   = isfield(tr, 'focal_distance_ep')   && ~isempty(tr.focal_distance_ep);
    has_bowl = isfield(tr, 'focal_distance_bowl') && ~isempty(tr.focal_distance_bowl);

    % Compute the geometry offset (bowl → exit plane) from type sub-struct.
    % Skip transducers that have no valid type or no matching sub-struct.
    tr_type = '';
    if isfield(tr, 'type') && (ischar(tr.type) || isstring(tr.type)) && ~isempty(tr.type)
        tr_type = char(tr.type);
    end
    focal_distance_offset = 0;
    if ~isempty(tr_type) && isfield(tr, tr_type)
        type_tr = tr.(tr_type);
        if isfield(type_tr, 'curv_radius_mm') && isfield(type_tr, 'dist_geom_ep_mm') && ...
                ~isempty(type_tr.curv_radius_mm) && ~isempty(type_tr.dist_geom_ep_mm) && ...
                isfinite(type_tr.curv_radius_mm)
            focal_distance_offset = type_tr.curv_radius_mm - type_tr.dist_geom_ep_mm;
        end
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
            warn('Transducer %d: trans_pos or focus_pos missing; cannot compute expected focal distance.', ti);
        elseif ~strcmp(manual_coordinate_space(parameters), 'voxel')
            % Positions are still RAS/MNI mm; preproc_transducer_placement
            % converts them to voxels and calls this function again
            disp('Positions not yet in T1 voxels; focal distance will be derived after placement.');
        else
            % trans_pos/focus_pos are T1 voxels: scale each axis by its voxel
            % size (non-integer and anisotropic T1 voxels)
            d_vox = tr.focus_pos(:)' - tr.trans_pos(:)';
            try
                t1_file = dir(fullfile(parameters.path.anat, sprintf(parameters.path.t1_pattern, parameters.subject_id)));
                t1_header = niftiinfo(fullfile(t1_file.folder, t1_file.name));
                voxel_mm = double(t1_header.PixelDimensions(1:numel(d_vox)));
                tr.focal_distance_bowl = norm(d_vox .* voxel_mm);
            catch
                warn('Could not load T1 header to determine resolution; assuming grid resolution.');
                tr.focal_distance_bowl = norm(d_vox) * parameters.grid.resolution_mm;
            end

            tr.focal_distance_ep = tr.focal_distance_bowl - focal_distance_offset;
        end
    end

    parameters.transducer(ti) = tr;
end

end
