function parameters = focal_distance_calculation(parameters)
% Compute expected focal distances for (multi-)transducer setup
%
% DESCRIPTION:
%   Ensures all transducers in parameters.transducer have a valid exp_FD_bowl field.
%   Supports legacy single-value propagation, transducer-specific values, and automatic derivation 
%   from T1-weighted MRI grid positions when needed. Essential preprocessing for acoustic simulations 
%   requiring precise geometric focus-targeting.
%
% INPUTS:
%   parameters  - (struct) Simulation configuration with:
%                   .transducer          - (1×N struct array) Transducer configurations where each element may contain:
%                     .exp_FD_bowl  - (scalar) [mm] Expected geometric focal distance (optional)
%                     .trans_pos        - (1×3 vector) Transducer position [voxels] in T1 grid (optional)
%                     .focus_pos  - (1×3 vector) Target focus position [voxels] in T1 grid (optional)
%                   .exp_FD_bowl  - (scalar) [mm] Legacy global focal distance (optional)
%                   .data_path           - (char) Path to subject data
%                   .t1_path_template    - (char) sprintf template for T1w filename e.g. 'sub-%03d_T1w.nii.gz'
%
% OUTPUTS:
%   parameters  - (struct) Updated with position information

for ti = 1:numel(parameters.transducer)
    tr = parameters.transducer(ti);

    if ~isfield(tr.position, 'exp_FD_bowl') || isempty(tr.position.exp_FD_bowl) || ...
            ~isfield(tr.position, 'exp_FD_ep') || isempty(tr.position.exp_FD_ep)
        
        disp('Expected focal distance (EP or bowl) not specified, trying to get it from transducer and target positions ...')
        
        % calculate distance from bowl (if not specified)
        if ~isfield(tr.position, 'exp_FD_bowl') || isempty(tr.position.exp_FD_bowl)

            if ~isfield(tr.position, 'trans_pos')  || isempty(tr.position.trans_pos) || ...
               ~isfield(tr.position, 'focus_pos') || isempty(tr.position.focus_pos)
                warning('Transducer %d: trans_pos or focus_pos missing; cannot compute expected focal distance.', ti);
            end

            % calculate grid distance between transducer bowl and focus
            focal_distance = norm(tr.position.focus_pos - tr.position.trans_pos);
    
            % scale grid distance to calculate mm
            % Usually: T1 planning image dimensions; in phantom/water this may be the grid size.
            try
                % get T1 resolution from the header
                t1_file = dir(fullfile(parameters.data_path, sprintf(parameters.t1_path_template, parameters.subject_id)));
                t1_header = niftiinfo(fullfile(t1_file.folder, t1_file.name));
                t1_resolution_mm = round(mean(t1_header.PixelDimensions(1:3)));
                tr.position.exp_FD_bowl = ...
                    focal_distance * t1_resolution_mm;
            catch
                warning('Could not load T1 header to determine resolution; cannot compute expected focal distance in mm. Assuming grid size.');
                tr.position.exp_FD_bowl = ...
                    focal_distance * parameters.grid.resolution_mm;
            end
        end

        % calculate focal distance offset
        % (between transducer bowl and exit plane for annular/curved arrays)
        if ~isfield(tr.(tr.type), 'focal_distance_offset') || isempty(tr.(tr.type).focal_distance_offset)
            if strcmp(tr.type, 'annular')
                tr.annular.focal_distance_offset = ...
                    tr.annular.curv_radius_mm - ...
                    tr.annular.dist_to_plane_mm;
            elseif strcmp(tr.type, 'matrix') && tr.matrix.is_curved == 1
                tr.matrix.focal_distance_offset = ...
                    tr.matrix.curv_radius_mm - ...
                    tr.matrix.dist_to_plane_mm;
            else
                tr.(tr.type).focal_distance_offset = 0;
            end
        end

        % calculate focal distance (from exit plane)
        tr.position.exp_FD_ep = ...
            tr.position.exp_FD_bowl-...
            tr.(tr.type).focal_distance_offset;

    end

    % Update transducer configuration
    parameters.transducer(ti) = tr;

end

end
