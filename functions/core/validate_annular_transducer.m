function tr = validate_annular_transducer(tr, t_i)
% VALIDATE_ANNULAR_TRANSDUCER Validates configuration of an annular transducer.
%
% This function checks that all required fields for an annular transducer
% definition are present in the input transducer structure.
%
% Input:
%   tr  - Struct containing a single transducer definition.
%   t_i - Index of the transducer in the configuration (used for error messages).
%
% Output:
%   (None)

    assert(isfield(tr.array_shape, 'annular'), ...
        'Transducer %i; Appropriate configuration for annular transducer is missing.', t_i);
    
    annular_tr = tr.array_shape.annular;
    
    % Validate required geometric parameters
    assert(isfield(annular_tr, 'Elements_ID_mm'), ...
        'Transducer %i; Missing Elements_ID_mm for annular transducer. This defines inner diameters of elements.', t_i);
    
    assert(isfield(annular_tr, 'Elements_OD_mm'), ...
        'Transducer %i; Missing Elements_OD_mm for annular transducer. This defines outer diameters of elements.', t_i);
    
    assert(isfield(annular_tr, 'n_elements'), ...
        'Transducer %i; Missing n_elements parameter for annular transducer. This defines element count.', t_i);

    tr.n_elements = annular_tr.n_elements;
    
    assert(numel(annular_tr.Elements_ID_mm) == annular_tr.n_elements, ...
        'Transducer %i; Elements_ID_mm length must match n_elements.', t_i);
            
    assert(numel(annular_tr.Elements_OD_mm) == annular_tr.n_elements, ...
        'Transducer %i; Elements_OD_mm length must match n_elements.', t_i);

    % Validate inner/outer diameter ordering
    assert(all(annular_tr.Elements_OD_mm > annular_tr.Elements_ID_mm), ...
        'Transducer %i; Outer diameter must be larger than inner diameter for all elements.', t_i);
    
    % Validate curvature definition
    assert(isfield(annular_tr, 'curv_radius_mm'), ...
        'Transducer %i; Missing curv_radius_mm field for annular transducer. Please specify radius of curvature.', t_i);

    tr.curv_radius_mm = annular_tr.curv_radius_mm;
    
    % 3D steering unavailable for annular arrays → align with focus
    tr.align_transducer_with_focus = true;
                        
    % Calculate distance to transducer plane if not provided
    if ~isfield(tr, 'dist_to_plane_mm')
        assert(annular_tr.curv_radius_mm > max(annular_tr.Elements_OD_mm)/2, ...
            'Transducer %i; curv_radius_mm must exceed aperture radius.', t_i);
            
        tr.dist_to_plane_mm = sqrt(annular_tr.curv_radius_mm^2 - ...
                                (max(annular_tr.Elements_OD_mm) / 2)^2);

        fprintf('Transducer %i; Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                t_i, tr.dist_to_plane_mm);
    else
        tr.dist_to_plane_mm = annular_tr.dist_to_plane_mm;
    end

    % Evaluate source phase expressions if stored as cell arrays
    if isfield(annular_tr, 'source_phase_rad') && iscell(annular_tr.source_phase_rad)
        for p_i = 1:numel(annular_tr.source_phase_rad)
            if ~isnumeric(annular_tr.source_phase_rad{p_i})
                annular_tr.source_phase_rad{p_i} = eval(annular_tr.source_phase_rad{p_i});
            end
        end
        annular_tr.source_phase_rad = cell2mat(annular_tr.source_phase_rad);
    end 

    if isfield(annular_tr, 'source_phase_deg') && iscell(annular_tr.source_phase_deg)
        for p_i = 1:numel(annular_tr.source_phase_deg)
            if ~isnumeric(annular_tr.source_phase_deg{p_i})
                annular_tr.source_phase_deg{p_i} = eval(annular_tr.source_phase_deg{p_i});
            end
        end
        annular_tr.source_phase_deg = cell2mat(annular_tr.source_phase_deg);
    end 

    % Validate that source phase is set in radians or degrees
    if ~isfield(annular_tr, 'source_phase_rad') && isfield(annular_tr, 'source_phase_deg')
        annular_tr.source_phase_rad = deg2rad(annular_tr.source_phase_deg);
    elseif ~isfield(annular_tr, 'source_phase_deg') && isfield(annular_tr, 'source_phase_rad')
        annular_tr.source_phase_deg = rad2deg(annular_tr.source_phase_rad);
    elseif ~isfield(annular_tr, 'source_phase_rad') && ~isfield(annular_tr, 'source_phase_deg')
        error('Transducer %i; Phase must be specified as source_phase_rad or source_phase_deg.', t_i);
    end

    % Encode updated annular transducer field
    tr.array_shape.annular = annular_tr;
   
end