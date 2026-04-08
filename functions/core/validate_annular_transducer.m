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

    assert(isfield(tr, 'annular'), ...
        'Transducer %i; Appropriate configuration for annular transducer is missing.', t_i);

    annular_tr = tr.annular;
    
    % Validate required geometric parameters
    assert(isfield(annular_tr, 'elem_id_mm'), ...
        'Transducer %i; Missing elem_id_mm for annular transducer. This defines inner diameters of elements.', t_i);
    
    assert(isfield(annular_tr, 'elem_od_mm'), ...
        'Transducer %i; Missing elem_od_mm for annular transducer. This defines outer diameters of elements.', t_i);
    
    assert(isfield(annular_tr, 'elem_n'), ...
        'Transducer %i; Missing elem_n parameter for annular transducer. This defines element count.', t_i);
    
    assert(numel(annular_tr.elem_id_mm) == annular_tr.elem_n, ...
        'Transducer %i; elem_id_mm length must match elem_n.', t_i);
            
    assert(numel(annular_tr.elem_od_mm) == annular_tr.elem_n, ...
        'Transducer %i; elem_od_mm length must match elem_n.', t_i);

    % Validate inner/outer diameter ordering
    assert(all(annular_tr.elem_od_mm > annular_tr.elem_id_mm), ...
        'Transducer %i; Outer diameter must be larger than inner diameter for all elements.', t_i);
    
    % Validate curvature definition
    assert(isfield(annular_tr, 'curv_radius_mm'), ...
        'Transducer %i; Missing curv_radius_mm field for annular transducer. Please specify radius of curvature.', t_i);

    % 3D steering unavailable for annular arrays → align with focus
    tr.align_to_focus = true;

    % Calculate distance to transducer plane if not provided
    if ~isfield(annular_tr, 'dist_geom_ep_mm') || isempty(annular_tr.dist_geom_ep_mm)
        assert(annular_tr.curv_radius_mm > max(annular_tr.elem_od_mm)/2, ...
            'Transducer %i; curv_radius_mm must exceed aperture radius.', t_i);

        annular_tr.dist_geom_ep_mm = sqrt(annular_tr.curv_radius_mm^2 - ...
                                (max(annular_tr.elem_od_mm) / 2)^2);

        fprintf('Transducer %i; Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                t_i, annular_tr.dist_geom_ep_mm);
    end

    % Evaluate source phase expressions if stored as cell arrays
    if isfield(annular_tr, 'elem_phase_rad') && iscell(annular_tr.elem_phase_rad)
        for p_i = 1:numel(annular_tr.elem_phase_rad)
            if ~isnumeric(annular_tr.elem_phase_rad{p_i})
                annular_tr.elem_phase_rad{p_i} = eval(annular_tr.elem_phase_rad{p_i});
            end
        end
        annular_tr.elem_phase_rad = cell2mat(annular_tr.elem_phase_rad);
    end 

    if isfield(annular_tr, 'elem_phase_deg') && iscell(annular_tr.elem_phase_deg)
        for p_i = 1:numel(annular_tr.elem_phase_deg)
            if ~isnumeric(annular_tr.elem_phase_deg{p_i})
                annular_tr.elem_phase_deg{p_i} = eval(annular_tr.elem_phase_deg{p_i});
            end
        end
        annular_tr.elem_phase_deg = cell2mat(annular_tr.elem_phase_deg);
    end 

    % Validate that source phase is set in radians or degrees
    if ~isfield(annular_tr, 'elem_phase_rad') && isfield(annular_tr, 'elem_phase_deg')
        annular_tr.elem_phase_rad = deg2rad(annular_tr.elem_phase_deg);
    elseif ~isfield(annular_tr, 'elem_phase_deg') && isfield(annular_tr, 'elem_phase_rad')
        annular_tr.elem_phase_deg = rad2deg(annular_tr.elem_phase_rad);
    elseif ~isfield(annular_tr, 'elem_phase_rad') && ~isfield(annular_tr, 'elem_phase_deg')
        error('Transducer %i; Phase must be specified as elem_phase_rad or elem_phase_deg.', t_i);
    end 

    % Ensure source amplitude matches number of transducer elements
    if numel(annular_tr.elem_amp) == 1 && annular_tr.elem_n > 1
        annular_tr.elem_amp = repmat(annular_tr.elem_amp, [1, annular_tr.elem_n]);
    end

    % Encode updated annular transducer field
    tr.annular = annular_tr;
   
end