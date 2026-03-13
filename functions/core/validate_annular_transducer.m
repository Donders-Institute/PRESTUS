function validate_annular_transducer(tr, t_i)
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
    
    % Validate curvature definition
    assert(isfield(annular_tr, 'curv_radius_mm'), ...
        'Transducer %i; Missing curv_radius_mm field for annular transducer. Please specify radius of curvature.', t_i);
    
    % Validate source phase definition
    assert(isfield(annular_tr, 'source_phase_deg'), ...
        'Transducer %i;: Missing source_phase_deg field for annular transducer. Please specify phases.', t_i);
    
end