function parameters = load_transducer_parameters(parameters)

% LOAD_TRANSDUCER_PARAMETERS Validates and normalizes transducer parameters.
%
% This function processes the 'transducer' field of the input parameters struct.
% It validates required fields, ensures backward compatibility with legacy 
% annular-only configurations, and normalizes parameters into a consistent format 
% for simulation use. For matrix transducers, it validates element geometry, 
% curvature, steering, grid definitions, and optional Clover or sparse grid setups. 
% For annular transducers, it validates element diameters, curvature, and phases. 
% The function also calculates derived values such as source phases in radians, 
% element counts, and distances to the transducer plane.
%
% Input:
%   parameters - Struct containing simulation parameters, expected to include
%                a 'transducer' field with one or more transducer definitions.
%
% Output:
%   parameters - Struct with validated and normalized transducer parameters,
%                including additional derived fields such as:
%                * n_elements, Elements_ID_mm, Elements_OD_mm
%                * source_phase_rad, source_phase_deg
%                * curv_radius_mm, dist_to_plane_mm
%                * n_elem_row, n_elem_col (for matrix grids)
%                * sparsity_factor, clover setup fields, etc.

    if isfield(parameters, 'transducer')
    
        for t_i = 1:numel(parameters.transducer)
    
            tr = parameters.transducer(t_i);
    
            % ---------------------------------------------------------------------
            % Validate and initialize transducer parameters
            % ---------------------------------------------------------------------
                
            % Detect legacy configurations where only annular transducers are defined
            if isfield(tr, 'array_shape')
                % Supported transducer geometries: matrix and annular arrays
    
                % Ensure the array_shape.type field is defined
                assert(isfield(tr.array_shape, 'type'),...
                    'Transducer %i; Missing type field. Please specify either "matrix" or "annular".', t_i);
    
                switch tr.array_shape.type
    
                    case 'matrix'
                        [parameters, tr] = validate_matrix_transducer(parameters, tr, t_i);
    
                    case 'annular'
                        tr = validate_annular_transducer(tr, t_i);
    
                    otherwise
                        error('Transducer %i; Element shape option "%s" is not implemented.', ...
                            tr.array_shape.type);
                end
            else
                tr.array_shape.type = 'annular';
            end
    
            if strcmp(tr.array_shape.type, 'annular')

                annular_tr = tr.array_shape.annular;

                assert(numel(annular_tr.Elements_ID_mm) == annular_tr.n_elements, ...
                    'Transducer %i; Elements_ID_mm length must match n_elements.', t_i);
                
                assert(numel(annular_tr.Elements_OD_mm) == annular_tr.n_elements, ...
                    'Transducer %i; Elements_OD_mm length must match n_elements.', t_i);

                % Validate inner/outer diameter ordering
                assert(all(annular_tr.Elements_OD_mm > annular_tr.Elements_ID_mm), ...
                    'Transducer %i; Outer diameter must be larger than inner diameter for all elements.', t_i);
    
                % Ensure source phase is set in radians or degrees
                if ~isfield(tr,'source_phase_rad') && isfield(tr,'source_phase_deg')
                    tr.source_phase_rad = deg2rad(tr.source_phase_deg);
                elseif ~isfield(tr,'source_phase_rad')
                    error('Transducer %i; Phase must be specified as source_phase_rad or source_phase_deg.', t_i);
                end
    
                % Calculate distance to transducer plane if not provided
                if ~isfield(tr, 'dist_to_plane_mm')
                    assert(tr.curv_radius_mm > max(tr.Elements_OD_mm)/2, ...
                        'Transducer %i; curv_radius_mm must exceed aperture radius.', t_i);
                    
                    tr.dist_to_plane_mm = sqrt(tr.curv_radius_mm^2 - ...
                        (max(tr.Elements_OD_mm) / 2)^2);

                    fprintf('Transducer %i; Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                        t_i, tr.dist_to_plane_mm);
                end
                
                assert(isfield(tr,'source_amp'), ...
                    'Transducer %i; Missing source_amp field.', t_i);

                % Ensure source amplitude matches number of transducer elements
                if numel(tr.source_amp) == 1 && tr.n_elements > 1
                    tr.source_amp = repmat(tr.source_amp, [1, tr.n_elements]);
                end
    
                % Evaluate source phase expressions if stored as cell arrays
                if iscell(tr.source_phase_rad)
                    for p_i = 1:numel(tr.source_phase_rad)
                        if ~isnumeric(tr.source_phase_rad{p_i})
                            tr.source_phase_rad{p_i} = eval(tr.source_phase_rad{p_i});
                        end
                    end
                    tr.source_phase_rad = cell2mat(tr.source_phase_rad);
                end
    
                % Ensure source phase degrees are calculated from radians if not provided
                if ~isfield(tr, 'source_phase_deg')
                    tr.source_phase_deg = rad2deg(tr.source_phase_rad);
                end
            end
            
            parameters.transducer(t_i) = tr;
        
        end
    else
        % Warn user about missing transducer information
        assert(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No'), ...
            'Exiting');
    end

end