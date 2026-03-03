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

        new_transducers = struct([]);
        
        for t_i = 1:numel(parameters.transducer)
    
            tr = parameters.transducer(t_i);
    
            % ---------------------------------------------------------------------
            % Validate and initialize transducer parameters
            % ---------------------------------------------------------------------
            
            % Detect legacy configurations where only annular transducers are defined
            if ~isfield(tr, 'array_shape')

                % Create a clean structure for the new format
                new_tr = struct();
                new_tr.array_shape.type = 'annular';
                new_tr.array_shape.annular = struct();  % initialize annular sub-struct;

                if isfield(tr, 'n_elements')
                    new_tr.array_shape.annular.n_elements = tr.n_elements;
                end

                if isfield(tr, 'Elements_ID_mm')
                    new_tr.array_shape.annular.Elements_ID_mm = tr.Elements_ID_mm;
                end

                if isfield(tr, 'Elements_OD_mm')
                    new_tr.array_shape.annular.Elements_OD_mm = tr.Elements_OD_mm;
                end

                if isfield(tr, 'curv_radius_mm')
                    new_tr.array_shape.annular.curv_radius_mm = tr.curv_radius_mm;
                end

                if isfield(tr, 'dist_to_plane_mm')
                    new_tr.array_shape.annular.dist_to_plane_mm = tr.dist_to_plane_mm;
                end

                if isfield(tr, 'source_phase_deg')
                    new_tr.array_shape.annular.source_phase_deg = tr.source_phase_deg;
                end

                % Replace old transducer completely
                tr = new_tr;
            end

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
    
            % Validate general parameters
            % Ensure source phase is set in radians or degrees
            if ~isfield(tr, 'source_phase_rad') && isfield(tr, 'source_phase_deg')
                tr.source_phase_rad = deg2rad(tr.source_phase_deg);
            elseif ~isfield(tr, 'source_phase_rad')
                error('Transducer %i; Phase must be specified as source_phase_rad or source_phase_deg.', t_i);
            end
			
			% Calculate distance between target and ep/bowl is not provided
			if ~isfield(parameters, 'expected_focal_distance_bowl') || ~isfield(parameters, 'expected_focal_distance_ep')
				parameters = focal_distance_calculation(parameters);
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
            
            if t_i == 1
                new_transducers = tr;
            else
                new_transducers(t_i) = tr;
            end
        end

        parameters.transducer = new_transducers;

    else
        % Warn user about missing transducer information
        assert(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No'), ...
            'Exiting');
    end

end