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
            if ~isfield(tr, 'type') || isempty(tr.type)

                % Create a clean structure for the new format
                new_tr = struct();
                new_tr.type = 'annular';
                new_tr.annular = struct();  % initialize annular sub-struct;

                if isfield(tr, 'n_elements')
                    new_tr.annular.n_elements = tr.n_elements;
                end

                if isfield(tr, 'Elements_ID_mm')
                    new_tr.annular.Elements_ID_mm = tr.Elements_ID_mm;
                end

                if isfield(tr, 'Elements_OD_mm')
                    new_tr.annular.Elements_OD_mm = tr.Elements_OD_mm;
                end

                if isfield(tr, 'curv_radius_mm')
                    new_tr.annular.curv_radius_mm = tr.curv_radius_mm;
                end

                if isfield(tr, 'dist_to_plane_mm')
                    new_tr.annular.dist_to_plane_mm = tr.dist_to_plane_mm;
                end

                if isfield(tr, 'source_amp')
                    new_tr.source_amp = tr.source_amp;
                end

                if isfield(tr, 'source_phase_deg')
                    new_tr.source_phase_deg = tr.source_phase_deg;
                end
				
                if isfield(tr, 'source_freq_hz')
                    new_tr.source_freq_hz = tr.source_freq_hz;
                end

                new_tr.position = struct();
                if isfield(tr, 'trans_pos')
                    new_tr.position.trans_pos = tr.trans_pos;
                end

                if isfield(tr, 'focus_pos')
                    new_tr.position.focus_pos = tr.focus_pos;
                end

                if isfield(tr, 'expected_focal_distance_ep')
                    new_tr.position.exp_FD_ep = tr.expected_focal_distance_ep;
                end

                if isfield(tr, 'expected_focal_distance_bowl')
                    new_tr.position.exp_FD_bowl = tr.expected_focal_distance_bowl;
                end

                % Replace old transducer completely
                tr = new_tr;
            end

            % Supported transducer geometries: matrix and annular arrays

            % Ensure the type field is defined
            assert(isfield(tr, 'type'),...
                'Transducer %i; Missing type field. Please specify either "matrix" or "annular".', t_i);

            switch tr.type

                case 'matrix'
                    [parameters, tr] = validate_matrix_transducer(parameters, tr, t_i);

                case 'annular'
                    tr = validate_annular_transducer(tr, t_i);

                otherwise
                    error('Transducer %i; Element shape option "%s" is not implemented.', ...
                        tr.type);
            end

             % Ensure distance between target and ep/bowl is provided
            parameters = focal_distance_calculation(parameters);
			
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