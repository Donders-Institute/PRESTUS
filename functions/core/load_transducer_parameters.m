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
%                * elem_n, elem_id_mm, elem_od_mm
%                * elem_phase_rad, elem_phase_deg
%                * curv_radius_mm, dist_geom_ep_mm
%                * elem_n_row, elem_n_col (for matrix grids)
%                * sparsity_factor, clover setup fields, etc.

    if isfield(parameters, 'transducer')

        new_transducers = struct([]);
        
        for t_i = 1:numel(parameters.transducer)
    
            tr = parameters.transducer(t_i);
    
            % ---------------------------------------------------------------------
            % Validate and initialize transducer parameters
            % ---------------------------------------------------------------------

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

            if t_i == 1
                new_transducers = tr;
            else
                new_transducers(t_i) = tr;
            end
        end

        parameters.transducer = new_transducers;

        % Warn if transducers have mismatched frequencies
        all_freqs = [parameters.transducer.freq_hz];
        if numel(unique(all_freqs)) > 1
            warning('Transducers have different frequencies: %s Hz. This is not fully supported.', ...
                num2str(unique(all_freqs)));
        end

        % Legacy: top-level expected_focal_distance_mm was ep distance
        if isfield(parameters, 'expected_focal_distance_mm') && ~isempty(parameters.expected_focal_distance_mm)
            for ti = 1:numel(parameters.transducer)
                if ~isfield(parameters.transducer(ti), 'focal_distance_ep') || ...
                        isempty(parameters.transducer(ti).focal_distance_ep)
                    parameters.transducer(ti).focal_distance_ep = parameters.expected_focal_distance_mm;
                end
            end
        end

        % Calculate focal distances (ep ↔ bowl) for all transducers
        parameters = focal_distance_calculation(parameters);

    else
        % Warn user about missing transducer information
        assert(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No'), ...
            'Exiting');
    end

end