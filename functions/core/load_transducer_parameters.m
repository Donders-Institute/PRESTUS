function parameters = load_transducer_parameters(parameters)
% LOAD_TRANSDUCER_PARAMETERS  Validate, normalize, and derive transducer parameters from config
%
% Iterates over all entries in parameters.transducer, dispatching to
% validate_matrix_transducer or validate_annular_transducer depending on
% tr.type. Derived scalar fields (elem_n, elem_phase_rad, focal distances)
% are computed and stored back into each transducer struct. Warns when
% multiple transducers have mismatched frequencies. Absorbs any legacy top-level
% expected_focal_distance_mm into per-transducer focal_distance_ep.
%
% Use as:
%   parameters = load_transducer_parameters(parameters)
%
% Input:
%   parameters - PRESTUS config containing a 'transducer' array of structs
%                with at minimum a 'type' field ('matrix' or 'annular')
%
% Output:
%   parameters - identical to input but with each transducer entry validated
%                and extended with derived fields including:
%                elem_n, elem_id_mm, elem_od_mm, elem_phase_rad, elem_phase_deg,
%                curv_radius_mm, dist_geom_ep_mm, focal_distance_ep,
%                focal_distance_bowl; matrix grids also get elem_n_row, elem_n_col,
%                sparsity_factor, clover fields
%
% See also: LOAD_PARAMETERS, VALIDATE_MATRIX_TRANSDUCER, VALIDATE_ANNULAR_TRANSDUCER,
%           FOCAL_DISTANCE_CALCULATION

arguments
    parameters (1,1) struct
end

    if isfield(parameters, 'transducer')

        % Check whether any transducer has been configured (non-empty type).
        % The default config contains a template entry with type: [] — skip
        % the entire validation block when no type has been set.
        has_type = arrayfun(@(t) isfield(t, 'type') && (ischar(t.type) || isstring(t.type)) && ~isempty(t.type), ...
                            parameters.transducer);
        if ~any(has_type)
            return;
        end

        new_transducers = struct([]);

        for t_i = 1:numel(parameters.transducer)

            tr = parameters.transducer(t_i);

            % yaml.loadFile returns MATLAB string objects; coerce to char for
            % consistent field access and switch/strcmp comparisons downstream.
            if isstring(tr.type)
                tr.type = char(tr.type);
            end

            % ---------------------------------------------------------------------
            % Validate and initialize transducer parameters
            % ---------------------------------------------------------------------

            % Supported transducer geometries: matrix and annular arrays

            % Ensure the type field is defined and is a non-empty string
            assert(isfield(tr, 'type') && ischar(tr.type) && ~isempty(tr.type),...
                'Transducer %i; Missing or empty type field. Please specify either "matrix" or "annular".', t_i);

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