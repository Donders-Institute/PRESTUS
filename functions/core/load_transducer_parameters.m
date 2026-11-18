function parameters = load_transducer_parameters(parameters)
% LOAD_TRANSDUCER_PARAMETERS  Validate, normalize, and derive transducer parameters from config
%
% Iterates over all entries in parameters.transducer, dispatching to
% validate_matrix_transducer or validate_annular_transducer depending on
% tr.type. Derived scalar fields (elem_n, elem_phase_rad, focal distances)
% are computed and stored back into each transducer struct. Warns when a
% single transducer lists several frequencies. A legacy top-level
% expected_focal_distance_mm is an error (its reference point was ambiguous).
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
        has_serial = arrayfun(@(t) isfield(t, 'serial') && (ischar(t.serial) || isstring(t.serial)) && ~isempty(t.serial), ...
                              parameters.transducer);
        if ~any(has_type) && ~any(has_serial)
            return;
        end

        % ── Serial-based auto-resolution ─────────────────────────────────────
        % When a transducer entry carries a 'serial' field, load geometry from
        % the equipment library and, if phases/amplitude are absent, resolve
        % them from the calibration library.  Inline fields always take
        % precedence over library values.
        if any(has_serial)
            % Your own libraries (transducer.library_path, path.transducer_library,
            % PRESTUS_TRANSDUCER_LIBRARY) may hold equipment entries too.
            equip_param = load_equipment_config([], transducer_library_folders(parameters));
            library_path = fullfile(get_prestus_path(), 'config', 'calibration');
            % Resolve each serial-bearing transducer.  Because resolution may
            % add fields (e.g. type, freq_hz) that were absent from the
            % original struct array element, we rebuild the array from scratch
            % to avoid MATLAB's heterogeneous-struct-array assignment error.
            resolved = cell(1, numel(parameters.transducer));
            for t_i = 1:numel(parameters.transducer)
                resolved{t_i} = parameters.transducer(t_i);
            end
            for t_i = find(has_serial)
                before = resolved{t_i};
                resolved{t_i} = resolve_transducer_from_serial( ...
                    resolved{t_i}, equip_param, library_path, t_i, ...
                    transducer_library_folders(parameters, t_i));
                parameters = record_library_drive(parameters, t_i, before, resolved{t_i});
            end
            % Collect all field names across all elements, then rebuild.
            all_fields = cellfun(@fieldnames, resolved, 'UniformOutput', false);
            all_fields = unique(vertcat(all_fields{:}));
            rebuilt    = struct();
            for t_i = 1:numel(resolved)
                for fi = 1:numel(all_fields)
                    f = all_fields{fi};
                    if isfield(resolved{t_i}, f)
                        rebuilt(t_i).(f) = resolved{t_i}.(f);
                    else
                        rebuilt(t_i).(f) = [];
                    end
                end
            end
            parameters.transducer = rebuilt;
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
                    % no drive given: a nominal one, scaled to the target ISPPA
                    % (the phases steer to focus_pos geometrically)
                    if isfield(tr, 'matrix') && isstruct(tr.matrix) && ...
                            (~isfield(tr.matrix, 'elem_p_amp') || isempty(tr.matrix.elem_p_amp)) && ...
                            ~(isfield(tr.matrix, 'elem_amp') && ~isempty(tr.matrix.elem_amp)) && ...
                            isfield(tr, 'target_isppa_wcm2') && ~isempty(tr.target_isppa_wcm2)
                        tr.matrix.elem_p_amp = 100e3;
                        fprintf(['[transducer %d] Nominal drive (100 kPa per element); results are ' ...
                                 'scaled to target_isppa_wcm2.\n'], t_i);
                    end
                    [parameters, tr] = validate_matrix_transducer(parameters, tr, t_i);

                case 'annular'
                    % drive from the focal depth setting and target ISPPA where
                    % phases / amplitudes are not given (custom geometry)
                    c_water = 1500;
                    if isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'water') ...
                            && isfield(parameters.medium_properties.water, 'sound_speed')
                        c_water = parameters.medium_properties.water.sound_speed;
                    end
                    tr = fill_annular_drive(tr, c_water, t_i);
                    tr = validate_annular_transducer(tr, t_i);

                otherwise
                    error('prestus:transducer:unknownType', ...
                        'Transducer %i; Element shape option "%s" is not implemented.', ...
                        t_i, tr.type);
            end

            if t_i == 1
                new_transducers = tr;
            else
                new_transducers(t_i) = tr;
            end
        end

        parameters.transducer = new_transducers;

        % Transducers at different frequencies run as separate single-frequency
        % simulations (ASYNC_TRANSDUCER_PIPELINE), each on its own medium:
        %   async    - separate firing windows (intensities alternate)
        %   coherent - firing together (IS_SIMULTANEOUS_MULTIFREQ): complex fields
        %              summed per frequency (transducers at one frequency
        %              interfere), intensities summed across frequencies.
        % A single transducer with several frequencies gets a warning.
        all_freqs = [parameters.transducer.freq_hz];
        if numel(unique(all_freqs)) > 1 && numel(parameters.transducer) == 1
            warn('prestus:multiFreq', ...
                'Transducer has multiple drive frequencies (%s Hz). Verify MI uses the intended frequency.', ...
                num2str(unique(all_freqs)));
        end

        % Legacy expected_focal_distance_mm: its reference point (bowl apex or
        % exit plane) differed between code paths, so it cannot be converted
        if isfield(parameters, 'expected_focal_distance_mm') && ~isempty(parameters.expected_focal_distance_mm)
            error('prestus:legacyConfigKey', ...
                ['expected_focal_distance_mm is no longer read: set transducer.focal_distance_ep ' ...
                 '(from the exit plane) or transducer.focal_distance_bowl (from the bowl apex).']);
        end

        % Calculate focal distances (ep ↔ bowl) for all transducers
        parameters = focal_distance_calculation(parameters);

    else
        % Warn user about missing transducer information
        assert(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No'), ...
            'Exiting');
    end

end
function parameters = record_library_drive(parameters, t_i, before, after)
% parameters.state.library_drive{t_i}: which annular drive fields the
% calibration library filled in (phases, amplitudes), so that a run at another
% focal depth (e.g. TIME_REVERSAL_POSE_PARAMETERS) can clear and re-resolve
% them. A later call finds them filled and keeps the earlier record.
    filled = @(f) ~annular_set(before, f) && annular_set(after, f);
    lib = struct('phases', filled('elem_phase_deg'), 'amp', filled('elem_p_amp'));
    if ~isfield(parameters, 'state') || ~isfield(parameters.state, 'library_drive')
        parameters.state.library_drive = {};
    end
    if numel(parameters.state.library_drive) >= t_i && isstruct(parameters.state.library_drive{t_i})
        old = parameters.state.library_drive{t_i};
        lib.phases = lib.phases || old.phases;
        lib.amp    = lib.amp || old.amp;
    end
    parameters.state.library_drive{t_i} = lib;
end

function tf = annular_set(tr, f)
    tf = isfield(tr, 'annular') && isstruct(tr.annular) && isfield(tr.annular, f) && ...
        ~isempty(tr.annular.(f)) && isnumeric(tr.annular.(f)) && ~all(isnan(tr.annular.(f)(:)));
end
