function elem_phase_deg = set_real_phases(phase_table, tran, focal_distance_ep, parameters)
% SET_REAL_PHASES  Look up or compute manufacturer-specified steering phases for a focal depth
%
% Computes per-element phases for an annular transducer at a given focal depth.
% When phase_table is empty ([]) or not provided, phases are derived purely from
% the ring geometry (elem_id_mm, elem_od_mm, curv_radius_mm) in the equipment
% YAML — no external file required.  This is the default path for all new
% transducers.
%
% When a phase_table is supplied:
%   Imasonic (.ini struct)       — element positions from .ini, path-length delays
%   Sonic Concepts (table)       — phases interpolated from manufacturer steer table
%
% Use as:
%   elem_phase_deg = set_real_phases(phase_table, tran, focal_distance_ep, parameters)
%   elem_phase_deg = set_real_phases([],          tran, focal_distance_ep, parameters)
%
% Input:
%   phase_table       - phase data table, ini struct, or [] for geometry-only mode
%   tran              - transducer struct with manufact, n_elem, transducer fields
%   focal_distance_ep - desired focal depth relative to exit plane [mm]
%   parameters        - PRESTUS config; uses medium_properties.water.sound_speed [m/s]
%
% Output:
%   elem_phase_deg - [1xN] element phases [°]
%
% See also: COMPUTE_PHASES, GENERATE_TRAN_INI_FROM_GEOMETRY, SAVE_OPTIMIZED_VALUES, CALIBRATION_TRANSDUCER

arguments
    phase_table
    tran              (1,1) struct
    focal_distance_ep (1,1) {mustBeNumeric}
    parameters        (1,1) struct
end

    use_table = ~isempty(phase_table);

    if ~use_table
        % Geometry-only mode: no external file required.
        tran_ini_data = generate_tran_ini_from_geometry(tran);
        init_phases = compute_phases(parameters.medium_properties.water.sound_speed, ...
            tran, focal_distance_ep, tran_ini_data);

    elseif isequal(tran.manufact, "Sonic Concepts")
        % Interpolate phases from manufacturer steer table.
        distance = phase_table.Distance;
        phase_mat = table2array(phase_table(:, 2:end-1));

        [~, i_lo] = min(abs(distance - focal_distance_ep));
        if distance(i_lo) == focal_distance_ep || numel(distance) == 1
            init_phases = phase_mat(i_lo, :);
        else
            % Linear interpolation between the two nearest rows.
            dists_sorted = sort(distance);
            idx_lo = find(dists_sorted <= focal_distance_ep, 1, 'last');
            idx_hi = find(dists_sorted >= focal_distance_ep, 1, 'first');
            if idx_lo == idx_hi
                init_phases = phase_mat(idx_lo, :);
            else
                t = (focal_distance_ep - dists_sorted(idx_lo)) / ...
                    (dists_sorted(idx_hi) - dists_sorted(idx_lo));
                init_phases = (1 - t) * phase_mat(idx_lo, :) + t * phase_mat(idx_hi, :);
            end
        end

    elseif isequal(tran.manufact, "Imasonic")
        % Element positions from .ini file, then path-length delays.
        init_phases = compute_phases(parameters.medium_properties.water.sound_speed, ...
            tran, focal_distance_ep, phase_table);

    else
        % Unknown manufacturer with a phase_table supplied — use geometry fallback
        % and warn so the caller knows the table was ignored.
        warning('set_real_phases:unknownManufacturer', ...
            'Unknown manufacturer ''%s'': ignoring phase_table, using geometry-only calculation.', ...
            tran.manufact);
        tran_ini_data = generate_tran_ini_from_geometry(tran);
        init_phases = compute_phases(parameters.medium_properties.water.sound_speed, ...
            tran, focal_distance_ep, tran_ini_data);
    end

    % Adjust phases for virtual elements if necessary
    if tran.n_elem ~= tran.transducer.annular.elem_n
        % Map real elements to virtual elements
        original_pos = linspace(1, tran.n_elem, tran.n_elem);

        % Convert phases to radians, unwrap, and interpolate
        unwrapped_phases = unwrap(init_phases * pi / 180) * 180 / pi;
        virtual_pos = linspace(1, tran.n_elem, tran.transducer.annular.elem_n);
        virtual_phases = interp1(original_pos, unwrapped_phases, virtual_pos, 'linear');

        % Convert back to degrees and apply modulo operation
        init_phases = mod(virtual_phases, 360);
    end

    % Output source phase
    elem_phase_deg = init_phases;

end
