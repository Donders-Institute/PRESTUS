function elem_phase_deg = set_real_phases(phase_table, tran, focal_distance_ep, parameters)
% SET_REAL_PHASES  Look up or compute manufacturer-specified steering phases for a focal depth
%
% For Sonic Concepts transducers, reads phase values directly from phase_table.
% For Imasonic transducers, calls COMPUTE_PHASES. Interpolates to virtual
% element count when the number of real and virtual elements differs.
%
% Use as:
%   elem_phase_deg = set_real_phases(phase_table, tran, focal_distance_ep, parameters)
%
% Input:
%   phase_table       - phase data table or ini struct (format depends on tran.manufact)
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

    % Extract or interpolate phases for the given focal depth
    if isequal(tran.manufact, "Sonic Concepts")
        % Extract phase values directly from the table
        distance = phase_table.Distance;
        phases = phase_table(:, 2:end-1); % Remove non-phase columns

        foc_index = find(distance == focal_distance_ep, 1);
        init_phases = table2array(phases(foc_index, :));
    elseif isequal(tran.manufact, "Imasonic")
        % Compute phases based on the focal depth
        init_phases = compute_phases(parameters.medium_properties.water.sound_speed, tran, focal_distance_ep, phase_table);
    else
        % Generic or other manufacturer: derive element positions analytically
        % from ring geometry (elem_id_mm, elem_od_mm, curv_radius_mm) in the
        % transducer YAML.  No external phase table file is required.
        tran_ini_data = generate_tran_ini_from_geometry(tran);
        init_phases = compute_phases(parameters.medium_properties.water.sound_speed, tran, focal_distance_ep, tran_ini_data);
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
