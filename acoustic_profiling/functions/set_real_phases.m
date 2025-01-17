function sim_param = set_real_phases(phase_table, tran, focus_wrt_exit_plane, sim_param, SOUND_SPEED_WATER)
    % Set the real phases of the transducer for a given focal depth.
    %
    % Arguments:
    % - phase_table: Additional data for phase calculations.
    % - tran: Transducer structure containing manufacturer and element details.
    % - focus_wrt_exit_plane: Desired focal depth with respect to the transducer exit plane [mm].
    % - sim_param: Simulation parameters structure.
    % - SOUND_SPEED_WATER: Speed of sound in water [m/s].
    %
    % Returns:
    % - sim_param: Updated simulation parameters with the new phase values.

    % Extract or interpolate phases for the given focal depth
    if isequal(tran.manufact, "Sonic Concepts")
        % Extract phase values directly from the table
        distance = phase_table.Distance;
        phases = phase_table(:, 2:end-1); % Remove non-phase columns
        
        foc_index = find(distance == focus_wrt_exit_plane, 1);
        init_phases = table2array(phases(foc_index, :));
    elseif isequal(tran.manufact, "Imasonic")
        % Compute phases based on the focal depth
        init_phases = compute_phases(SOUND_SPEED_WATER, tran, focus_wrt_exit_plane, phase_table);
    else
        error('Unsupported transducer manufacturer: %s', tran.manufact);
    end

    % Adjust phases for virtual elements if necessary
    if tran.n_elem ~= tran.prestus.transducer.n_elements
        % Map real elements to virtual elements
        original_pos = linspace(1, tran.n_elem, tran.n_elem);

        % Convert phases to radians, unwrap, and interpolate
        unwrapped_phases = unwrap(init_phases * pi / 180) * 180 / pi;
        virtual_pos = linspace(1, tran.n_elem, tran.prestus.transducer.n_elements);
        virtual_phases = interp1(original_pos, unwrapped_phases, virtual_pos, 'linear');

        % Convert back to degrees and apply modulo operation
        init_phases = mod(virtual_phases, 360);
    end

    % Update simulation parameters with the calculated phases
    sim_param.transducer.source_phase_deg = init_phases;

end
