function elem_phase_deg = set_real_phases(tran, focal_distance_ep, parameters)
% SET_REAL_PHASES  Compute geometric steering phases for a focal depth
%
% Returns COMPUTE_PHASES for the ring geometry (elem_id_mm, elem_od_mm,
% curv_radius_mm) in the equipment YAML: one phase per ring
% (transducer.annular.elem_n).
%
% Use as:
%   elem_phase_deg = set_real_phases(tran, focal_distance_ep, parameters)
%
% Input:
%   tran              - equipment transducer struct (n_elem, transducer.annular)
%   focal_distance_ep - desired focal depth relative to exit plane [mm]
%   parameters        - PRESTUS config; uses medium_properties.water.sound_speed [m/s]
%
% Output:
%   elem_phase_deg - [1xN] element phases [°]
%
% See also: COMPUTE_PHASES, SAVE_OPTIMIZED_VALUES, CALIBRATION_TRANSDUCER

arguments
    tran              (1,1) struct
    focal_distance_ep (1,1) {mustBeNumeric}
    parameters        (1,1) struct
end

    init_phases = compute_phases(parameters.medium_properties.water.sound_speed, ...
        tran, focal_distance_ep);

    % Adjust phases for virtual elements if necessary.
    % The geometry path returns one phase per annular ring (elem_n), so
    % numel(init_phases) == elem_n and this block is normally skipped.
    if tran.n_elem ~= tran.transducer.annular.elem_n && ...
            numel(init_phases) == tran.n_elem
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
