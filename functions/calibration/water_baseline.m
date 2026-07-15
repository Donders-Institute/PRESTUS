function freefield_isppa_wcm2 = water_baseline(parameters, kgrid, source, sensor, source_labels)
% WATER_BASELINE  Measure free-water ISPPA at the current source amplitude
%
% Runs a fast homogeneous water simulation using the current grid, source,
% and sensor, and returns the peak free-water ISPPA [W/cm²] at the focal
% sensor. This value is stored as provenance alongside the cached acoustic
% results so that any target_isppa_wcm2 can later be achieved by post-hoc
% linear pressure scaling, without re-running the full skull simulation.
%
% The simulation backend matches parameters.acoustic_backend so that the
% free-water baseline is produced by the same solver as the tissue simulation.
%   'jwave'  → jwave_wrapper with homogeneous water medium
%   default  → acoustic_simulation (k-Wave) — the original behaviour
%
% No amplitude parameters are modified. This function is purely a measurement.
%
% Use as:
%   freefield_isppa_wcm2 = water_baseline(parameters, kgrid, source, sensor)
%   freefield_isppa_wcm2 = water_baseline(parameters, kgrid, source, sensor, source_labels)
%
% Input:
%   parameters    - (1,1) PRESTUS config struct; must contain:
%                     medium_properties.water.*   - sound speed, density, etc.
%                     grid.pml_size               - PML thickness
%                     simulation.code_type        - k-Wave backend
%                     acoustic_backend            - (optional) 'jwave' to use jwave
%   kgrid         - kWaveGrid object from source_sensor_setup
%   source        - k-Wave source struct from source_sensor_setup
%   sensor        - k-Wave sensor struct from source_sensor_setup (unused for jwave)
%   source_labels - (optional) element label grid; passed to jwave_wrapper for provenance
%
% Output:
%   freefield_isppa_wcm2 - peak free-water ISPPA [W/cm²] at current elem_p_amp
%
% See also: PRESTUS_PIPELINE, MULTI_ISPPA_PIPELINE, ACOUSTIC_WRAPPER, JWAVE_WRAPPER

arguments
    parameters    (1,1) struct
    kgrid         (1,1)
    source        (1,1) struct
    sensor        (1,1) struct
    source_labels           = []
end

    fprintf('========================================\n');
    fprintf('FREE-WATER BASELINE MEASUREMENT\n');
    fprintf('========================================\n\n');

    %% Build homogeneous water medium
    wp = parameters.medium_properties.water;
    water_medium = struct( ...
        'sound_speed',          wp.sound_speed, ...
        'density',              wp.density, ...
        'alpha_coeff',          wp.alpha_coeff, ...
        'alpha_power',          wp.alpha_power);

    %% Shared baseline parameters (no interactive display, distinct output affix)
    bp = parameters;
    bp.simulation.interactive = 0;
    bp.simulation.debug       = 0;
    bp.io.output_affix        = [parameters.io.output_affix '_baseline'];

    use_jwave = isfield(parameters, 'acoustic_backend') && ...
                strcmpi(parameters.acoustic_backend, 'jwave');

    if use_jwave
        %% jwave baseline path
        % Use the same jwave solver as the tissue simulation so the free-water
        % ISPPA is produced by an identical forward model.
        bp.simulation.medium = 'water';   % tag the output file as a water run
        fprintf('Running free-water baseline simulation (jwave / %s)...\n', ...
            char(get_jwave_solver(bp)));
        sensor_data_baseline = jwave_wrapper(bp, kgrid, water_medium, source, source_labels);
    else
        %% k-Wave baseline path (default)
        input_args = struct( ...
            'PMLInside', false, ...
            'PMLSize',   resolve_pml_size(parameters.grid.pml_size, kgrid), ...
            'PlotPML',   false, ...
            'PlotSim',   false);

        fprintf('Running free-water baseline simulation (%s)...\n', bp.simulation.code_type);
        sensor_data_baseline = acoustic_simulation(kgrid, water_medium, source, sensor, ...
                                                   input_args, bp);
    end

    %% Compute ISPPA from peak pressure
    rho = wp.density;
    c   = wp.sound_speed;
    peak_pa = double(max(sensor_data_baseline.p_max_all(:)));

    if peak_pa <= 0
        error(['water_baseline: water simulation returned zero or negative peak ' ...
               'pressure. Check source setup and sensor coverage of the focal region.']);
    end

    freefield_isppa_wcm2 = (peak_pa^2 / (2 * rho * c)) * 1e-4;
    fprintf('Free-water baseline ISPPA: %.2f W/cm² (peak pressure: %.0f Pa)\n\n', ...
            freefield_isppa_wcm2, peak_pa);
end


function solver = get_jwave_solver(parameters)
% Return the configured jwave solver string, defaulting to 'helmholtz'.
if isfield(parameters, 'jwave') && isfield(parameters.jwave, 'solver')
    solver = parameters.jwave.solver;
else
    solver = 'helmholtz';
end
end
