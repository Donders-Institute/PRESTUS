function parameters = calibration_amplitude_scaling(parameters, kgrid, source, sensor)
% CALIBRATION_AMPLITUDE_SCALING  Scale source amplitude to a target free-water pressure
%
% Runs a fast homogeneous water simulation using the current grid, source,
% and sensor, measures the peak pressure across all sensor points, and
% linearly scales source_amp_pa so that the free-water focal peak matches
% calibration.target_pressure_mpa.
%
% This stage is inserted into the pipeline after source setup and before the
% main acoustic simulation, so the main sim automatically uses the calibrated
% amplitude. It is only active when modules.run_amplitude_calibration = 1.
%
% Use as:
%   parameters = calibration_amplitude_scaling(parameters, kgrid, source, sensor)
%
% Input:
%   parameters - (1,1) PRESTUS config struct; must contain:
%                  calibration.target_isppa_wcm2    - target free-water ISPPA [W/cm²]
%                  medium_properties.water.*         - sound speed, density, etc.
%                  simulation.code_type             - k-Wave backend
%   kgrid      - kWaveGrid object from source_sensor_setup
%   source     - k-Wave source struct from source_sensor_setup
%   sensor     - k-Wave sensor struct from source_sensor_setup
%
% Output:
%   parameters - updated struct with source_amp_pa scaled to hit target
%                pressure; calibration.measured_pressure_mpa and
%                calibration.amplitude_scale_factor are recorded
%
% See also: PRESTUS_PIPELINE, SOURCE_SENSOR_SETUP, ACOUSTIC_SIMULATION

arguments
    parameters (1,1) struct
    kgrid      (1,1)
    source     (1,1) struct
    sensor     (1,1) struct
end

    %% Validate
    if ~isfield(parameters, 'calibration') || ...
            ~isfield(parameters.calibration, 'target_isppa_wcm2') || ...
            isempty(parameters.calibration.target_isppa_wcm2)
        error(['calibration.target_isppa_wcm2 must be set for amplitude calibration. ' ...
               'Provide a scalar value in W/cm² (e.g. calibration.target_isppa_wcm2: 30).']);
    end

    % Convert target ISPPA [W/cm²] to peak pressure [Pa] via p = sqrt(2 * I * 1e4 * rho * c)
    rho = parameters.medium_properties.water.density;
    c   = parameters.medium_properties.water.sound_speed;
    target_pa = sqrt(2 * parameters.calibration.target_isppa_wcm2 * 1e4 * rho * c);

    fprintf('========================================\n');
    fprintf('AMPLITUDE CALIBRATION (water sim)\n');
    fprintf('Target: %.1f W/cm² ISPPA (%.0f Pa)\n', parameters.calibration.target_isppa_wcm2, target_pa);
    fprintf('========================================\n\n');

    %% Build homogeneous water medium (scalar — k-Wave accepts uniform media)
    wp = parameters.medium_properties.water;
    water_medium = struct( ...
        'sound_speed',        wp.sound_speed, ...
        'density',            wp.density, ...
        'alpha_coeff',        wp.alpha_coeff, ...
        'alpha_power',        wp.alpha_power, ...
        'thermal_conductivity', 0, ...   % stripped inside acoustic_simulation
        'specific_heat',        0, ...   % stripped inside acoustic_simulation
        'perfusion_coeff',      0);      % stripped inside acoustic_simulation

    %% Run calibration water simulation
    % Disable plotting and C++ disk I/O; force matlab_cpu for speed
    cal_parameters = parameters;
    cal_parameters.simulation.code_type  = 'matlab_cpu';
    cal_parameters.simulation.interactive = 0;
    cal_parameters.simulation.debug       = 0;
    % Prevent C++ backends from writing HDF5 files into the output directory
    cal_parameters.io.output_affix = [parameters.io.output_affix '_cal'];

    input_args = struct('PMLInside', false, ...
                        'PMLSize',   parameters.grid.pml_size, ...
                        'PlotPML',   false, ...
                        'PlotSim',   false);

    fprintf('Running calibration water simulation (matlab_cpu)...\n');
    sensor_data_cal = acoustic_simulation(kgrid, water_medium, source, sensor, ...
                                          input_args, cal_parameters);

    %% Measure peak free-water pressure across all sensor points
    % sensor_data.p_max_all is [1 x N_sensor_points] — temporal max at each point
    peak_pa = double(max(sensor_data_cal.p_max_all(:)));
    % Convert measured peak pressure to ISPPA [W/cm²]: I = p² / (2 * rho * c) * 1e-4
    measured_isppa = (peak_pa^2 / (2 * rho * c)) * 1e-4;
    fprintf('Measured free-water ISPPA: %.2f W/cm² (%.0f Pa)\n', measured_isppa, peak_pa);

    if peak_pa <= 0
        error(['Calibration water simulation returned zero or negative peak pressure. ' ...
               'Check source setup and sensor coverage of the focal region.']);
    end

    %% Compute and apply linear amplitude scale factor
    scale = target_pa / peak_pa;
    fprintf('Amplitude scale factor: %.4f (%.1f%%)\n', scale, scale*100);

    % Scale source_amp_pa in every transducer
    for ti = 1:numel(parameters.transducer)
        t = parameters.transducer(ti);
        t_type = t.type;
        if isfield(t.(t_type), 'elem_amp')
            parameters.transducer(ti).(t_type).elem_amp = ...
                parameters.transducer(ti).(t_type).elem_amp * scale;
        end
    end

    % Record calibration results for report / traceability
    parameters.calibration.measured_isppa_wcm2 = measured_isppa;
    parameters.calibration.amplitude_scale_factor = scale;

    fprintf('Calibration complete. Source amplitude scaled by %.4f.\n\n', scale);
end
