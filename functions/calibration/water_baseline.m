function [freefield_isppa_wcm2, freefield_axial] = water_baseline(parameters, kgrid, source, sensor, source_labels)
% WATER_BASELINE  Measure free-water ISPPA at the current source amplitude
%
% Runs a fast homogeneous water simulation using the current grid (k-Wave:
% cropped to the source and beam, grid.water_baseline_grid), source,
% and sensor, and returns the peak free-water ISPPA [W/cm²] at the focal
% sensor. This value is stored as provenance alongside the cached acoustic
% results so that any target_isppa_wcm2 can later be achieved by post-hoc
% linear pressure scaling, without re-running the full skull simulation.
%
% The simulation backend matches parameters.acoustic_backend so that the
% free-water baseline is produced by the same solver as the tissue simulation.
%   'jwave'  → jwave_wrapper with homogeneous water medium
%   default  → acoustic_simulation (k-Wave)
%
% No amplitude parameters are modified. This function is purely a measurement.
%
% Use as:
%   freefield_isppa_wcm2 = water_baseline(parameters, kgrid, source, sensor)
%   freefield_isppa_wcm2 = water_baseline(parameters, kgrid, source, sensor, source_labels)
%   [freefield_isppa_wcm2, freefield_axial] = water_baseline(...)
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
%   freefield_axial      - [1 x N_transducers] struct with axial_distance_bowl [mm]
%                          and axial_intensity [W/cm²]: the free-water intensity
%                          along each transducer's beam axis, at the current
%                          elem_p_amp (empty fields if the axis could not be sampled)
%
% See also: PRESTUS_PIPELINE, MULTI_ISPPA_PIPELINE, ACOUSTIC_WRAPPER, JWAVE_WRAPPER, SAMPLE_BEAM_AXIS

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
        % grid.water_baseline_grid 'beam' (default): only the box around the
        % source and the beam (WATER_BASELINE_BOX); in water nothing outside
        % it reflects back, so the free-field result is the same
        box_offset = zeros(1, kgrid.dim);
        if ~isfield(parameters.grid, 'water_baseline_grid') || isempty(parameters.grid.water_baseline_grid) || ...
                strcmpi(char(parameters.grid.water_baseline_grid), 'beam')
            [kgrid, source, sensor, box_offset] = water_baseline_box(parameters, kgrid, source, sensor);
        end
        pml_range = [];
        if isfield(parameters.grid, 'pml_size_range'); pml_range = parameters.grid.pml_size_range; end
        input_args = struct( ...
            'PMLInside', false, ...
            'PMLSize',   resolve_pml_size(parameters.grid.pml_size, kgrid, pml_range), ...
            'PlotPML',   false, ...
            'PlotSim',   false);
        if isfield(parameters.grid, 'pml_alpha') && ~isempty(parameters.grid.pml_alpha)
            input_args.PMLAlpha = parameters.grid.pml_alpha;
        end

        fprintf('Running free-water baseline simulation (%s)...\n', bp.simulation.code_type);
        sensor_data_baseline = acoustic_simulation(kgrid, water_medium, source, sensor, ...
                                                   input_args, bp);
    end

    %% Compute ISPPA from peak pressure
    % gather() first: p_max_all is a gpuArray under GPU backends, and double()
    % alone does not move it off the device; a cached gpuArray scalar fails to
    % reload in a later session (see cache_acoustic).
    rho = wp.density;
    c   = wp.sound_speed;
    peak_pa = double(gather(max(sensor_data_baseline.p_max_all(:))));

    if peak_pa <= 0
        error(['water_baseline: water simulation returned zero or negative peak ' ...
               'pressure. Check source setup and sensor coverage of the focal region.']);
    end

    freefield_isppa_wcm2 = (peak_pa^2 / (2 * rho * c)) * 1e-4;
    fprintf('Free-water baseline ISPPA: %.2f W/cm² (peak pressure: %.0f Pa)\n\n', ...
            freefield_isppa_wcm2, peak_pa);

    %% Free-water intensity along each beam axis
    % Distances are measured from the bowl apex, so the profile does not depend
    % on the grid frame (axisymmetric grid or acoustic FOV crop).
    freefield_axial = struct('axial_distance_bowl', {}, 'axial_intensity', {});
    I_water = double(gather(sensor_data_baseline.p_max_all)).^2 / (2 * rho * c) * 1e-4;
    if ~exist('box_offset', 'var'), box_offset = 0; end   % jwave: full grid
    for ti = 1:numel(parameters.transducer)
        freefield_axial(ti).axial_distance_bowl = [];
        freefield_axial(ti).axial_intensity     = [];
        try
            [freefield_axial(ti).axial_distance_bowl, freefield_axial(ti).axial_intensity] = ...
                sample_beam_axis(I_water, parameters.transducer(ti).trans_pos(:)' - box_offset, ...
                                 parameters.transducer(ti).focus_pos(:)' - box_offset, parameters.grid.resolution_mm);
        catch ME
            warning('water_baseline:axialProfile', ...
                'water_baseline: axial profile not extracted for transducer %d (%s).', ti, ME.message);
        end
    end
end


function solver = get_jwave_solver(parameters)
% Return the configured jwave solver string, defaulting to 'helmholtz'.
if isfield(parameters, 'jwave') && isfield(parameters.jwave, 'solver')
    solver = parameters.jwave.solver;
else
    solver = 'helmholtz';
end
end
